"""Regression snapshot of the entry points, for checking that a change is output-identical.

`run` executes, and hashes in a SHA-256 manifest:

- every registered model and `_ext` variant on every fatigue-life workbook, through
  `scripts/1_run_material.py` (workbook constants only), then the calibration and the comparison;
- `multifat-pipeline --models all` on every workbook (#73): calibrated plane constants (FS, FIN,
  LI, MGSE, MSWT), the fatigue-index path of the fatigue-limit workbooks, the sidecars and run
  reports; and one cross-batch transfer (42CrMo4_Reis designed on 42CrMo4_Santus).

Volatile fields of the pipeline outputs (dates, calibration ids, git state, timings) are left out
of their hashes, and the archived calibration copies, named by id, are not hashed.
`compare` diffs two snapshots; logs are compared after masking paths that legitimately
differ between runs (output directory, script and package locations, traceback lines).

    uv run python tools/regression_snapshot.py run OUT [--jobs 4] [--skip-calibrate]
        [--skip-pipeline] [--data-folder DIR]
    uv run python tools/regression_snapshot.py compare BASE NEW

Run from the repository root. Outputs go to OUT only; tracked Results are untouched.
The script stage takes about 20 minutes on a laptop and the pipeline stage about as long again
(`--skip-pipeline` leaves it out). `--data-folder` runs every workbook of another folder (e.g.
frozen inputs). `compare` checks CSVs on the base snapshot's columns, so columns appended by a
later revision are not differences. The materials come from the workbooks and their
classification, so adding a workbook needs no edit here.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

DATA_FOLDER = ROOT / "Experimental_Data"
# One cross-batch transfer of the pipeline stage: target and design workbooks (issue #71)
TRANSFER = ("14_42CrMo4_Reis.xlsx", "6_42CrMo4_Santus.xlsx")
# Pipeline outputs whose content changes from run to run without a numerical change
VOLATILE_KEYS = {"created_utc", "created_at", "calibration_id", "git", "seconds"}


def workbooks(data_folder: Path | None) -> list[Path]:
    from multifat.imp.schema import workbook_files

    return list(workbook_files(data_folder or DATA_FOLDER))


def life_workbooks(data_folder: Path | None) -> list[Path]:
    """Workbooks classified ``fatigue_life``: the ones the script stage can run."""
    from multifat.imp.classification import classify_dataset
    from multifat.imp.experimental_data import load_material_workbook

    out = []
    for path in workbooks(data_folder):
        try:
            ((material, block),) = load_material_workbook(path).items()
        except Exception:  # noqa: BLE001 - an unreadable workbook is left to the script stage
            out.append(path)
            continue
        if classify_dataset(block, material)["workflow"] == "fatigue_life":
            out.append(path)
    return out


def material_args(data_folder: Path | None) -> list[str]:
    """`--material` values of the script stage: the fatigue-life workbooks (paths)."""
    return [str(p) for p in life_workbooks(data_folder)]


def _without_volatile(value):
    if isinstance(value, dict):
        return {k: _without_volatile(v) for k, v in value.items() if k not in VOLATILE_KEYS}
    if isinstance(value, list):
        return [_without_volatile(v) for v in value]
    return value


def stable_bytes(out: Path, path: Path) -> bytes | None:
    """Content hashed for the snapshot file `path`; None for a file that is not hashed.

    Pipeline sidecars (.json) and run reports (.tsv) lose their volatile fields; the copies
    archived under ``calibrations/``, named by calibration id, are not hashed.
    """
    parts = path.relative_to(out).parts
    if parts[0] == "pipeline" and "calibrations" in parts:
        return None
    if parts[0] == "pipeline" and path.suffix == ".json":
        payload = _without_volatile(json.loads(path.read_text(encoding="utf-8")))
        return json.dumps(payload, sort_keys=True).encode()
    if parts[0] == "pipeline" and path.suffix == ".tsv":
        rows = list(csv.reader(path.read_text(encoding="utf-8").splitlines(), delimiter="\t"))
        keep = [i for i, h in enumerate(rows[0]) if h not in VOLATILE_KEYS] if rows else []
        return "\n".join("\t".join(r[i] for i in keep if i < len(r)) for r in rows).encode()
    return path.read_bytes()


def _models() -> list[str]:
    from multifat.models.registry import model_names

    bases = model_names()
    return bases + [b + "_ext" for b in bases]


def run(
    out: Path,
    jobs: int,
    skip_calibrate: bool,
    data_folder: Path | None = None,
    skip_pipeline: bool = False,
) -> None:
    out.mkdir(parents=True, exist_ok=True)
    (out / "logs").mkdir(exist_ok=True)
    env = dict(os.environ, MPLBACKEND="Agg", PYTHONHASHSEED="0")

    def call(args: list[str], log: str) -> int:
        p = subprocess.run(
            [sys.executable, *args], capture_output=True, text=True, env=env, cwd=ROOT
        )
        text = p.stdout + "\n--- stderr ---\n" + p.stderr
        (out / "logs" / log).write_text(text, encoding="utf-8")
        return p.returncode

    def one(job: tuple[str, str]) -> tuple[str, int]:
        model, mat = job
        args = ["scripts/1_run_material.py", "--model", model, "--material", mat]
        args += ["--results-dir", str(out / "runs"), "--no-plots"]
        name = Path(mat).stem
        return f"{model}|{name}", call(args, f"run_{model}_{name}.log")

    mats = material_args(data_folder)
    with ThreadPoolExecutor(jobs) as ex:
        codes = dict(ex.map(one, [(m, mat) for m in _models() for mat in mats]))
    if not skip_calibrate:
        cal = ["scripts/0_calibrate_material_params.py"]
        if data_folder is not None:
            cal += ["--data-folder", str(data_folder)]
        codes["calibrate"] = call(cal, "calibrate.log")
    compare_args = ["scripts/2_compare_results.py", "--results_dir", str(out / "runs")]
    codes["compare"] = call(compare_args + ["--out_dir", str(out / "comparison")], "compare.log")
    if not skip_pipeline:
        codes.update(pipeline_stage(out, call, jobs, data_folder))

    files = {}
    for f in sorted(out.rglob("*")):
        if f.is_file() and f.name != "manifest.json":
            data = stable_bytes(out, f)
            if data is not None:
                files[f.relative_to(out).as_posix()] = hashlib.sha256(data).hexdigest()
    manifest = {"codes": codes, "files": files}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1), encoding="utf-8")
    print("nonzero exit codes:", {k: v for k, v in codes.items() if v})
    print(len(files), "files")


def pipeline_stage(out: Path, call, jobs: int, data_folder: Path | None) -> dict:
    """`multifat-pipeline --models all` on every workbook, and one transfer; exit codes."""
    base = ["-m", "multifat.cli.pipeline", "--models", "all", "--no-compare", "--jobs", "1"]

    def one(path: Path) -> tuple[str, int]:
        dest = out / "pipeline" / path.stem
        args = [*base, str(path), "--results-dir", str(dest)]
        return f"pipeline|{path.stem}", call(args, f"pipeline_{path.stem}.log")

    paths = workbooks(data_folder)
    with ThreadPoolExecutor(jobs) as ex:
        codes = dict(ex.map(one, paths))
    folder = data_folder or DATA_FOLDER
    target, design = (folder / name for name in TRANSFER)
    if target.exists() and design.exists():
        dest = out / "pipeline" / f"{target.stem}_from_{design.stem}"
        args = [*base, str(target), "--design-material", design.stem, "--results-dir", str(dest)]
        codes[f"pipeline|{dest.name}"] = call(args, f"pipeline_{dest.name}.log")
    return codes


_MASKS = [
    (re.compile(r"[^\s'\"{}]*?[\\/]comparison[\\/]"), "<OUT>/comparison/"),
    (re.compile(r"[^\s'\"]*?[\\/]runs[\\/]?"), "<OUT>/runs/"),
    (re.compile(r'File "[^"]*"'), 'File "<PATH>"'),
    (re.compile(r"[^\s'\"]*?[\\/](src[\\/])?(multifat|fatigue)[\\/]"), "<PKG>/"),
    (re.compile(r"line \d+"), "line N"),
    (re.compile(r"\.py:\d+:"), ".py:N:"),
]


def _norm(text: str, root: Path | None = None) -> str:
    if root is not None:
        for form in {str(root), root.as_posix(), str(root).replace("\\", "\\\\")}:
            text = text.replace(form, "<OUT>")
    for rx, rep in _MASKS:
        text = rx.sub(rep, text)
    lines = text.splitlines(True)
    out = []
    in_tb = False
    for line in lines:
        if "Traceback (most recent call last):" in line:
            in_tb = True
            out.append(line)
        elif in_tb:
            if line.startswith("  File "):
                out.append(line)
            elif line.startswith("    ") or set(line.strip()) <= {"^", "~"}:
                continue
            else:
                in_tb = False
                out.append(line)
        else:
            out.append(line)
    return "".join(out)


def _same_on_base_columns(a: Path, b: Path) -> bool:
    """Same rows, and every base column unchanged; columns appended in `b` are ignored."""
    with open(a, newline="", encoding="utf-8") as fa, open(b, newline="", encoding="utf-8") as fb:
        ra, rb = list(csv.reader(fa)), list(csv.reader(fb))
    if len(ra) != len(rb) or not ra:
        return False
    n = len(ra[0])
    return all(x == y[:n] for x, y in zip(ra, rb))


def compare(base: Path, new: Path) -> int:
    a = json.loads((base / "manifest.json").read_text(encoding="utf-8"))
    b = json.loads((new / "manifest.json").read_text(encoding="utf-8"))
    status = 0
    if a["codes"] != b["codes"]:
        status = 1
        for k in sorted(set(a["codes"]) | set(b["codes"])):
            if a["codes"].get(k) != b["codes"].get(k):
                print("exit code differs:", k, a["codes"].get(k), b["codes"].get(k))
    fa, fb = a["files"], b["files"]
    outputs = sorted(k for k in set(fa) | set(fb) if not k.startswith("logs/"))
    changed = []
    for k in outputs:
        if fa.get(k) == fb.get(k):
            continue
        if k in fa and k in fb and k.endswith(".tex"):
            ta = _norm((base / k).read_text(encoding="utf-8"), base)
            tb = _norm((new / k).read_text(encoding="utf-8"), new)
            if ta == tb:
                continue
        if k in fa and k in fb and k.endswith(".csv"):
            if _same_on_base_columns(base / k, new / k):
                continue
        changed.append(k)
    logs = [k for k in fa if k.startswith("logs/") and k in fb]
    log_changed = [
        k
        for k in logs
        if _norm((base / k).read_text(encoding="utf-8"), base)
        != _norm((new / k).read_text(encoding="utf-8"), new)
    ]
    print(f"outputs: {len(outputs)} files, {len(changed)} differ (after masking .tex paths)")
    for k in changed[:40]:
        print("  ", k)
    print(f"logs: {len(logs)} files, {len(log_changed)} differ after masking")
    for k in log_changed[:20]:
        print("  ", k)
    return 1 if (changed or log_changed or status) else 0


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("out", type=Path)
    r.add_argument("--jobs", type=int, default=4)
    r.add_argument("--skip-calibrate", action="store_true")
    r.add_argument("--skip-pipeline", action="store_true", help="leave out the pipeline stage")
    r.add_argument("--data-folder", type=Path, default=None)
    c = sub.add_parser("compare")
    c.add_argument("base", type=Path)
    c.add_argument("new", type=Path)
    return ap


def main() -> None:
    args = build_parser().parse_args()
    if args.cmd == "run":
        folder = args.data_folder.resolve() if args.data_folder else None
        run(args.out.resolve(), args.jobs, args.skip_calibrate, folder, args.skip_pipeline)
    else:
        sys.exit(compare(args.base.resolve(), args.new.resolve()))


if __name__ == "__main__":
    main()
