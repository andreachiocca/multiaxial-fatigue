"""Check scientific dossier files and benchmark JSON without dependencies."""

import argparse
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1] / "docs/scientific-reviews"


def build_parser() -> argparse.ArgumentParser:
    return argparse.ArgumentParser(description=__doc__.splitlines()[0])


def main():
    build_parser().parse_args()
    reviews = [path for path in sorted(ROOT.iterdir()) if path.is_dir()]
    if not reviews:
        raise SystemExit("No scientific review dossiers found")
    for path in reviews:
        required = {
            "README.md",
            f"{path.name}.tex",
            "references.bib",
            "benchmarks.json",
            f"{path.name}.pdf",
        }
        missing = required - {item.name for item in path.iterdir()}
        if missing:
            raise SystemExit(f"{path.name}: missing {', '.join(sorted(missing))}")
        try:
            data = json.loads((path / "benchmarks.json").read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            raise SystemExit(f"{path.name}: invalid benchmarks.json: {exc}") from exc
        cases = data.get("benchmarks")
        if data.get("method_id") != path.name or not isinstance(cases, list) or not cases:
            raise SystemExit(f"{path.name}: method ID or benchmark list is invalid")
        ids = [case.get("id") for case in cases if isinstance(case, dict)]
        if (
            len(ids) != len(cases)
            or any(not isinstance(case_id, str) or not case_id for case_id in ids)
            or len(ids) != len(set(ids))
        ):
            raise SystemExit(f"{path.name}: benchmark IDs must be nonempty and unique")
    print(f"{len(reviews)} scientific review dossier(s); required files and benchmark JSON OK")


if __name__ == "__main__":
    main()
