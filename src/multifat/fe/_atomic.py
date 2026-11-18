"""Protect solver sources and existing destinations while publishing converted files."""

from contextlib import contextmanager
from pathlib import Path
import tempfile


@contextmanager
def atomic_destination(out_path, sources):
    """Yield a temporary sibling and replace the destination only on success."""
    out_path = Path(out_path)
    for source in map(Path, sources):
        if out_path.resolve() == source.resolve() or (
            out_path.exists() and source.exists() and out_path.samefile(source)
        ):
            raise ValueError("output destination must differ from every source file")
    handle = tempfile.NamedTemporaryFile(dir=out_path.parent, suffix=out_path.suffix, delete=False)
    temporary = Path(handle.name)
    handle.close()
    try:
        yield temporary
        temporary.replace(out_path)
    finally:
        temporary.unlink(missing_ok=True)
