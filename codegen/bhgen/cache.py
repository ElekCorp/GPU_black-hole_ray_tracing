"""Disk cache for the expensive symbolic derivations.

Deriving and simplifying the Kerr-Newman geodesic right-hand side takes minutes,
which is fine once but painful when iterating on the emitter or the event
declarations.  The cache key is a hash of the metric module's own source, so
editing the metric invalidates it automatically and editing anything else does
not.
"""

from __future__ import annotations

import hashlib
import inspect
import pathlib
import pickle
import sys

CACHE_DIR = pathlib.Path(__file__).resolve().parent.parent / ".cache"


def _key(spec, tag: str) -> pathlib.Path:
    mod = sys.modules.get(type(spec).__module__)
    try:
        src = inspect.getsource(sys.modules[spec.__module__])
    except Exception:
        src = ""
    if not src:
        for name, m in list(sys.modules.items()):
            if name.endswith(f"metrics.{spec.name}"):
                src = inspect.getsource(m)
                break
    h = hashlib.sha256(src.encode()).hexdigest()[:16]
    return CACHE_DIR / f"{spec.name}.{tag}.{h}.pkl"


def cached(spec, tag: str, produce):
    """Return (value, was_a_cache_hit)."""
    path = _key(spec, tag)
    if path.exists():
        try:
            return pickle.loads(path.read_bytes()), True
        except Exception:
            pass
    value = produce()
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    try:
        path.write_bytes(pickle.dumps(value))
    except Exception:
        pass
    return value, False
