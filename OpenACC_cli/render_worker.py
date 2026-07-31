"""Persistent render workers: dlopen libblackhole.so once per worker process
and keep it warm, instead of paying CUDA context init on every single frame
the way spawning ./main per-request did.

Two workers, not one, and each is its own OS process (multiprocessing.Process,
not a thread) so each gets its own CUDA context - deliberately, so a slow
background HQ render can run fully concurrently with interactive drag/settle
renders instead of serializing behind one shared lock. See INTERACTIVE_WORKER
and HQ_WORKER below.
"""

import ctypes
import multiprocessing as mp
import threading
from pathlib import Path

import numpy as np

LIB_PATH = Path("libblackhole.so").resolve()

# Explicit fork context: this project only targets linux-64 (pixi.toml), where
# fork is always available. Neither the parent process nor this module touches
# CUDA before the fork happens, so the classic "CUDA + fork" trap doesn't apply
# here - each worker dlopen's the library and initializes its own CUDA context
# only after it's already a separate process. Explicit > relying on whatever
# multiprocessing's default start method happens to be on a given platform/
# Python version (macOS's "spawn" default, tried during development, can't
# safely start a process from top-level module code at all).
_ctx = mp.get_context("fork")


def _worker_main(request_q, response_q) -> None:
    lib = ctypes.CDLL(str(LIB_PATH))
    lib.render_frame_f32.argtypes = [
        ctypes.c_double, ctypes.c_double, ctypes.c_double,  # r0 theta0 phi0
        ctypes.c_double, ctypes.c_double, ctypes.c_double,  # a Q rs
        ctypes.c_double, ctypes.c_double,                   # errormax de0
        ctypes.c_uint64, ctypes.c_uint64,                   # szeles magas
        ctypes.POINTER(ctypes.c_float),
    ]
    lib.render_frame_f32.restype = ctypes.c_int

    while True:
        job = request_q.get()
        if job is None:  # sentinel: shut down
            return
        job_id, r0, theta0, phi0, a, Q, rs, errormax, de0, szeles, magas = job
        buf = np.empty(3 * szeles * magas, dtype=np.float32)
        ptr = buf.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
        rc = lib.render_frame_f32(r0, theta0, phi0, a, Q, rs, errormax, de0, szeles, magas, ptr)
        response_q.put((job_id, rc, buf if rc == 0 else None))


class RenderWorker:
    """One persistent process holding one warm CUDA context. render() is
    thread-safe (calls from multiple Python threads serialize onto this one
    worker via _lock), but does NOT run concurrently with other calls to the
    same worker - use a second RenderWorker for that."""

    def __init__(self) -> None:
        self._request_q = _ctx.Queue()
        self._response_q = _ctx.Queue()
        self._process = _ctx.Process(target=_worker_main, args=(self._request_q, self._response_q), daemon=True)
        self._process.start()
        self._lock = threading.Lock()
        self._next_id = 0

    def render(self, r0: float, theta0: float, phi0: float, a: float, Q: float, rs: float,
               errormax: float, de0: float, szeles: int, magas: int) -> np.ndarray:
        with self._lock:
            self._next_id += 1
            job_id = self._next_id
            self._request_q.put((job_id, r0, theta0, phi0, a, Q, rs, errormax, de0, szeles, magas))
            got_id, rc, buf = self._response_q.get()
            assert got_id == job_id  # this worker only ever serves one caller at a time

        if rc != 0:
            raise RuntimeError(f"render_frame_f32 failed with code {rc} (likely a resolution aspect-ratio mismatch)")
        return buf.reshape(szeles, magas, 3)

    def shutdown(self) -> None:
        # daemon=True alone isn't reliable here: on a graceful SIGTERM,
        # uvicorn's own shutdown path doesn't reliably run Python's
        # atexit-based multiprocessing daemon cleanup (verified directly -
        # the worker processes were still alive afterward), so without this
        # explicit call every restart leaks a GPU-context-holding process.
        self._request_q.put(None)
        self._process.join(timeout=5)
        if self._process.is_alive():
            self._process.terminate()


# Separate processes/CUDA contexts so a slow HQ render never blocks interactive frames.
INTERACTIVE_WORKER = RenderWorker()
HQ_WORKER = RenderWorker()


def shutdown_all() -> None:
    INTERACTIVE_WORKER.shutdown()
    HQ_WORKER.shutdown()
