"""Run many simulations in parallel with a live progress bar you can speed up or slow down.

In a terminal:   +  more workers    -  fewer workers    p  pause / resume    q  stop early

Workers run at the lowest CPU priority (nice 19), and only as many simulations
as you've asked for are ever in flight, so you can watch your CPU temperature
and dial it up or down while it runs. Without a terminal (e.g. a CI log) it
prints a progress line every few seconds instead.
"""
from __future__ import annotations

import os
import select
import sys
import termios
import time
import tty
from collections import deque
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from contextlib import contextmanager
from typing import Any, Callable, Iterable


def _low_priority():
    try:
        os.nice(19)
    except OSError:
        pass


@contextmanager
def _keyboard(enabled: bool):
    """Single-key reads without Enter; always restores the terminal, even on Ctrl-C."""
    if not enabled:
        yield lambda: ""
        return
    fd = sys.stdin.fileno()
    old = termios.tcgetattr(fd)
    try:
        tty.setcbreak(fd)

        def read():
            # os.read on the raw fd: sys.stdin's buffer can hide keypresses from select()
            keys = ""
            while select.select([fd], [], [], 0)[0]:
                keys += os.read(fd, 64).decode(errors="ignore")
            return keys

        yield read
    finally:
        termios.tcsetattr(fd, termios.TCSADRAIN, old)


def _fmt(seconds: float) -> str:
    if seconds != seconds or seconds == float("inf"):
        return "--"
    seconds = int(round(seconds))
    return f"{seconds // 60}m{seconds % 60:02d}s" if seconds >= 60 else f"{seconds}s"


def run(fn: Callable[[Any], Any], tasks: Iterable[Any], workers: int = 2, label: str = "simulations",
        max_workers: int | None = None) -> list:
    """Map ``fn`` over ``tasks`` in worker processes. Returns results in task order (None if stopped early).

    ``fn`` must be a top-level (picklable) function.
    """
    tasks = list(tasks)
    n = len(tasks)
    max_workers = max_workers or os.cpu_count() or 2
    workers = max(1, min(workers, max_workers))
    results: list = [None] * n
    queue = list(range(n))[::-1]
    in_flight: dict = {}
    durations: deque = deque(maxlen=40)
    done, paused, stop = 0, False, False
    t_start = last_log = time.time()
    interactive = sys.stdin.isatty() and sys.stdout.isatty()

    def draw(final=False):
        nonlocal last_log
        remaining = n - done
        per_task = sum(durations) / len(durations) if durations else float("nan")
        eta = remaining * per_task / workers if not paused else float("inf")
        pct = done / n if n else 1.0
        width = 30
        bar = "#" * int(width * pct) + "-" * (width - int(width * pct))
        state = "PAUSED" if paused else f"workers {workers}/{max_workers}"
        line = (f"{label}: [{bar}] {done}/{n} {100 * pct:5.1f}%  {state}  "
                f"elapsed {_fmt(time.time() - t_start)}  left ~{_fmt(eta)}")
        if interactive:
            sys.stdout.write("\r\x1b[2K" + line + ("" if final else "   (+ faster, - slower, p pause, q quit)"))
            if final:
                sys.stdout.write("\n")
            sys.stdout.flush()
        elif final or time.time() - last_log > 5:
            print(line, flush=True)
            last_log = time.time()

    with _keyboard(interactive) as read_keys, ProcessPoolExecutor(max_workers=max_workers, initializer=_low_priority) as ex:
        try:
            while done < n and not stop:
                for k in read_keys():
                    if k in "+=":
                        workers = min(workers + 1, max_workers)
                    elif k in "-_":
                        workers = max(workers - 1, 1)
                    elif k in "pP":
                        paused = not paused
                    elif k in "qQ":
                        stop = True
                target = 0 if paused else workers
                while len(in_flight) < target and queue:
                    i = queue.pop()
                    in_flight[ex.submit(fn, tasks[i])] = (i, time.time())
                if in_flight:
                    finished, _ = wait(list(in_flight), timeout=0.25, return_when=FIRST_COMPLETED)
                    for fut in finished:
                        i, t_sub = in_flight.pop(fut)
                        results[i] = fut.result()
                        durations.append(time.time() - t_sub)
                        done += 1
                else:
                    time.sleep(0.25)
                draw()
        except KeyboardInterrupt:
            stop = True
        if stop:
            for fut in in_flight:
                fut.cancel()
            ex.shutdown(wait=True, cancel_futures=True)
    draw(final=True)
    if stop:
        print(f"stopped early: {done} of {n} {label} finished", flush=True)
    return results
