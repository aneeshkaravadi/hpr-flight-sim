"""Calibrate the drag model against a real flight.

Every drag model misses something (rail buttons, joints, paint, a slightly
crooked fin). The standard fix is to scale the drag coefficient until the
simulated apogee matches the altimeter, then trust the calibrated rocket for
the next flight on a different motor. That cross-prediction is the honest test.
"""
from __future__ import annotations

import csv
from collections.abc import Callable
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from .flight import Launch, simulate
from .rocket import Rocket

FT = 0.3048


def read_altitude_log(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """Time (s) and altitude above ground (m) from an altimeter CSV or an OpenRocket CSV export.

    Finds the first column whose header contains 'time' and the first containing 'alt'.
    A header mentioning 'ft' or 'feet' is converted to meters.
    """
    rows = Path(path).read_text().splitlines()
    header_idx = next(i for i, r in enumerate(rows) if "time" in r.lower() and "alt" in r.lower())
    header = [h.strip().lstrip("#").strip() for h in rows[header_idx].split(",")]
    ti = next(i for i, h in enumerate(header) if "time" in h.lower())
    ai = next(i for i, h in enumerate(header) if "alt" in h.lower())
    scale = FT if any(u in header[ai].lower() for u in ("ft", "feet")) else 1.0
    t, a = [], []
    for r in csv.reader(rows[header_idx + 1:]):
        if not r or r[0].strip().startswith("#"):
            continue
        try:
            t.append(float(r[ti]))
            a.append(float(r[ai]) * scale)
        except (ValueError, IndexError):
            continue
    return np.array(t), np.array(a)


def fit_cd_scale(make_rocket: Callable[[float], Rocket], launch: Launch, measured_apogee: float,
                 lo: float = 0.4, hi: float = 2.5) -> float:
    """Drag multiplier that makes the simulated apogee equal the measured one."""
    def err(scale):
        return simulate(make_rocket(scale), launch, ascent_only=True).apogee - measured_apogee
    return brentq(err, lo, hi, xtol=1e-4)
