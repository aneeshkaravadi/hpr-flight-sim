"""Rocket motors from RASP .eng thrust curves (the format thrustcurve.org publishes).

Propellant is assumed to burn in proportion to the impulse delivered, so the
propellant left at time t is  m_p0 * (1 - I(t) / I_total).  That is the
standard assumption in OpenRocket and RockSim, and it makes the mass flow
m_dot = m_p0 * F(t) / I_total, i.e. a constant effective exhaust velocity
c = I_total / m_p0.

Thrust at altitude: the curve is measured on a static test stand near sea
level. Higher up, the nozzle exhausts into lower ambient pressure and gains
pressure thrust,  F = F_curve + (p_ref - p_ambient) * A_exit,  while the mass
flow is set by the chamber and doesn't change. A .eng file doesn't include the
nozzle exit area, so it is an input (exit_area = 0 turns the correction off).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


@dataclass
class Motor:
    name: str
    diameter: float  # m
    length: float  # m
    prop_mass: float  # kg
    total_mass: float  # kg (loaded)
    t: np.ndarray  # s
    F: np.ndarray  # N
    impulse_scale: float = 1.0  # for Monte Carlo: scales thrust (and burn rate) uniformly
    exit_area: float = 0.0  # m^2, nozzle exit; 0 = no ambient-pressure correction
    p_ref: float = 101325.0  # Pa, ambient pressure the thrust curve was measured at
    _I: np.ndarray = field(init=False, repr=False)

    def __post_init__(self):
        self.t = np.asarray(self.t, float)
        self.F = np.asarray(self.F, float)
        if self.t[0] > 0:  # most .eng files start at the first nonzero point
            self.t = np.r_[0.0, self.t]
            self.F = np.r_[0.0, self.F]
        dt = np.diff(self.t)
        self._I = np.r_[0.0, np.cumsum(0.5 * (self.F[1:] + self.F[:-1]) * dt)]

    @classmethod
    def from_eng(cls, path: str | Path) -> Motor:
        lines = [ln.strip() for ln in Path(path).read_text().splitlines()]
        lines = [ln for ln in lines if ln and not ln.startswith(";")]
        name, dia_mm, len_mm, _delays, prop_kg, total_kg, *_maker = lines[0].split()
        pts = np.array([[float(v) for v in ln.split()[:2]] for ln in lines[1:] if len(ln.split()) >= 2])
        return cls(name, float(dia_mm) / 1000, float(len_mm) / 1000, float(prop_kg), float(total_kg),
                   pts[:, 0], pts[:, 1])

    @property
    def burn_time(self) -> float:
        return float(self.t[-1])

    @property
    def total_impulse(self) -> float:
        return float(self._I[-1]) * self.impulse_scale

    @property
    def case_mass(self) -> float:
        return self.total_mass - self.prop_mass

    @property
    def exhaust_velocity(self) -> float:
        return self.total_impulse / self.prop_mass

    def thrust(self, t: float, p_ambient: float | None = None) -> float:
        """Thrust at time t. Pass the ambient pressure to include the altitude correction."""
        if t < 0 or t >= self.t[-1]:
            return 0.0
        F = float(np.interp(t, self.t, self.F)) * self.impulse_scale
        if p_ambient is not None and self.exit_area > 0 and F > 0:
            F += (self.p_ref - p_ambient) * self.exit_area
        return F

    def prop_remaining(self, t: float) -> float:
        if t <= 0:
            return self.prop_mass
        if t >= self.t[-1]:
            return 0.0
        return self.prop_mass * (1.0 - float(np.interp(t, self.t, self._I)) / float(self._I[-1]))

    def mass_flow(self, t: float) -> float:
        """From the measured curve only: pressure thrust doesn't burn extra propellant."""
        return self.thrust(t) / self.exhaust_velocity
