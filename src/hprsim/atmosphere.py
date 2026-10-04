"""International Standard Atmosphere up to 20 km (troposphere + lower stratosphere)."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

G0 = 9.80665
R_AIR = 287.05287  # J/(kg K)
GAMMA = 1.4
T0, P0, LAPSE = 288.15, 101325.0, 0.0065
H_TROPOPAUSE = 11000.0
T_TROPOPAUSE = T0 - LAPSE * H_TROPOPAUSE
P_TROPOPAUSE = P0 * (T_TROPOPAUSE / T0) ** (G0 / (R_AIR * LAPSE))


@dataclass(frozen=True)
class Air:
    T: float  # K
    P: float  # Pa
    rho: float  # kg/m^3
    a: float  # speed of sound, m/s
    mu: float  # dynamic viscosity, Pa s

    @property
    def nu(self) -> float:
        return self.mu / self.rho


def isa(h: float, dT: float = 0.0) -> Air:
    """ISA at geopotential altitude h (m above sea level). dT shifts temperature (hot/cold day)."""
    h = float(np.clip(h, -500.0, 20000.0))
    if h <= H_TROPOPAUSE:
        T_std = T0 - LAPSE * h
        P = P0 * (T_std / T0) ** (G0 / (R_AIR * LAPSE))
    else:
        T_std = T_TROPOPAUSE
        P = P_TROPOPAUSE * np.exp(-G0 * (h - H_TROPOPAUSE) / (R_AIR * T_std))
    T = T_std + dT
    rho = P / (R_AIR * T)
    mu = 1.458e-6 * T**1.5 / (T + 110.4)  # Sutherland's law
    return Air(T, P, rho, float(np.sqrt(GAMMA * R_AIR * T)), mu)
