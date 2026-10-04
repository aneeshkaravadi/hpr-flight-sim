"""Subsonic drag by component buildup.

Follows the approach in S. Niskanen, "Development of an Open Source model
rocket simulation software" (Helsinki Univ. of Technology, 2009), the method
behind OpenRocket:

  skin friction   Cf from Reynolds number (turbulent, or roughness-limited, whichever is larger),
                  times wetted area with form factors (1 + 1/(2 f_B)) for the body and (1 + 2t/c) for fins
  base drag       0.12 + 0.13 M^2 on the base area (reduced by the motor's area while it is firing)
  fin pressure    rounded leading edges plus blunt trailing edges

Valid for subsonic flight (M < ~0.8). Real rockets have rail buttons, joints and
paint seams this doesn't see, which is why Rocket.cd_scale exists to calibrate
it against a real flight.
"""
from __future__ import annotations

import numpy as np

from .rocket import Rocket


def skin_friction(re: float, mach: float, roughness: float, length: float) -> float:
    if re < 1e4:
        cf = 1.48e-2
    else:
        cf = 1.0 / (1.50 * np.log(re) - 5.6) ** 2
    cf_rough = 0.032 * (roughness / length) ** 0.2
    return max(cf, cf_rough) * (1 - 0.1 * mach**2)


def base_drag_coefficient(mach: float) -> float:
    return 0.12 + 0.13 * mach**2


def drag_coefficient(rocket: Rocket, mach: float, reynolds: float, thrusting: bool) -> float:
    A = rocket.ref_area
    d = rocket.diameter
    f = rocket.fins
    L = rocket.length
    fineness = L / d
    cf = skin_friction(reynolds, mach, rocket.roughness, L)

    body_wet = rocket.nose.wetted_area + rocket.body.wetted_area
    fin_wet = 2 * f.count * f.planform_area
    mean_chord = (f.root_chord + f.tip_chord) / 2
    cd_friction = cf * ((1 + 1 / (2 * fineness)) * body_wet + (1 + 2 * f.thickness / mean_chord) * fin_wet) / A

    cd_b = base_drag_coefficient(mach)
    base_area = A - (np.pi * rocket.motor.diameter**2 / 4 if thrusting else 0.0)
    cd_base = cd_b * base_area / A

    # fins: rounded leading edge (compressibility term) + square trailing edge (acts like a small base)
    edge_area = f.count * f.thickness * f.span
    sweep_angle = np.arctan2(f.sweep, f.span)
    m = min(mach, 0.9)
    cd_le = ((1 - m**2) ** -0.417 - 1) * np.cos(sweep_angle) ** 2 * edge_area / A
    cd_te = cd_b * edge_area / A

    return rocket.cd_scale * (cd_friction + cd_base + cd_le + cd_te)
