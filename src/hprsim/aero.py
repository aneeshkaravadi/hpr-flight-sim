"""Drag by component buildup, from subsonic speeds through Mach ~1.2.

Follows the approach in S. Niskanen, "Development of an Open Source model
rocket simulation software" (Helsinki Univ. of Technology, 2009), the method
behind OpenRocket:

  skin friction   Cf from Reynolds number (turbulent, or roughness-limited, whichever is larger),
                  times wetted area with form factors (1 + 1/(2 f_B)) for the body and (1 + 2t/c) for fins;
                  compressibility (1 - 0.1 M^2) below Mach 1, and the supersonic corrections above it
  base drag       0.12 + 0.13 M^2 below Mach 1 and 0.25 / M above, on the base area
                  (less the motor's own area while it is firing, because the exhaust fills it)
  fin edges       by cross-section: a square leading edge feels stagnation pressure, a rounded one
                  much less; a square trailing edge acts like a small base, a rounded one like half
                  of one (my assumption), and an airfoil's sharp trailing edge like none

Not included: the nose cone's own pressure (wave) drag, which is small below
Mach ~0.8 for smooth noses but grows through Mach 1. Compared with OpenRocket
it is the missing piece above Mach 0.9 (see the README). Real rockets also have
rail buttons, joints and paint seams this doesn't see, which is why
Rocket.cd_scale exists to calibrate it against a real flight.
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
    if mach < 1.0:
        return max(cf, cf_rough) * (1 - 0.1 * mach**2)
    return max(cf / (1 + 0.15 * mach**2) ** 0.58, cf_rough / (1 + 0.18 * mach**2))


def base_drag_coefficient(mach: float) -> float:
    return 0.12 + 0.13 * mach**2 if mach < 1.0 else 0.25 / mach


def stagnation_drag_coefficient(mach: float) -> float:
    """Blunt face square to the flow: 0.85 times the stagnation pressure ratio q_stag / q."""
    if mach <= 1.0:
        ratio = 1 + mach**2 / 4 + mach**4 / 40
    else:
        ratio = 1.84 - 0.76 / mach**2 + 0.166 / mach**4 + 0.035 / mach**6
    return 0.85 * ratio


def rounded_edge_drag_coefficient(mach: float) -> float:
    """Rounded leading edge; the three pieces meet at Mach 0.9 and 1."""
    if mach < 0.9:
        return (1 - mach**2) ** -0.417 - 1
    if mach < 1.0:
        return 1 - 1.785 * (mach - 0.9)
    return 1.214 - 0.502 / mach**2 + 0.1095 / mach**4


TRAILING_EDGE = {"square": 1.0, "rounded": 0.5, "airfoil": 0.0}  # share of base drag on the trailing edge


def drag_components(rocket: Rocket, mach: float, reynolds: float, thrusting: bool) -> dict:
    """Each drag coefficient (on A_ref) before the calibration multiplier: friction, base and fin pressure."""
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

    # fin edges: the leading edge by its shape (scaled by the sweep), the trailing edge as a share of base drag
    edge_area = f.count * f.thickness * f.span
    sweep_angle = np.arctan2(f.sweep, f.span)
    le = stagnation_drag_coefficient(mach) if f.cross_section == "square" else rounded_edge_drag_coefficient(mach)
    cd_le = le * np.cos(sweep_angle) ** 2 * edge_area / A
    cd_te = TRAILING_EDGE[f.cross_section] * cd_b * edge_area / A
    return {"friction": cd_friction, "base": cd_base, "fin_edges": cd_le + cd_te}


def drag_coefficient(rocket: Rocket, mach: float, reynolds: float, thrusting: bool) -> float:
    return rocket.cd_scale * sum(drag_components(rocket, mach, reynolds, thrusting).values())
