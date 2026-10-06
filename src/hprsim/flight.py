"""6-DOF ascent, 3-DOF descent.

World frame: x east, y north, z up, origin at the launch rail (ground level).
Body frame: x_b along the rocket axis toward the nose.

Ascent (rigid body, quaternion attitude):
  m dv/dt  = thrust * x_b + gravity + drag (opposite the airspeed) + normal force at the CP
  I dw/dt  = M - w x (I w)
  normal force   N = q A CN_alpha * alpha, acting at the CP, pushing the body toward the crossflow
  pitch damping  aerodynamic  0.5 rho V A sum CN_alpha,i (x_i - x_cg)^2
                 jet          m_dot (x_nozzle - x_cg)^2        (Barrowman / Mandell)
  Roll is not modeled (no fin cant).

On the rail the rocket only slides along the rail. After apogee the recovery
system turns it into a point mass hanging under a drag area (Cd * A). A canopy
doesn't open instantly: its drag area grows with the square of the distance it
has traveled since deployment, reaching full size after fill_constant canopy
diameters, and the peak force while it fills is the opening shock.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .aero import drag_coefficient
from .atmosphere import G0, isa
from .rocket import Rocket


@dataclass
class Launch:
    rail_length: float = 1.83  # m of rail the rocket travels before it is free (6 ft rail)
    rail_tilt_deg: float = 0.0  # from vertical
    rail_azimuth_deg: float = 0.0  # compass direction the rail leans toward (0 = north, 90 = east)
    site_elevation: float = 200.0  # m above sea level
    wind_speed: float = 0.0  # m/s at 10 m above ground, 1/7 power-law profile
    wind_from_deg: float = 270.0  # where the wind comes from (meteorological convention)
    dT: float = 0.0  # temperature offset from ISA, K
    gravity: float = G0

    def wind(self, z: float) -> np.ndarray:
        if self.wind_speed == 0:
            return np.zeros(3)
        s = self.wind_speed * (max(z, 1.0) / 10.0) ** (1 / 7)
        toward = np.radians(self.wind_from_deg + 180.0)
        return np.array([s * np.sin(toward), s * np.cos(toward), 0.0])

    @property
    def rail_direction(self) -> np.ndarray:
        e, a = np.radians(self.rail_tilt_deg), np.radians(self.rail_azimuth_deg)
        return np.array([np.sin(e) * np.sin(a), np.sin(e) * np.cos(a), np.cos(e)])


@dataclass
class Recovery:
    """Dual deploy by default. drogue_cda=None or main_altitude=None means the main opens at apogee."""

    drogue_cda: float | None = 0.057  # m^2 (12 in drogue); None = no drogue
    main_cda: float = 1.41  # m^2 (60 in main)
    main_altitude: float | None = 150.0  # m AGL; None = main opens at apogee
    drogue_diameter: float = 0.305  # m, nominal (sets how far it travels while filling)
    main_diameter: float = 1.524  # m
    fill_constant: float = 8.0  # canopy diameters traveled while filling (assumed; varies by canopy type)


@dataclass
class Flight:
    t: np.ndarray
    pos: np.ndarray  # (n, 3)
    vel: np.ndarray  # (n, 3)
    tilt_deg: np.ndarray  # body axis from vertical (NaN in descent)
    aoa_deg: np.ndarray
    mach: np.ndarray
    margin: np.ndarray  # static margin, calibers (NaN in descent)
    events: dict = field(default_factory=dict)

    @property
    def apogee(self) -> float:
        return self.events["apogee_m"]


# ---------------------------------------------------------------- quaternions (w, x, y, z)

def quat_mul(a, b):
    aw, ax, ay, az = a
    bw, bx, by, bz = b
    return np.array([aw * bw - ax * bx - ay * by - az * bz,
                     aw * bx + ax * bw + ay * bz - az * by,
                     aw * by - ax * bz + ay * bw + az * bx,
                     aw * bz + ax * by - ay * bx + az * bw])


def quat_to_mat(q):
    w, x, y, z = q
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
                     [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
                     [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)]])


def quat_from_to(a, b):
    """Shortest rotation taking unit vector a to unit vector b."""
    a, b = a / np.linalg.norm(a), b / np.linalg.norm(b)
    c = np.cross(a, b)
    d = float(np.dot(a, b))
    if d < -0.999999:
        axis = np.cross(a, [0, 1, 0]) if abs(a[0]) > 0.9 else np.cross(a, [1, 0, 0])
        axis /= np.linalg.norm(axis)
        return np.array([0.0, *axis])
    q = np.array([1 + d, *c])
    return q / np.linalg.norm(q)


# ---------------------------------------------------------------- dynamics

def _aero(rocket: Rocket, launch: Launch, t, r, v, Rm, xcg, air=None):
    """Aerodynamic force (world), moment about the CG (body), damping coefficient, and diagnostics."""
    if air is None:
        air = isa(launch.site_elevation + r[2], launch.dT)
    v_air = v - launch.wind(r[2])
    u = Rm.T @ v_air
    V = float(np.linalg.norm(u))
    force, moment = np.zeros(3), np.zeros(3)
    if V < 0.1:
        return force, moment, 0.0, 0.0, 0.0
    mach = V / air.a
    cd = drag_coefficient(rocket, mach, V * rocket.length / air.nu, rocket.motor.thrust(t) > 0)
    qd = 0.5 * air.rho * V * V
    A = rocket.ref_area
    force -= qd * A * cd * v_air / V
    cna, xcp, S1, S2 = rocket.stability(min(mach, 0.9))
    ut = np.array([0.0, u[1], u[2]])
    utm = float(np.linalg.norm(ut))
    alpha = np.arctan2(utm, u[0])
    if utm > 1e-12:
        n_body = -qd * A * cna * alpha * ut / utm
        force += Rm @ n_body
        moment += _cross((xcg - xcp, 0.0, 0.0), n_body)
    damping = 0.5 * air.rho * V * A * (S2 - 2 * xcg * S1 + xcg**2 * cna)
    return force, moment, damping, alpha, mach


def _deriv_free(rocket: Rocket, launch: Launch, t, y):
    r, v, q, w = y[0:3], y[3:6], y[6:10], y[10:13]
    q = q / np.linalg.norm(q)
    Rm = quat_to_mat(q)
    m, xcg, I_p, I_r = rocket.mass_props(t)
    air = isa(launch.site_elevation + r[2], launch.dT)
    force = np.array([0.0, 0.0, -m * launch.gravity]) + rocket.motor.thrust(t, air.P) * Rm[:, 0]
    fa, moment, damp, _, _ = _aero(rocket, launch, t, r, v, Rm, xcg, air)
    force += fa
    damp += rocket.motor.mass_flow(t) * (rocket.nozzle_station - xcg) ** 2
    moment[1:] -= damp * w[1:]
    inertia = np.array([I_r, I_p, I_p])
    w_dot = (moment - _cross(w, inertia * w)) / inertia
    w_dot[0] = 0.0
    q_dot = 0.5 * quat_mul(q, np.array([0.0, *w]))
    return np.concatenate([v, force / m, q_dot, w_dot])


def _deriv_rail(rocket: Rocket, launch: Launch, t, y, u_rail):
    r, v = y[0:3], y[3:6]
    q = y[6:10] / np.linalg.norm(y[6:10])
    Rm = quat_to_mat(q)
    m, xcg, _, _ = rocket.mass_props(t)
    air = isa(launch.site_elevation + r[2], launch.dT)
    force = np.array([0.0, 0.0, -m * launch.gravity]) + rocket.motor.thrust(t, air.P) * Rm[:, 0]
    fa, _, _, _, _ = _aero(rocket, launch, t, r, v, Rm, xcg, air)
    a = float(np.dot(force + fa, u_rail)) / m
    if a < 0 and np.dot(v, u_rail) <= 1e-9:
        a = 0.0  # still sitting on the rail: thrust hasn't beaten weight yet
    return np.concatenate([v, a * u_rail, np.zeros(4), np.zeros(3)])


def _deriv_descent(launch: Launch, m, cda_of, t, y):
    """State: position, velocity, and the air-relative distance traveled since the current canopy deployed."""
    r, v = y[0:3], y[3:6]
    air = isa(launch.site_elevation + r[2], launch.dT)
    v_air = v - launch.wind(r[2])
    speed = float(np.linalg.norm(v_air))
    force = np.array([0.0, 0.0, -m * launch.gravity]) - 0.5 * air.rho * speed * v_air * cda_of(y[6])
    return np.concatenate([v, force / m, [speed]])


def _canopy(base_cda, full_cda, diameter, fill_constant):
    """Drag area as a function of distance traveled since deployment: grows as distance^2 until filled."""
    fill = fill_constant * diameter

    def cda_of(s):
        if fill <= 0 or s >= fill:
            return full_cda
        return base_cda + (full_cda - base_cda) * (max(s, 0.0) / fill) ** 2
    return cda_of, fill


def descend(launch: Launch, recovery: Recovery, m: float, t: float, y: np.ndarray, ev: dict, rec: dict | None = None,
            dt: float = 0.05, dt_fill: float = 0.002, t_end: float | None = None):
    """Point mass under the recovery system, from state y = (position, velocity) at time t until landing.

    Records the peak force while each canopy fills (<canopy>_opening_force_N and _g) in ``ev``.
    """
    use_drogue = bool(recovery.drogue_cda) and recovery.main_altitude is not None
    stage = "drogue" if use_drogue else "main"
    if use_drogue:
        cda_of, fill = _canopy(0.0, recovery.drogue_cda, recovery.drogue_diameter, recovery.fill_constant)
    else:
        cda_of, fill = _canopy(0.0, recovery.main_cda, recovery.main_diameter, recovery.fill_constant)
    yd = np.concatenate([np.asarray(y, float)[0:6], [0.0]])
    peak = 0.0

    def finish_stage():
        ev[f"{stage}_opening_force_N"] = peak
        ev[f"{stage}_opening_g"] = peak / (m * G0)

    while True:
        filling = yd[6] < fill
        h = dt_fill if filling else dt
        y_new = _rk4(lambda tt, yy: _deriv_descent(launch, m, cda_of, tt, yy), t, yd, h)
        if filling:
            air = isa(launch.site_elevation + y_new[2], launch.dT)
            v_air = y_new[3:6] - launch.wind(y_new[2])
            peak = max(peak, 0.5 * air.rho * float(v_air @ v_air) * cda_of(y_new[6]))
            if y_new[6] >= fill:
                finish_stage()
        if stage == "drogue" and y_new[2] <= recovery.main_altitude:
            if yd[6] < fill:
                finish_stage()
            stage, peak = "main", 0.0
            cda_of, fill = _canopy(recovery.drogue_cda, recovery.main_cda, recovery.main_diameter,
                                   recovery.fill_constant)
            ev["main_time"] = t + h
            ev["drogue_descent_speed"] = float(-y_new[5])
            y_new[6] = 0.0
        if y_new[2] <= 0:
            frac = yd[2] / (yd[2] - y_new[2])
            y_land = yd + frac * (y_new - yd)
            t += frac * h
            if y_land[6] < fill:
                finish_stage()
            ev.update(landing_time=t, landing_xy=y_land[0:2].copy(), landing_distance_m=float(np.hypot(*y_land[0:2])),
                      landing_speed=float(np.linalg.norm(y_land[3:6] - launch.wind(0.0))))
            if rec is not None:
                rec_descent(rec, y_land, t)
            return t, y_land[0:6]
        t += h
        yd = y_new
        if rec is not None:
            rec_descent(rec, yd, t)
        if t > 1200 or (t_end is not None and t >= t_end):
            return t, yd[0:6]


def _cross(a, b):
    """3-vector cross product (np.cross has a lot of overhead for tiny arrays)."""
    return np.array([a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]])


def _rk4(f, t, y, h):
    k1 = f(t, y)
    k2 = f(t + h / 2, y + h / 2 * k1)
    k3 = f(t + h / 2, y + h / 2 * k2)
    k4 = f(t + h, y + h * k3)
    return y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


# ---------------------------------------------------------------- simulation

def simulate(rocket: Rocket, launch: Launch = None, recovery: Recovery = None, dt_burn: float = 0.002,
             dt_coast: float = 0.01, dt_descent: float = 0.05, state0: np.ndarray | None = None, t0: float = 0.0,
             t_end: float | None = None, ascent_only: bool = False) -> Flight:
    """Fly the rocket. ``state0``/``t0`` start a free flight from an arbitrary state (used in tests)."""
    launch = launch or Launch()
    recovery = recovery or Recovery()
    u_rail = launch.rail_direction
    if state0 is None:
        q0 = quat_from_to(np.array([1.0, 0.0, 0.0]), u_rail)
        y = np.concatenate([np.zeros(3), np.zeros(3), q0, np.zeros(3)])
        phase = "rail"
    else:
        y = np.asarray(state0, float).copy()
        phase = "free"
    t = t0
    burn = rocket.motor.burn_time
    ev: dict = {}
    rec = {k: [] for k in ("t", "pos", "vel", "tilt", "aoa", "mach", "margin")}

    def record(y, t):
        rec["t"].append(t)
        rec["pos"].append(y[0:3].copy())
        rec["vel"].append(y[3:6].copy())
        Rm = quat_to_mat(y[6:10] / np.linalg.norm(y[6:10]))
        rec["tilt"].append(np.degrees(np.arccos(np.clip(Rm[2, 0], -1, 1))))
        _, xcg, _, _ = rocket.mass_props(t)
        _, _, _, alpha, mach = _aero(rocket, launch, t, y[0:3], y[3:6], Rm, xcg)
        rec["aoa"].append(np.degrees(alpha))
        rec["mach"].append(mach)
        rec["margin"].append(rocket.static_margin(t, min(mach, 0.9)))

    record(y, t)
    max_speed = max_acc = 0.0
    # ---------------- ascent
    while True:
        if t_end is not None and t >= t_end:
            break
        h = dt_burn if t < burn else dt_coast
        if phase == "rail":
            y_new = _rk4(lambda tt, yy: _deriv_rail(rocket, launch, tt, yy, u_rail), t, y, h)
            s_old, s_new = float(np.dot(y[0:3], u_rail)), float(np.dot(y_new[0:3], u_rail))
            if s_new >= launch.rail_length:
                # interpolate to the instant the rocket actually leaves the rail
                frac = (launch.rail_length - s_old) / (s_new - s_old)
                y_exit = y + frac * (y_new - y)
                ev["rail_exit_time"] = t + frac * h
                ev["rail_exit_speed"] = float(np.linalg.norm(y_exit[3:6]))
                w_rel = y_exit[3:6] - launch.wind(y_exit[2])
                ev["rail_exit_aoa_deg"] = float(np.degrees(np.arccos(np.clip(np.dot(w_rel, u_rail) / np.linalg.norm(w_rel), -1, 1))))
                phase = "free"
        else:
            y_new = _rk4(lambda tt, yy: _deriv_free(rocket, launch, tt, yy), t, y, h)
            y_new[6:10] /= np.linalg.norm(y_new[6:10])
        acc = np.linalg.norm(y_new[3:6] - y[3:6]) / h
        max_acc = max(max_acc, acc)
        if phase == "free" and t > 0.1 and y[5] > 0 >= y_new[5] and state0 is None:
            frac = y[5] / (y[5] - y_new[5])
            t_ap = t + frac * h
            y_ap = y + frac * (y_new - y)
            ev.update(apogee_time=t_ap, apogee_m=float(y_ap[2]), apogee_drift_m=float(np.hypot(*y_ap[0:2])),
                      apogee_xy=y_ap[0:2].copy())
            t, y = t_ap, y_ap
            record(y, t)
            break
        t += h
        y = y_new
        max_speed = max(max_speed, float(np.linalg.norm(y[3:6])))
        record(y, t)
        if t > 300:
            break
    ev["max_speed"] = max_speed
    ev["max_accel_g"] = max_acc / G0
    ev["burnout_time"] = burn
    ev["max_mach"] = float(np.nanmax(rec["mach"]))
    if ascent_only or state0 is not None:
        return _pack(rec, ev)

    # ---------------- descent (point mass under parachutes)
    m_dry, _, _, _ = rocket.mass_props(burn + 1.0)
    descend(launch, recovery, m_dry, t, y[0:6], ev, rec, dt=dt_descent)
    return _pack(rec, ev)


def rec_descent(rec, y, t):
    rec["t"].append(t)
    rec["pos"].append(y[0:3].copy())
    rec["vel"].append(y[3:6].copy())
    for k in ("tilt", "aoa", "margin", "mach"):
        rec[k].append(np.nan)


def _pack(rec, ev) -> Flight:
    return Flight(np.array(rec["t"]), np.array(rec["pos"]), np.array(rec["vel"]), np.array(rec["tilt"]),
                  np.array(rec["aoa"]), np.array(rec["mach"]), np.array(rec["margin"]), ev)
