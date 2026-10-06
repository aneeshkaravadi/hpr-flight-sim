"""Every test compares against an answer worked out independently of the simulator."""
import json
from pathlib import Path

import numpy as np
import pytest

from hprsim import aero, atmosphere, flight, motor, rocket

ROOT = Path(__file__).resolve().parents[1]
ROCKET = ROOT / "rockets" / "example_3in.toml"
META = json.loads((ROOT / "data" / "motors" / "thrustcurve_metadata.json").read_text())
G = atmosphere.G0


def constant_thrust_motor(F=300.0, burn=2.0, prop=1e-6, case=0.3):
    return motor.Motor("const", 0.038, 0.3, prop, prop + case, [0.0, 1e-4, burn, burn + 1e-4], [0.0, F, F, 0.0])


# ---------------------------------------------------------------- atmosphere and motors

def test_isa_table_values():
    assert atmosphere.isa(0).rho == pytest.approx(1.2250, rel=1e-4)
    a5 = atmosphere.isa(5000)
    assert a5.T == pytest.approx(255.65, abs=1e-6)
    assert a5.P == pytest.approx(54020, rel=1e-3)  # US Standard Atmosphere 1976
    a11 = atmosphere.isa(11000)
    assert a11.T == pytest.approx(216.65, abs=1e-6)
    assert a11.P == pytest.approx(22632, rel=1e-3)
    assert atmosphere.isa(0).a == pytest.approx(340.29, abs=0.05)


@pytest.mark.parametrize("name", ["AeroTech_H128W", "AeroTech_I284W", "AeroTech_J420R"])
def test_eng_files_match_published_impulse(name):
    """The .eng files are ~25-point samples of the certification data, so allow 3%."""
    m = motor.Motor.from_eng(ROOT / "data" / "motors" / f"{name}.eng")
    assert m.total_impulse == pytest.approx(META[name]["totImpulseNs"], rel=0.03)


def test_mass_flow_burns_exactly_the_propellant():
    m = motor.Motor.from_eng(ROOT / "data" / "motors" / "AeroTech_I284W.eng")
    t = np.linspace(0, m.burn_time, 20001)
    burned = np.trapezoid([m.mass_flow(x) for x in t], t)
    assert burned == pytest.approx(m.prop_mass, rel=1e-3)
    assert m.prop_remaining(m.burn_time) == 0.0


def test_pressure_thrust_correction_only_while_burning():
    m = motor.Motor.from_eng(ROOT / "data" / "motors" / "AeroTech_I284W.eng")
    m.exit_area = np.pi * 0.015**2 / 4
    p = atmosphere.isa(1500.0).P
    assert m.thrust(1.0, p) - m.thrust(1.0) == pytest.approx((m.p_ref - p) * m.exit_area, rel=1e-12)
    assert m.thrust(m.burn_time + 0.1, p) == 0.0
    t = np.linspace(0, m.burn_time, 20001)  # pressure thrust doesn't burn extra propellant
    assert np.trapezoid([m.mass_flow(x) for x in t], t) == pytest.approx(m.prop_mass, rel=1e-3)


def test_pressure_thrust_adds_the_expected_speed():
    """No gravity or drag, constant mass: extra speed = (A_e / m) * integral of (p_ref - p) dt while burning."""
    Ae = np.pi * 0.02**2 / 4
    plain, corrected = constant_thrust_motor(F=400.0, burn=1.0), constant_thrust_motor(F=400.0, burn=1.0)
    corrected.exit_area = Ae
    launch = flight.Launch(gravity=0.0, site_elevation=3000.0)
    runs = [flight.simulate(rocket.load(ROCKET, motor=m, cd_scale=0.0), launch, t_end=1.0, ascent_only=True)
            for m in (plain, corrected)]
    dv = np.linalg.norm(runs[1].vel[-1]) - np.linalg.norm(runs[0].vel[-1])
    mass = rocket.load(ROCKET, motor=corrected).mass_props(0)[0]
    t, z = runs[1].t, runs[1].pos[:, 2]
    dp = [(corrected.p_ref - atmosphere.isa(3000.0 + zi).P) * (corrected.thrust(ti) > 0) for ti, zi in zip(t, z)]
    assert dv == pytest.approx(Ae / mass * np.trapezoid(dp, t), rel=2e-3)


def test_toml_nozzle_exit_diameter_sets_exit_area(tmp_path):
    text = ROCKET.read_text().replace('overhang = 0.0', 'overhang = 0.0\nnozzle_exit_diameter = 0.016')
    text = text.replace('file = "../data/', f'file = "{ROOT}/data/')
    f = tmp_path / "r.toml"
    f.write_text(text)
    assert rocket.load(f).motor.exit_area == pytest.approx(np.pi * 0.016**2 / 4)


# ---------------------------------------------------------------- Barrowman

def test_conical_nose_cp_is_two_thirds_length():
    n = rocket.NoseCone("conical", 0.3, 0.08, 0.1)
    cna, xcp = n.cnalpha_cp()
    assert cna == 2.0
    assert xcp == pytest.approx(0.2, rel=1e-4)


@pytest.mark.parametrize("shape, expected", [("ogive", 0.466), ("vonkarman", 0.5)])
def test_nose_cp_matches_barrowman_table(shape, expected):
    n = rocket.NoseCone(shape, 0.38, 0.0787, 0.1)
    assert n.cnalpha_cp()[1] / 0.38 == pytest.approx(expected, rel=0.04)


def test_parabolic_series_nose_cp_is_seven_fifteenths():
    """r = R (2u - u^2): volume = 8/15 pi R^2 L, so CP = L - V/A = 7/15 L.

    (Not Barrowman's tabulated 'parabola', which is the r ~ sqrt(x) shape with CP = L/2.)
    """
    n = rocket.NoseCone("parabolic", 0.38, 0.0787, 0.1)
    assert n.cnalpha_cp()[1] == pytest.approx(0.38 * 7 / 15, rel=1e-4)


@pytest.mark.parametrize("shape, param, cp_fraction", [
    ("ellipsoid", None, 1 / 3),  # V = 2/3 pi R^2 L
    ("power", 0.5, 0.5),  # V = pi R^2 L / (2n + 1), so CP = 2n / (2n + 1) L
    ("power", 0.75, 0.6),
    ("parabolic", 0.5, None),
    ("haack", 0.0, 0.5),  # Von Karman: V = pi R^2 L / 2
])
def test_nose_shapes_cp_from_volume(shape, param, cp_fraction):
    n = rocket.NoseCone(shape, 0.4, 0.1, 0.1, param=param)
    if cp_fraction is None:  # parabolic series, K = 1/2: integrate r^2 = R^2 (2u - u^2/2)^2 / (3/2)^2 by hand
        cp_fraction = 1 - (4 / 3 - 1 / 2 + 1 / 20) / (9 / 4)
    assert n.cnalpha_cp()[1] == pytest.approx(cp_fraction * 0.4, rel=1e-4)


def test_rectangular_fin_cp_is_quarter_chord_and_slope_matches_hand_calc():
    c, s, d, n = 0.10, 0.08, 0.08, 4
    f = rocket.FinSet(n, c, c, s, 0.0, 0.003, 0.1, station=1.0, body_diameter=d)
    cna, xcp = f.cnalpha_cp(d)
    assert xcp == pytest.approx(1.0 + c / 4, abs=1e-12)
    interference = 1 + (d / 2) / (s + d / 2)
    by_hand = interference * 4 * n * (s / d) ** 2 / (1 + np.sqrt(1 + (s / c) ** 2))
    assert cna == pytest.approx(by_hand, rel=1e-12)


def test_cone_planform_area_and_centroid():
    """Side view of a cone is a triangle: area R * L, centroid 2/3 of the way back."""
    n = rocket.NoseCone("conical", 0.3, 0.08, 0.1)
    area, x = n.planform
    assert area == pytest.approx(0.04 * 0.3, rel=1e-6)
    assert x == pytest.approx(0.2, rel=1e-6)


def test_body_lift_force_and_moment_follow_galejs():
    """Body lift adds q * K * A_plan * sin^2(alpha) against the crossflow, acting at the planform centroid."""
    with_lift, without = rocket.load(ROCKET), rocket.load(ROCKET, body_lift_k=0.0)
    launch = flight.Launch(site_elevation=0.0)
    V, alpha = 50.0, np.radians(20.0)
    v = V * np.array([np.cos(alpha), np.sin(alpha), 0.0])  # body axis along world x, crossflow along +y
    xcg = with_lift.mass_props(0.0)[1]
    f1, m1, *_ = flight._aero(with_lift, launch, 1.0, np.zeros(3), v, np.eye(3), xcg)
    f0, m0, *_ = flight._aero(without, launch, 1.0, np.zeros(3), v, np.eye(3), xcg)
    q = 0.5 * atmosphere.isa(0.0).rho * V**2
    a_plan, x_plan = with_lift.planform
    N = q * 1.1 * a_plan * np.sin(alpha) ** 2
    assert f1 - f0 == pytest.approx([0.0, -N, 0.0], abs=1e-9 * N)
    assert (m1 - m0)[2] == pytest.approx(-(xcg - x_plan) * N, rel=1e-9)


def test_body_lift_moves_the_cp_forward_at_high_angle():
    r = rocket.load(ROCKET)
    xcp = r.stability(0.0)[1]
    assert r.cp_at(np.radians(0.01)) == pytest.approx(xcp, abs=1e-3)  # its share grows with alpha, so tiny here
    assert r.cp_at(np.radians(20.0)) < xcp - 0.05  # the planform centroid is well ahead of the fins


def test_fin_slope_compressibility_acts_through_the_aspect_ratio():
    """beta = sqrt(1 - M^2) multiplies the 2 L_F / (Cr + Ct) term, which grows the slope far less than 1/beta would."""
    c, s, d, n, M = 0.10, 0.08, 0.08, 4, 0.6
    f = rocket.FinSet(n, c, c, s, 0.0, 0.003, 0.1, station=1.0, body_diameter=d)
    beta = np.sqrt(1 - M**2)
    interference = 1 + (d / 2) / (s + d / 2)
    by_hand = interference * 4 * n * (s / d) ** 2 / (1 + np.sqrt(1 + (beta * s / c) ** 2))
    assert f.cnalpha_cp(d, M)[0] == pytest.approx(by_hand, rel=1e-12)
    assert f.cnalpha_cp(d, M)[0] / f.cnalpha_cp(d, 0.0)[0] < 1 / beta


def test_example_rocket_is_stable_and_margin_grows_during_burn():
    r = rocket.load(ROCKET)
    assert 1.0 < r.static_margin(0.0) < r.static_margin(r.motor.burn_time + 1)


def test_example_drag_is_in_the_usual_range():
    """Typical subsonic HPR drag coefficients from OpenRocket are about 0.4-0.6."""
    r = rocket.load(ROCKET)
    cd = aero.drag_coefficient(r, 0.3, 1.0e7, thrusting=False)
    assert 0.4 < cd < 0.6


# ---------------------------------------------------------------- flight physics

def test_vertical_flight_without_drag_matches_kinematics():
    F, burn = 300.0, 2.0
    r = rocket.load(ROCKET, motor=constant_thrust_motor(F, burn), cd_scale=0.0)
    m = r.mass_props(0)[0]
    a = F / m - G
    v_b, h_b = a * burn, 0.5 * a * burn**2
    fl = flight.simulate(r, flight.Launch(), ascent_only=True)
    assert fl.apogee == pytest.approx(h_b + v_b**2 / (2 * G), rel=2e-3)
    assert fl.events["rail_exit_speed"] == pytest.approx(np.sqrt(2 * a * 1.83), rel=2e-3)


def test_burnout_speed_matches_rocket_equation():
    """No gravity, no drag: delta-v = c ln(m0 / mf) (Tsiolkovsky)."""
    mot = constant_thrust_motor(F=500.0, burn=2.0, prop=0.4)
    r = rocket.load(ROCKET, motor=mot, cd_scale=0.0)
    m0, mf = r.mass_props(0)[0], r.mass_props(3)[0]
    fl = flight.simulate(r, flight.Launch(gravity=0.0), t_end=2.2, ascent_only=True)
    v = np.linalg.norm(fl.vel[-1])
    assert v == pytest.approx(mot.exhaust_velocity * np.log(m0 / mf), rel=2e-3)


def test_pitch_oscillation_frequency_and_damping():
    """Coasting at constant speed: the body oscillates about the flight path.

    Linearized short-period motion (DERIVATIONS.md section 7), with K = q A CN_alpha (x_cp - x_cg)
    and Z = q A CN_alpha / (m V):  alpha'' + (C/I + Z) alpha' + (K/I + C Z / I) alpha = 0,
    so it decays at sigma = C / (2 I) + Z / 2 and rings at omega_d = sqrt(K/I + C Z/I - sigma^2).
    """
    r = rocket.load(ROCKET, cd_scale=0.0)
    launch = flight.Launch(gravity=0.0, site_elevation=0.0)
    V, z = 150.0, 1000.0
    t0 = r.motor.burn_time + 1.0  # motor burnt out, mass constant
    y0 = np.concatenate([[0, 0, z], [V, 0, 0], [1.0, 0, 0, 0], [0, 0.2, 0]])
    fl = flight.simulate(r, launch, state0=y0, t0=t0, t_end=t0 + 3.0, dt_coast=0.0005)
    theta = 90.0 - fl.tilt_deg  # pitch of the body above the horizontal, degrees
    gamma = np.degrees(np.arctan2(fl.vel[:, 2], fl.vel[:, 0]))
    aoa = theta - gamma
    crossings = np.where(np.diff(np.sign(aoa)) != 0)[0]
    period = 2 * np.mean(np.diff(fl.t[crossings]))

    air = atmosphere.isa(z)
    m, xcg, I, _ = r.mass_props(t0)
    cna, xcp, S1, S2 = r.stability(V / air.a)
    q, A = 0.5 * air.rho * V**2, r.ref_area
    K = q * A * cna * (xcp - xcg)
    c_aero = 0.5 * air.rho * V * A * (S2 - 2 * xcg * S1 + xcg**2 * cna)
    Z = q * A * cna / (m * V)
    sigma = c_aero / (2 * I) + Z / 2
    omega_d = np.sqrt(K / I + c_aero * Z / I - sigma**2)
    assert period == pytest.approx(2 * np.pi / omega_d, rel=0.02)

    peaks = [np.abs(aoa[a:b]).max() for a, b in zip(crossings[:-1], crossings[1:])]
    decay = np.log(peaks[0] / peaks[4]) / (4 * period / 2)
    assert decay == pytest.approx(sigma, rel=0.1)


def test_canted_fins_spin_up_to_the_steady_roll_rate():
    """Coasting at constant speed: forcing q A Kf cant balances damping q (p / V) Kd at p = A Kf cant V / Kd,
    approached with time constant I_roll V / (q Kd)."""
    r = rocket.load(ROCKET, cd_scale=0.0)
    r.fins.cant = np.radians(1.0)
    launch = flight.Launch(gravity=0.0, site_elevation=0.0)
    V, z = 100.0, 1000.0
    t0 = r.motor.burn_time + 1.0
    y0 = np.concatenate([[0, 0, z], [V, 0, 0], [1.0, 0, 0, 0], [0, 0, 0]])
    fl = flight.simulate(r, launch, state0=y0, t0=t0, t_end=t0 + 1.0, dt_coast=0.001)
    air = atmosphere.isa(z)
    kf, kd = r.roll_aero(V / air.a)
    p_ss = r.ref_area * kf * r.fins.cant * V / kd
    tau = r.mass_props(t0)[3] * V / (0.5 * air.rho * V**2 * kd)
    assert tau < 0.1  # so one second is many time constants
    assert fl.roll_rate[-1] == pytest.approx(p_ss, rel=1e-3)
    k = np.searchsorted(fl.t, t0 + tau)  # first-order response: 1 - 1/e of the way there after one tau
    assert fl.roll_rate[k] == pytest.approx(p_ss * (1 - np.exp(-(fl.t[k] - t0) / tau)), rel=0.02)


def test_thrust_misalignment_torques_the_rocket_about_the_cg():
    """At rest (no aero), a thrust line tilted by eps at the nozzle gives force F (cos eps, sin eps, 0) and
    yaw moment (x_cg - x_nozzle) F sin eps, which pushes the nose the other way."""
    eps = np.radians(0.5)
    r = rocket.load(ROCKET, cd_scale=0.0, thrust_misalignment_deg=0.5)
    launch = flight.Launch(gravity=0.0)
    t = 0.5
    y = np.concatenate([[0, 0, 100.0], [0, 0, 0], [1.0, 0, 0, 0], [0, 0, 0]])  # body axis along world x
    d = flight._deriv_free(r, launch, t, y)
    m, xcg, I_p, _ = r.mass_props(t)
    F = r.motor.thrust(t, atmosphere.isa(launch.site_elevation + 100.0).P)
    assert d[3:6] == pytest.approx(F * np.array([np.cos(eps), np.sin(eps), 0.0]) / m, rel=1e-9)
    assert d[12] == pytest.approx((xcg - r.nozzle_station) * F * np.sin(eps) / I_p, rel=1e-9)
    assert d[12] < 0


def test_rocket_weathercocks_into_the_wind():
    r = rocket.load(ROCKET)
    fl = flight.simulate(r, flight.Launch(wind_speed=6.0, wind_from_deg=270.0))  # wind from the west
    assert fl.events["apogee_xy"][0] < -10.0  # apogee is upwind (west)
    assert fl.events["landing_xy"][0] > fl.events["apogee_xy"][0]  # then it drifts back downwind


def test_descent_reaches_terminal_velocity_under_main():
    r = rocket.load(ROCKET)
    rec = flight.Recovery()
    fl = flight.simulate(r, flight.Launch(site_elevation=0.0), rec)
    m = r.mass_props(100)[0]
    v_t = np.sqrt(2 * m * G / (atmosphere.isa(0).rho * rec.main_cda))
    assert fl.events["landing_speed"] == pytest.approx(v_t, rel=0.02)


def test_opening_shock_matches_the_closed_form():
    """Gravity-free canopy whose drag area grows as (x / L)^2 over the fill distance L:
    m v dv/dx = -1/2 rho v^2 C (x/L)^2  gives  v = v0 exp(-k x^3),  k = rho C / (6 m L^2),
    so the force peaks at x* = (3k)^(-1/3), or at full inflation when x* > L (a heavy payload).
    """
    launch = flight.Launch(gravity=0.0, site_elevation=0.0)
    rho = atmosphere.isa(1000.0).rho
    C, D, n, v0 = 1.41, 1.524, 8.0, 30.0
    L = n * D
    for m in (2.0, 200.0):
        rec = flight.Recovery(drogue_cda=None, main_altitude=None, main_cda=C, main_diameter=D, fill_constant=n)
        ev = {}
        flight.descend(launch, rec, m, 0.0, np.array([0.0, 0.0, 1000.0, v0, 0.0, 0.0]), ev, t_end=5.0, dt_fill=5e-4)
        k = rho * C / (6 * m * L**2)
        x = min((1 / (3 * k)) ** (1 / 3), L)
        expected = 0.5 * rho * C * v0**2 * (x / L) ** 2 * np.exp(-2 * k * x**3)
        assert ev["main_opening_force_N"] == pytest.approx(expected, rel=2e-3)


def test_time_step_is_converged():
    r = rocket.load(ROCKET)
    a = flight.simulate(r, flight.Launch(wind_speed=4.0), ascent_only=True).apogee
    b = flight.simulate(r, flight.Launch(wind_speed=4.0), dt_burn=0.001, dt_coast=0.0025, ascent_only=True).apogee
    assert a == pytest.approx(b, rel=1e-3)


# ---------------------------------------------------------------- calibration

def test_drag_calibration_recovers_a_known_multiplier(tmp_path):
    """Fly a 'real' rocket with 15% more drag, log altitude in feet, and fit it back."""
    from hprsim import calibrate
    launch = flight.Launch(wind_speed=3.0)
    truth = flight.simulate(rocket.load(ROCKET, cd_scale=1.15), launch, ascent_only=True)
    log = tmp_path / "flight.csv"
    rows = ["Time (s),Altitude (ft)"] + [f"{t:.3f},{z / calibrate.FT:.2f}" for t, z in zip(truth.t, truth.pos[:, 2])]
    log.write_text("\n".join(rows))
    t, alt = calibrate.read_altitude_log(log)
    assert alt.max() == pytest.approx(truth.apogee, rel=1e-4)
    scale = calibrate.fit_cd_scale(lambda s: rocket.load(ROCKET, cd_scale=s), launch, float(alt.max()))
    assert scale == pytest.approx(1.15, rel=5e-3)


def test_reads_openrocket_csv_export(tmp_path):
    from hprsim import calibrate
    f = tmp_path / "or.csv"
    f.write_text("# OpenRocket export\n# Time (s),Altitude (m),Vertical velocity (m/s)\n0,0,0\n1.5,120.0,90\n# Event APOGEE\n9,800.5,0\n")
    t, alt = calibrate.read_altitude_log(f)
    assert list(t) == [0, 1.5, 9] and alt.max() == pytest.approx(800.5)
