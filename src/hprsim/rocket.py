"""Rocket geometry, mass properties and Barrowman stability.

Stations are measured along the body from the nose tip, positive toward the tail
(the OpenRocket convention). Normal-force slopes are per radian, referenced to
the body cross-section A_ref = pi d^2 / 4.

Barrowman's method (J. Barrowman, "The Practical Calculation of the
Aerodynamic Characteristics of Slender Finned Vehicles", 1967) adds up the
normal-force slope of the nose and the fins; the center of pressure (CP) is
their slope-weighted average position. That linear model ignores the body's own
lift, which grows as sin^2(alpha) and matters at the angles a rocket sees leaving
a short rail in wind. It is added separately (Galejs, "Wind Instability: What
Barrowman Left Out", 1999):  CN_body = K (A_planform / A_ref) sin^2(alpha),  with
K = 1.1, acting at the centroid of the planform area.
"""
from __future__ import annotations

import tomllib
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path

import numpy as np

from .motor import Motor


# ---------------------------------------------------------------- components

@dataclass
class NoseCone:
    """Shapes follow OpenRocket's definitions. ``param`` is the shape parameter where one applies:

    conical, ogive (tangent), ellipsoid
    power      r = R u^n                          (param = n, 0 < n <= 1)
    parabolic  r = R (2u - K u^2) / (2 - K)       (param = K; K = 1 is the full parabola)
    haack      theta = arccos(1 - 2u),  r = R / sqrt(pi) * sqrt(theta - sin(2 theta) / 2 + C sin^3 theta)
               (param = C; 0 is Von Karman, 1/3 is LV-Haack). "vonkarman" is haack with C = 0.
    with u = x / L measured from the tip.
    """

    shape: str
    length: float
    diameter: float
    mass: float
    station: float = 0.0  # tip position
    param: float | None = None

    def radius(self, x):
        """Profile radius at distance x from the tip."""
        L, R = self.length, self.diameter / 2
        x = np.clip(np.asarray(x, float), 0.0, L)
        u = x / L
        if self.shape == "conical":
            return R * u
        if self.shape == "ogive":
            rho = (R**2 + L**2) / (2 * R)
            return np.sqrt(rho**2 - (L - x) ** 2) + R - rho
        if self.shape == "ellipsoid":
            return R * np.sqrt(np.clip(2 * u - u**2, 0.0, None))
        if self.shape == "power":
            return R * u ** (0.5 if self.param is None else self.param)
        if self.shape == "parabolic":
            K = 1.0 if self.param is None else self.param
            return R * (2 * u - K * u**2) / (2 - K)
        if self.shape in ("haack", "vonkarman"):
            C = 0.0 if (self.shape == "vonkarman" or self.param is None) else self.param
            theta = np.arccos(1 - 2 * u)
            return R / np.sqrt(np.pi) * np.sqrt(theta - np.sin(2 * theta) / 2 + C * np.sin(theta) ** 3)
        raise ValueError(f"unknown nose cone shape {self.shape!r}")

    def _profile(self, n=2001):
        """Profile points and trapezoid-rule weights along the arc length (so sum(f * ds) integrates f ds)."""
        x = np.linspace(0, self.length, n)
        r = self.radius(x)
        seg = np.hypot(np.diff(x), np.diff(r))
        ds = np.r_[seg / 2, 0.0] + np.r_[0.0, seg / 2]
        return x, r, ds

    @cached_property
    def volume(self) -> float:
        x, r, _ = self._profile()
        return float(np.trapezoid(np.pi * r**2, x))

    @cached_property
    def wetted_area(self) -> float:
        x, r, ds = self._profile()
        return float(np.sum(2 * np.pi * r * ds))

    @cached_property
    def planform(self) -> tuple[float, float]:
        """Side-view (projected) area and the station of its centroid."""
        x, r, _ = self._profile()
        area = float(np.trapezoid(2 * r, x))
        return area, self.station + float(np.trapezoid(2 * r * x, x)) / area

    def cnalpha_cp(self):
        """Slender-body theory: CN_alpha = 2, CP = L - V / A_base for any profile."""
        A_base = np.pi * (self.diameter / 2) ** 2
        return 2.0, self.station + self.length - self.volume / A_base

    def mass_props(self):
        """Thin uniform shell: mass spread over the surface area."""
        x, r, ds = self._profile()
        dm = 2 * np.pi * r * ds
        dm = self.mass * dm / dm.sum()
        xc = float(np.sum(x * dm) / self.mass)
        I_t = float(np.sum((r**2 / 2 + (x - xc) ** 2) * dm))
        return self.mass, self.station + xc, I_t, float(np.sum(r**2 * dm))


@dataclass
class BodyTube:
    length: float
    diameter: float
    mass: float
    station: float  # fore end

    @property
    def wetted_area(self) -> float:
        return np.pi * self.diameter * self.length

    def mass_props(self):
        R = self.diameter / 2
        return self.mass, self.station + self.length / 2, self.mass * (R**2 / 2 + self.length**2 / 12), self.mass * R**2


@dataclass
class FinSet:
    count: int
    root_chord: float
    tip_chord: float
    span: float  # semi-span, root to tip
    sweep: float  # axial distance from root leading edge to tip leading edge
    thickness: float
    mass: float  # all fins together
    station: float  # root leading edge
    body_diameter: float
    cant: float = 0.0  # rad, each fin's cant angle; positive spins the rocket positively about its nose axis
    cross_section: str = "rounded"  # "square", "rounded" or "airfoil" (sets the edge drag)

    def __post_init__(self):
        if self.cross_section not in ("square", "rounded", "airfoil"):
            raise ValueError(f"fin cross_section must be square, rounded or airfoil, not {self.cross_section!r}")

    @property
    def mid_chord_length(self) -> float:
        Cr, Ct, s, xt = self.root_chord, self.tip_chord, self.span, self.sweep
        return float(np.hypot(s, xt + Ct / 2 - Cr / 2))

    @property
    def planform_area(self) -> float:
        return self.span * (self.root_chord + self.tip_chord) / 2

    def cnalpha_cp(self, ref_diameter: float, mach: float = 0.0):
        Cr, Ct, s, xt, n = self.root_chord, self.tip_chord, self.span, self.sweep, self.count
        R = self.body_diameter / 2
        interference = 1 + R / (s + R)
        # Compressibility enters through the effective aspect ratio (Diederich's semi-empirical form, the one
        # OpenRocket uses), not as a 1/beta factor on the whole slope; capped at M = 0.8 (subsonic theory)
        beta = np.sqrt(1 - min(mach, 0.8) ** 2)
        cna = interference * 4 * n * (s / ref_diameter) ** 2 / (
            1 + np.sqrt(1 + (beta * 2 * self.mid_chord_length / (Cr + Ct)) ** 2))
        xf = xt / 3 * (Cr + 2 * Ct) / (Cr + Ct) + (Cr + Ct - Cr * Ct / (Cr + Ct)) / 6
        return float(cna), self.station + xf

    def roll_coefficients(self, ref_diameter: float, mach: float = 0.0):
        """Roll moment = q A_ref Kf cant - q (p / V) Kd, returning (Kf, Kd).

        Forcing (Barrowman): each canted fin makes normal force CNa1 * cant at its mean aerodynamic chord,
        a distance y_MAC + r_body from the axis. Damping (strip theory, Barrowman): rolling at p, a fin strip
        at radius r meets the air at an extra angle p r / V and pushes back with the 2-D slope 2 pi / beta.
        """
        Cr, Ct, s, n = self.root_chord, self.tip_chord, self.span, self.count
        rt = self.body_diameter / 2
        beta = np.sqrt(1 - min(mach, 0.8) ** 2)
        cna1 = 8 * (s / ref_diameter) ** 2 / (1 + np.sqrt(1 + (beta * 2 * self.mid_chord_length / (Cr + Ct)) ** 2))
        y_mac = s / 3 * (Cr + 2 * Ct) / (Cr + Ct)
        y = np.linspace(0.0, s, 401)
        chord = Cr - (Cr - Ct) * y / s
        kd = n * 2 * np.pi / beta * float(np.trapezoid(chord * (rt + y) ** 2, y))
        return n * cna1 * (y_mac + rt), kd

    def mass_props(self):
        Cr, Ct, xt, s = self.root_chord, self.tip_chord, self.sweep, self.span
        x_bar = (Cr**2 + Cr * Ct + Ct**2 + xt * (Cr + 2 * Ct)) / (3 * (Cr + Ct))  # planform centroid
        r_bar = self.body_diameter / 2 + s * (Cr + 2 * Ct) / (3 * (Cr + Ct))
        chord = (Cr + Ct) / 2
        I_t = self.mass * (r_bar**2 / 2 + chord**2 / 12)
        return self.mass, self.station + x_bar, I_t, self.mass * r_bar**2


@dataclass
class PointMass:
    name: str
    mass: float
    station: float
    inertia: float = 0.0  # transverse, about its own CG (0 for a true point mass)
    roll_inertia: float = 0.0

    def mass_props(self):
        return self.mass, self.station, self.inertia, self.roll_inertia


# ---------------------------------------------------------------- rocket

@dataclass
class Rocket:
    name: str
    nose: NoseCone
    body: BodyTube
    fins: FinSet
    masses: list[PointMass]
    motor: Motor
    motor_aft_station: float  # where the motor's nozzle end sits
    roughness: float = 60e-6  # m, regular paint
    cd_scale: float = 1.0  # multiplier on the drag model (calibrated against flight data)
    body_lift_k: float = 1.1  # Galejs body-lift constant; 0 turns body lift off
    thrust_misalignment_deg: float = 0.0  # angle between the thrust line and the body axis
    thrust_misalignment_azimuth_deg: float = 0.0  # which way it tilts, measured from body +y toward +z
    dry_mass_scale: float = 1.0  # for Monte Carlo
    _dry: tuple = field(init=False, repr=False)
    _aero_cache: dict = field(init=False, repr=False, default_factory=dict)

    def __post_init__(self):
        comps = [self.nose, self.body, self.fins, *self.masses]
        m = s1 = s2 = icm = iroll = 0.0
        for c in comps:
            mi, xi, Ii, Ir = c.mass_props()
            mi, Ii, Ir = mi * self.dry_mass_scale, Ii * self.dry_mass_scale, Ir * self.dry_mass_scale
            m += mi
            s1 += mi * xi
            s2 += mi * xi**2
            icm += Ii
            iroll += Ir
        self._dry = (m, s1, s2, icm, iroll)

    # -- geometry
    @property
    def diameter(self) -> float:
        return self.body.diameter

    @property
    def ref_area(self) -> float:
        return np.pi * self.diameter**2 / 4

    @property
    def length(self) -> float:
        return self.nose.length + self.body.length

    @property
    def nozzle_station(self) -> float:
        return self.motor_aft_station

    @cached_property
    def planform(self) -> tuple[float, float]:
        """Side-view area of nose + body and the station of its centroid (where body lift acts)."""
        a_n, x_n = self.nose.planform
        a_b = self.body.diameter * self.body.length
        x_b = self.body.station + self.body.length / 2
        return a_n + a_b, (a_n * x_n + a_b * x_b) / (a_n + a_b)

    def body_lift_cn(self, alpha: float) -> float:
        """Body normal-force coefficient (on A_ref) at angle of attack alpha, rad."""
        return self.body_lift_k * self.planform[0] / self.ref_area * np.sin(alpha) ** 2

    def cp_at(self, alpha: float, mach: float = 0.0) -> float:
        """Center of pressure including body lift at a finite angle of attack (rad)."""
        cna, xcp, _, _ = self.stability(mach)
        cn_lin, cn_body = cna * alpha, self.body_lift_cn(alpha)
        if cn_lin + cn_body == 0:
            return xcp
        return (cn_lin * xcp + cn_body * self.planform[1]) / (cn_lin + cn_body)

    # -- mass properties
    def dry_mass_props(self):
        """Mass and CG station without the motor."""
        m, s1, *_ = self._dry
        return m, s1 / m

    def mass_props(self, t: float):
        """Total mass, CG station, pitch inertia about the CG, roll inertia."""
        m, s1, s2, icm, iroll = self._dry
        mo = self.motor
        Lm, Rm = mo.length, mo.diameter / 2
        xm = self.motor_aft_station - Lm / 2
        mc, mp = mo.case_mass, mo.prop_remaining(t)
        # case: thin tube; propellant: solid cylinder (BATES grains are close enough)
        m_tot = m + mc + mp
        s1 += (mc + mp) * xm
        s2 += (mc + mp) * xm**2
        icm += mc * (Rm**2 / 2 + Lm**2 / 12) + mp * (3 * Rm**2 + Lm**2) / 12
        iroll += mc * Rm**2 + mp * Rm**2 / 2
        x_cg = s1 / m_tot
        I_pitch = icm + s2 - m_tot * x_cg**2
        return m_tot, x_cg, I_pitch, iroll

    # -- aerodynamics
    def stability(self, mach: float = 0.0):
        """Total CN_alpha, CP station, and the sums needed for pitch damping."""
        key = round(mach, 2)
        if key not in self._aero_cache:
            parts = [self.nose.cnalpha_cp(), self.fins.cnalpha_cp(self.diameter, key)]
            cna = sum(p[0] for p in parts)
            xcp = sum(p[0] * p[1] for p in parts) / cna
            S1 = sum(p[0] * p[1] for p in parts)
            S2 = sum(p[0] * p[1] ** 2 for p in parts)
            self._aero_cache[key] = (cna, xcp, S1, S2)
        return self._aero_cache[key]

    def roll_aero(self, mach: float = 0.0):
        """(Kf, Kd) for the fin set at this Mach number (see FinSet.roll_coefficients)."""
        key = ("roll", round(mach, 2))
        if key not in self._aero_cache:
            self._aero_cache[key] = self.fins.roll_coefficients(self.diameter, key[1])
        return self._aero_cache[key]

    @cached_property
    def thrust_direction(self) -> np.ndarray:
        """Unit thrust vector in the body frame (x along the axis toward the nose)."""
        e, a = np.radians(self.thrust_misalignment_deg), np.radians(self.thrust_misalignment_azimuth_deg)
        return np.array([np.cos(e), np.sin(e) * np.cos(a), np.sin(e) * np.sin(a)])

    def static_margin(self, t: float = 0.0, mach: float = 0.0) -> float:
        """(CP - CG) / diameter, in calibers. Positive = stable."""
        _, xcg, _, _ = self.mass_props(t)
        return (self.stability(mach)[1] - xcg) / self.diameter


# ---------------------------------------------------------------- loading

def load(path: str | Path, motor: Motor | None = None, motor_file: str | Path | None = None, **overrides) -> Rocket:
    """Build a Rocket from a TOML description (see rockets/example_3in.toml) or an OpenRocket .ork file."""
    path = Path(path)
    if path.suffix.lower() == ".ork":
        from . import ork
        cfg = ork.read(path)
    else:
        cfg = tomllib.loads(path.read_text())
    return from_config(cfg, path.parent, motor, motor_file, **overrides)


def find_motor(designation: str, *dirs: Path) -> Path:
    """The .eng file whose name contains the motor designation (e.g. 'J420R' -> AeroTech_J420R.eng)."""
    for d in dirs:
        for f in sorted(Path(d).glob("*.eng")) if Path(d).is_dir() else []:
            if designation and designation.lower() in f.stem.lower():
                return f
    raise FileNotFoundError(f"no .eng file for motor {designation!r}: download it from thrustcurve.org into "
                            f"data/motors/ or pass motor_file=")


def from_config(cfg: dict, base_dir: Path, motor: Motor | None = None, motor_file: str | Path | None = None,
                **overrides) -> Rocket:
    d = cfg["body"]["diameter"]
    nose = NoseCone(cfg["nose"]["shape"], cfg["nose"]["length"], d, cfg["nose"]["mass"],
                    param=cfg["nose"].get("param"))
    body = BodyTube(cfg["body"]["length"], d, cfg["body"]["mass"], station=nose.length)
    f = cfg["fins"]
    fin_station = nose.length + body.length - f["root_chord"] - f.get("aft_offset", 0.0)
    fins = FinSet(f["count"], f["root_chord"], f["tip_chord"], f["span"], f["sweep"], f["thickness"], f["mass"],
                  fin_station, d, cant=np.radians(f.get("cant_deg", 0.0)),
                  cross_section=f.get("cross_section", "rounded"))
    masses = [PointMass(k, v["mass"], v["station"], v.get("inertia", 0.0), v.get("roll_inertia", 0.0))
              for k, v in cfg.get("mass", {}).items()]
    mcfg = cfg["motor"]
    if motor is None:
        if motor_file:
            mf = motor_file
        elif mcfg.get("file"):
            mf = base_dir / mcfg["file"]
        else:
            mf = find_motor(mcfg.get("designation", ""), base_dir, base_dir / "data" / "motors",
                            base_dir.parent / "data" / "motors", Path.cwd() / "data" / "motors")
        motor = Motor.from_eng(mf)
    d_exit = mcfg.get("nozzle_exit_diameter", 0.0)
    if d_exit and not motor.exit_area:
        motor.exit_area = np.pi * d_exit**2 / 4
    aft = mcfg.get("aft_station", nose.length + body.length + mcfg.get("overhang", 0.0))
    kwargs = {"roughness": cfg.get("roughness", 60e-6),
              "thrust_misalignment_deg": mcfg.get("misalignment_deg", 0.0),
              "thrust_misalignment_azimuth_deg": mcfg.get("misalignment_azimuth_deg", 0.0)}
    kwargs.update(overrides)
    return Rocket(cfg.get("name", "rocket"), nose, body, fins, masses, motor, aft, **kwargs)
