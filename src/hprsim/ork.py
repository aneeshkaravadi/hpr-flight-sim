"""Read OpenRocket designs (.ork).

An .ork file is XML, zipped by newer OpenRocket versions and gzipped by older
ones. This reads the parts hprsim models:

  nose cone        shape, length, wall thickness and material
  body tubes       a stack of same-diameter tubes, merged into one for the aerodynamics
  fins             one trapezoidal fin set
  everything else  couplers, bulkheads, rings, inner tubes, mass components, parachutes,
                   shock cords, launch lugs and rail buttons, each a point mass at its own CG
  parachutes       Cd, diameter, and when each one deploys
  motor            the default configuration's motor, and where its nozzle ends up

Masses are computed from each part's material and geometry the way OpenRocket
does it, then the mass overrides saved in the file are applied (an override
that includes subcomponents scales the whole subtree, keeping its CG).
Transitions between diameters and secant ogives aren't supported yet.

read() returns the same structure as a rocket TOML file, so rocket.load() can
fly an .ork directly and examples/import_ork.py can write it out as TOML.
"""
from __future__ import annotations

import gzip
import io
import re
import xml.etree.ElementTree as ET
import zipfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .rocket import NoseCone

AUTO_CHUTE_CD = 0.8  # what OpenRocket uses when a parachute's Cd is left on "auto"
FINISH_ROUGHNESS = {"rough": 500e-6, "unfinished": 150e-6, "normal": 60e-6, "smooth": 20e-6, "polished": 2e-6}
TUBE_LIKE = {"tubecoupler", "innertube", "launchlug"}
PACKED = {"parachute", "streamer", "shockcord", "masscomponent"}


def _parse(path) -> ET.Element:
    raw = Path(path).read_bytes()
    if raw[:2] == b"PK":
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            raw = z.read(next(n for n in z.namelist() if n.endswith((".ork", ".xml"))))
    elif raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    return ET.fromstring(raw)


def _num(el: ET.Element, tag: str, default=None):
    """A number, or ("auto", last computed value) for auto-sized dimensions like 'auto 0.025'."""
    s = el.findtext(tag)
    if s is None:
        return default
    s = s.strip()
    if s.lower().startswith("auto"):
        rest = s[4:].strip()
        return ("auto", float(rest) if rest else None)
    return float(s)


def _fixed(v, fallback):
    if isinstance(v, tuple):
        return fallback if fallback is not None else v[1]
    return fallback if v is None else v


def _density(el: ET.Element, tag: str = "material") -> float:
    m = el.find(tag)
    return float(m.get("density")) if m is not None and m.get("density") else 0.0


def _placement(el: ET.Element):
    a = el.find("axialoffset")
    if a is not None:
        return (a.get("method") or "after").lower(), float(a.text or 0.0)
    p = el.find("position")
    if p is not None:
        return (p.get("type") or "after").lower(), float(p.text or 0.0)
    return "after", 0.0


@dataclass
class Part:
    """One mass item: a component's own mass (subcomponents are separate parts)."""

    name: str
    kind: str
    mass: float
    cg: float  # station of its CG, m from the nose tip
    inertia: float = 0.0  # transverse, about its own CG
    roll_inertia: float = 0.0
    tags: dict = field(default_factory=dict)


class _Reader:
    def __init__(self, root: ET.Element):
        self.root = root
        self.warnings: list[str] = []
        self.parts: list[Part] = []
        self.chutes: list[dict] = []
        self.mount = None  # (aft station of the mount, overhang, motor element)
        rk = root.find("rocket")
        self.rocket = rk
        default = [m.get("configid") for m in rk.findall("motorconfiguration") if m.get("default") == "true"]
        self.config_id = default[0] if default else None

    # -- the stack along the axis: nose cone, body tubes
    def read(self):
        stages = self.rocket.findall("subcomponents/stage")
        if len(stages) != 1:
            raise ValueError(f"{len(stages)} stages found; only single-stage rockets are supported")
        stage = stages[0]
        items = list(stage.find("subcomponents"))
        tubes = [c for c in items if c.tag == "bodytube"]
        if any(c.tag == "transition" for c in items):
            raise ValueError("transitions (boattails or changes in diameter) aren't supported yet")
        fixed_r = [r for r in (_num(t, "radius") for t in tubes) if not isinstance(r, tuple) and r is not None]
        auto_r = [r[1] for r in (_num(t, "radius") for t in tubes) if isinstance(r, tuple) and r[1]]
        R = max(fixed_r) if fixed_r else max(auto_r)
        if fixed_r and (max(fixed_r) - min(fixed_r)) > 1e-6:
            raise ValueError("body tubes of different diameters aren't supported yet")
        self.radius = R

        x = 0.0
        self.nose = None
        self.tube_spans = []
        stage_start = len(self.parts)
        for c in items:
            if c.tag == "nosecone":
                self.nose = self._nose(c, x)
                x += self.nose["length"]
            elif c.tag == "bodytube":
                L = float(c.findtext("length"))
                start = len(self.parts)
                self._tube(c, x, L)
                self._children(c, x, L, self.radius - float(c.findtext("thickness") or 0.0))
                if c.find("motormount") is not None and self.mount is None:
                    mm = c.find("motormount")
                    self.mount = (x + L, float(mm.findtext("overhang") or 0.0), mm)
                self._override(c, start)
                self.tube_spans.append((x, L, c.findtext("finish") or "normal"))
                x += L
            else:
                self.warnings.append(f"skipped {c.tag} '{c.findtext('name')}' at the top level")
        if self.nose is None:
            raise ValueError("no nose cone found")
        self.length = x
        if (stage.findtext("overridesubcomponentsmass") or "false").lower() == "true":
            self._override(stage, stage_start)
        elif stage.findtext("overridemass") is not None:
            self.warnings.append("stage mass override without subcomponents ignored")
        return self

    def _nose(self, c, x):
        L = float(c.findtext("length"))
        shape, param = self._nose_shape(c)
        start = len(self.parts)
        cone = NoseCone(shape, L, 2 * self.radius, 1.0, station=x, param=param)
        t = float(c.findtext("thickness") or 0.0)
        filled = (c.findtext("filled") or "false").lower() == "true"
        rho = _density(c)
        mass = rho * (cone.volume if filled or t <= 0 else cone.wetted_area * t)
        _, cg, inertia, roll = NoseCone(shape, L, 2 * self.radius, max(mass, 1e-12), station=x, param=param).mass_props()
        self.parts.append(Part(c.findtext("name") or "nose cone", "nose", mass, cg, inertia, roll,
                               {"shape": shape, "param": param, "length": L}))
        self._children(c, x, L, self.radius - t)
        self._override(c, start)
        return {"shape": shape, "param": param, "length": L}

    def _nose_shape(self, c):
        shape = (c.findtext("shape") or "ogive").strip().lower()
        p = c.findtext("shapeparameter")
        p = float(p) if p is not None else None
        if shape == "ogive":
            if p is None or abs(p - 1.0) < 1e-9:
                return "ogive", None
            if p == 0.0:
                return "conical", None
            raise ValueError("secant ogives (ogive shape parameter below 1) aren't supported yet")
        if shape in ("conical", "ellipsoid"):
            return shape, None
        if shape in ("power", "parabolic", "haack"):
            return shape, p
        raise ValueError(f"nose cone shape {shape!r} isn't supported")

    def _tube(self, c, x, L):
        R = self.radius
        t = float(c.findtext("thickness") or 0.0)
        m = _density(c) * np.pi * (R**2 - (R - t) ** 2) * L
        self.parts.append(Part(c.findtext("name") or "body tube", "tube", m, x + L / 2,
                               m * (R**2 / 2 + L**2 / 12), m * R**2))

    # -- everything attached to (or inside) a nose cone or body tube
    def _children(self, parent, px, plen, inner_r):
        sub = parent.find("subcomponents")
        if sub is None:
            return
        inner_tubes = [c for c in sub if c.tag == "innertube"]
        mount_r = max((float(_fixed(_num(c, "outerradius"), 0.0)) for c in inner_tubes), default=0.0)
        prev_end = px
        for c in sub:
            start = len(self.parts)
            own_len = self._length(c)
            method, off = _placement(c)
            if method == "top":
                cx = px + off
            elif method == "bottom":
                cx = px + plen - own_len + off
            elif method == "middle":
                cx = px + (plen - own_len) / 2 + off
            elif method == "absolute":
                cx = off
            else:
                cx = prev_end + off
            prev_end = cx + own_len
            self._component(c, cx, own_len, inner_r, mount_r)
            self._override(c, start)

    def _length(self, c):
        if c.tag == "trapezoidfinset":
            return float(c.findtext("rootchord"))
        if c.tag in PACKED:
            return float(_fixed(_num(c, "packedlength"), 0.0))
        if c.tag == "railbutton":
            return 0.0
        return float(_fixed(_num(c, "length"), 0.0))

    def _component(self, c, x, L, inner_r, mount_r):
        name, kind, rho = c.findtext("name") or c.tag, c.tag, _density(c)
        n = int(float(c.findtext("instancecount") or 1))
        if kind == "trapezoidfinset":
            self._fins(c, x)
            return
        if kind in TUBE_LIKE:
            ro = float(_fixed(_num(c, "outerradius"), inner_r))
            t = float(c.findtext("thickness") or 0.0)
            m = n * rho * np.pi * (ro**2 - max(ro - t, 0.0) ** 2) * L
            self.parts.append(Part(name, kind, m, x + L / 2, m * (ro**2 / 2 + L**2 / 12), m * ro**2))
            if kind == "innertube":
                mm = c.find("motormount")
                if mm is not None:
                    self.mount = (x + L, float(mm.findtext("overhang") or 0.0), mm)
            self._children(c, x, L, ro - t)
            return
        if kind == "bulkhead":
            ro = float(_fixed(_num(c, "outerradius"), inner_r))
            m = n * rho * np.pi * ro**2 * L
            self.parts.append(Part(name, kind, m, x + L / 2, m * (ro**2 / 4 + L**2 / 12), m * ro**2 / 2))
            return
        if kind == "centeringring":
            ro = float(_fixed(_num(c, "outerradius"), inner_r))
            ri = float(_fixed(_num(c, "innerradius"), mount_r))
            m = n * rho * np.pi * (ro**2 - ri**2) * L
            self.parts.append(Part(name, kind, m, x + L / 2, m * ((ro**2 + ri**2) / 4 + L**2 / 12),
                                   m * (ro**2 + ri**2) / 2))
            return
        if kind == "masscomponent":
            m = n * float(c.findtext("mass") or 0.0)
        elif kind == "parachute":
            D = float(c.findtext("diameter"))
            m = rho * np.pi * D**2 / 4 + _density(c, "linematerial") * int(c.findtext("linecount") or 0) * float(
                c.findtext("linelength") or 0.0)
            cd = _num(c, "cd")
            cd = AUTO_CHUTE_CD if isinstance(cd, tuple) or cd is None else cd
            alt = c.findtext("deployaltitude")
            self.chutes.append({"name": name, "cd": cd, "diameter": D, "event": (c.findtext("deployevent") or "").lower(),
                                "altitude": float(alt) if alt else None})
        elif kind == "shockcord":
            m = _density(c) * float(c.findtext("cordlength") or 0.0)
        elif kind == "streamer":
            m = rho * float(c.findtext("striplength") or 0.0) * float(c.findtext("stripwidth") or 0.0)
        elif kind == "railbutton":
            od = float(c.findtext("outerdiameter") or 0.0)
            m = n * rho * np.pi * od**2 / 4 * float(c.findtext("height") or 0.0)
        else:
            self.warnings.append(f"skipped {kind} '{name}'")
            return
        self.parts.append(Part(name, kind, m, x + L / 2, m * L**2 / 12))

    def _fins(self, c, x):
        if hasattr(self, "fins"):
            raise ValueError("more than one fin set; only one is supported")
        n = int(c.findtext("fincount"))
        Cr, Ct = float(c.findtext("rootchord")), float(c.findtext("tipchord"))
        s, xt, t = float(c.findtext("height")), float(c.findtext("sweeplength")), float(c.findtext("thickness"))
        tab = float(c.findtext("tabheight") or 0.0) * float(c.findtext("tablength") or 0.0)
        m = _density(c) * n * t * (s * (Cr + Ct) / 2 + tab)
        self.fins = {"count": n, "root_chord": Cr, "tip_chord": Ct, "span": s, "sweep": xt, "thickness": t, "station": x}
        x_bar = (Cr**2 + Cr * Ct + Ct**2 + xt * (Cr + 2 * Ct)) / (3 * (Cr + Ct))
        self.parts.append(Part(c.findtext("name") or "fins", "fins", m, x + x_bar))

    # -- overrides apply to a component and, optionally, everything inside it
    def _override(self, c, start):
        om = c.findtext("overridemass")
        if om is None:
            return
        target = float(om)
        if (c.findtext("overridesubcomponentsmass") or "false").lower() == "true":
            subtree = self.parts[start:]
            total = sum(p.mass for p in subtree)
            k = target / total if total > 0 else 0.0
            for p in subtree:
                p.mass, p.inertia, p.roll_inertia = p.mass * k, p.inertia * k, p.roll_inertia * k
        else:
            own = next((p for p in self.parts[start:]), None)
            if own is not None and own.mass > 0:
                k = target / own.mass
                own.mass, own.inertia, own.roll_inertia = target, own.inertia * k, own.roll_inertia * k
            elif own is not None:
                own.mass = target
        if c.findtext("overridecg") is not None:
            self.warnings.append(f"CG override on '{c.findtext('name')}' ignored")

    # -- assemble the TOML-shaped config
    def config(self) -> dict:
        if not hasattr(self, "fins"):
            raise ValueError("no trapezoidal fin set found")
        nose = next(p for p in self.parts if p.kind == "nose")
        fins = next(p for p in self.parts if p.kind == "fins")
        masses, used = {}, set()
        for p in self.parts:
            if p.kind in ("nose", "fins") or p.mass <= 0:
                continue
            key = re.sub(r"[^a-z0-9]+", "_", p.name.lower()).strip("_") or p.kind
            k, i = key, 2
            while k in used:
                k, i = f"{key}_{i}", i + 1
            used.add(k)
            masses[k] = {"mass": p.mass, "station": p.cg, "inertia": p.inertia, "roll_inertia": p.roll_inertia}
        body_len = self.length - self.nose["length"]
        f = self.fins
        nose_cfg = {"shape": self.nose["shape"], "length": self.nose["length"], "mass": nose.mass}
        if self.nose["param"] is not None:
            nose_cfg["param"] = self.nose["param"]
        finish = max(self.tube_spans, key=lambda s: s[1])[2] if self.tube_spans else "normal"
        cfg = {
            "name": self.rocket.findtext("name") or "imported rocket",
            "roughness": FINISH_ROUGHNESS.get(finish.lower(), 60e-6),
            "nose": nose_cfg,
            "body": {"diameter": 2 * self.radius, "length": body_len, "mass": 0.0},
            "fins": {"count": f["count"], "root_chord": f["root_chord"], "tip_chord": f["tip_chord"],
                     "span": f["span"], "sweep": f["sweep"], "thickness": f["thickness"], "mass": fins.mass,
                     "aft_offset": self.length - f["root_chord"] - f["station"]},
            "mass": masses,
            "motor": self._motor(),
            "recovery": self._recovery(),
        }
        return cfg

    def _motor(self):
        if self.mount is None:
            raise ValueError("no motor mount found")
        aft, overhang, mm = self.mount
        motors = mm.findall("motor")
        pick = [m for m in motors if m.get("configid") == self.config_id] or motors
        out = {"aft_station": aft + overhang}
        if pick:
            m = pick[0]
            out.update(designation=m.findtext("designation") or "", manufacturer=m.findtext("manufacturer") or "")
        return out

    def _recovery(self):
        chutes = sorted(self.chutes, key=lambda c: c["diameter"])
        if not chutes:
            self.warnings.append("no parachutes found; using hprsim's default recovery")
            return {}

        def cda(c):
            return c["cd"] * np.pi * c["diameter"] ** 2 / 4
        main = chutes[-1]
        rec = {"main_cda": cda(main), "main_diameter": main["diameter"]}
        if len(chutes) >= 2 and main["event"] == "altitude" and main["altitude"] is not None:
            drogue = chutes[-2]
            rec.update(drogue_cda=cda(drogue), drogue_diameter=drogue["diameter"], main_altitude=main["altitude"])
        else:
            rec.update(drogue_cda=None, main_altitude=None)
            if main["event"] == "altitude":
                self.warnings.append("a single parachute deploying at an altitude is modeled as opening at apogee")
        return rec


def read(path) -> dict:
    """Parse an .ork file into hprsim's rocket config (the same structure as a rocket TOML file)."""
    r = _Reader(_parse(path)).read()
    cfg = r.config()
    cfg["_warnings"] = r.warnings
    return cfg


def to_toml(cfg: dict) -> str:
    """Write a config from read() as a rocket TOML file."""
    def val(v):
        if v is None:
            return None
        if isinstance(v, str):
            return '"' + v.replace('"', '\\"') + '"'
        if isinstance(v, bool):
            return "true" if v else "false"
        return repr(float(v)) if isinstance(v, (float, np.floating)) else str(v)

    lines = [f"# Imported from OpenRocket: {cfg['name']}", f"name = {val(cfg['name'])}", f"roughness = {val(cfg['roughness'])}", ""]
    for table in ("nose", "body", "fins", "motor", "recovery"):
        lines.append(f"[{table}]")
        for k, v in cfg[table].items():
            if val(v) is not None:
                lines.append(f"{k} = {val(v)}")
        lines.append("")
    for name, m in cfg["mass"].items():
        lines.append(f"[mass.{name}]")
        lines += [f"{k} = {val(v)}" for k, v in m.items()]
        lines.append("")
    return "\n".join(lines)
