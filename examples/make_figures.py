"""Regenerate the figures and numbers in the README.

    python examples/make_figures.py                       # everything, 2 workers
    python examples/make_figures.py --only montecarlo     # just one section
    python examples/make_figures.py --workers 4 --mc 300

Every simulation goes through hprsim.runner: one progress bar for the whole job,
lowest CPU priority, and + / - keys to change how many run at once while it goes.
Sections: thrust, baseline, weathercock, landing, montecarlo, altitude, recovery.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from hprsim import flight, motor, rocket, runner

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "docs" / "figures"
RESULTS = ROOT / "docs" / "results.json"
ROCKET = ROOT / "rockets" / "example_3in.toml"
MOTORS = {name: ROOT / "data" / "motors" / f"AeroTech_{name}.eng" for name in ("H128W", "I284W", "J420R")}
SECTIONS = ("thrust", "baseline", "weathercock", "landing", "montecarlo", "altitude", "recovery")
SITES = (0.0, 500.0, 1000.0, 1500.0, 2000.0)  # launch site elevation, m
EXITS = (0.0, 0.012, 0.016, 0.020)  # nozzle exit diameter, m (0 = curve as certified)
WINDS = np.arange(0, 10.1, 1.0)
RAILS = {"4 ft": 1.22, "6 ft": 1.83, "8 ft": 2.44}
plt.rcParams.update({"figure.dpi": 140, "axes.grid": True, "grid.alpha": 0.3, "axes.spines.top": False,
                     "axes.spines.right": False, "font.size": 10})


def fly(spec: dict) -> dict:
    """One simulation described by a plain dict (so it can be sent to a worker process)."""
    mot = motor.Motor.from_eng(MOTORS[spec["motor"]])
    mot.impulse_scale = spec.get("impulse_scale", 1.0)
    mot.exit_area = np.pi * spec.get("exit_diameter", 0.0) ** 2 / 4
    r = rocket.load(ROCKET, motor=mot, cd_scale=spec.get("cd_scale", 1.0),
                    dry_mass_scale=spec.get("dry_mass_scale", 1.0))
    fl = flight.simulate(r, flight.Launch(**spec.get("launch", {})), flight.Recovery(**spec.get("recovery", {})),
                         ascent_only=spec.get("ascent_only", False))
    out = {"events": fl.events}
    if spec.get("keep_trajectory"):
        out.update(t=fl.t, pos=fl.pos, vel=fl.vel, margin=fl.margin)
    return out


def save(fig, name):
    fig.tight_layout()
    fig.savefig(FIG / name)
    plt.close(fig)


def monte_carlo_specs(n: int, cd_sigma: float, seed: int):
    rng = np.random.default_rng(seed)
    specs, inputs = [], []
    for _ in range(n):
        x = {"drag": rng.normal(1.0, cd_sigma), "dry mass": rng.normal(1.0, 0.02),
             "motor impulse": rng.normal(1.0, 0.03), "wind speed": rng.uniform(0.0, 8.0),
             "air temperature": rng.normal(0.0, 8.0)}
        specs.append({"motor": "I284W", "cd_scale": x["drag"], "dry_mass_scale": x["dry mass"],
                      "impulse_scale": x["motor impulse"], "ascent_only": True,
                      "launch": {"wind_speed": x["wind speed"], "wind_from_deg": rng.uniform(0, 360),
                                 "dT": x["air temperature"]}})
        inputs.append(x)
    return specs, inputs


def sensitivity(inputs, y):
    """Standardized regression coefficients; their squares approximate each input's share of the variance."""
    keys = list(inputs[0])
    X = np.array([[x[k] for k in keys] for x in inputs])
    Xs = (X - X.mean(0)) / X.std(0)
    ys = (y - y.mean()) / y.std()
    beta, *_ = np.linalg.lstsq(np.c_[np.ones(len(ys)), Xs], ys, rcond=None)
    return dict(zip(keys, beta[1:]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=",".join(SECTIONS), help="comma-separated sections")
    ap.add_argument("--workers", type=int, default=2, help="simulations at once (change live with + / -)")
    ap.add_argument("--mc", type=int, default=300, help="Monte Carlo flights per case")
    args = ap.parse_args()
    only = set(args.only.split(","))
    FIG.mkdir(parents=True, exist_ok=True)
    results = json.loads(RESULTS.read_text()) if RESULTS.exists() else {}

    # ---- build every simulation up front so one progress bar covers the whole job
    jobs: dict[str, list] = {}
    if "baseline" in only:
        jobs["baseline"] = [{"motor": m, "keep_trajectory": True} for m in MOTORS]
    if "weathercock" in only:
        jobs["weathercock"] = [{"motor": m, "ascent_only": True, "launch": {"wind_speed": float(w)}}
                               for m in ("H128W", "I284W") for w in WINDS]
        jobs["rails"] = [{"motor": "H128W", "ascent_only": True, "launch": {"rail_length": L, "wind_speed": 8.0}}
                         for L in RAILS.values()]
    if "landing" in only:
        jobs["landing"] = [{"motor": "I284W", "launch": {"wind_speed": float(w)}, "recovery": rec}
                           for w in WINDS for rec in ({}, {"drogue_cda": None, "main_altitude": None})]
    if "montecarlo" in only:
        jobs["mc_raw"], mc_raw_in = monte_carlo_specs(args.mc, 0.10, seed=1)
        jobs["mc_cal"], mc_cal_in = monte_carlo_specs(args.mc, 0.02, seed=2)
    if "recovery" in only:
        jobs["recovery"] = [{"motor": "I284W"}, {"motor": "I284W", "recovery": {"drogue_cda": None, "main_altitude": None}}]
    if "altitude" in only:
        jobs["altitude"] = [{"motor": "I284W", "ascent_only": True, "exit_diameter": d,
                             "launch": {"site_elevation": h}} for h in SITES for d in EXITS]
    flat = [s for v in jobs.values() for s in v]
    out = runner.run(fly, flat, workers=args.workers) if flat else []
    if any(o is None for o in out):
        print("stopped early: figures not updated")
        return
    res, i = {}, 0
    for k, v in jobs.items():
        res[k], i = out[i:i + len(v)], i + len(v)

    # ---- 1. thrust curves (no simulation)
    if "thrust" in only:
        fig, ax = plt.subplots(figsize=(6.5, 3.6))
        for name, path in MOTORS.items():
            m = motor.Motor.from_eng(path)
            ax.plot(m.t, m.F, label=f"{name}: {m.total_impulse:.0f} N s, {m.burn_time:.2f} s")
        ax.set_xlabel("time (s)")
        ax.set_ylabel("thrust (N)")
        ax.set_title("Certified thrust curves (thrustcurve.org, public domain)", fontsize=10)
        ax.legend(fontsize=8)
        save(fig, "thrust_curves.png")

    # ---- 2. baseline flights + stability
    if "baseline" in only:
        results["baseline"] = {}
        fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
        fig2, axm = plt.subplots(figsize=(6.5, 3.8))
        for name, o in zip(MOTORS, res["baseline"]):
            e = o["events"]
            r = rocket.load(ROCKET, motor=motor.Motor.from_eng(MOTORS[name]))
            results["baseline"][name] = {k: round(float(e[k]), 2) for k in
                                         ("apogee_m", "apogee_time", "max_speed", "max_mach", "rail_exit_speed",
                                          "max_accel_g", "landing_time")}
            results["baseline"][name]["margin_liftoff_cal"] = round(r.static_margin(0.0), 2)
            results["baseline"][name]["margin_burnout_cal"] = round(r.static_margin(r.motor.burn_time + 0.1), 2)
            asc = o["t"] <= e["apogee_time"]
            axes[0].plot(o["t"], o["pos"][:, 2], label=f"{name}: apogee {e['apogee_m']:.0f} m")
            axes[1].plot(o["t"][asc], np.linalg.norm(o["vel"][asc], axis=1), label=name)
            axm.plot(o["t"][asc], o["margin"][asc], label=name)
        axes[0].set_xlabel("time (s)")
        axes[0].set_ylabel("altitude above pad (m)")
        axes[1].set_xlabel("time (s)")
        axes[1].set_ylabel("speed (m/s)")
        for ax in axes:
            ax.legend(fontsize=8)
        fig.suptitle("Same 3-inch rocket, three motors (no wind, dual deploy at 150 m)", fontsize=10)
        save(fig, "flight_profiles.png")
        axm.axhline(1.0, color="k", ls=":", lw=1)
        axm.text(0.2, 1.05, "1 caliber", fontsize=8)
        axm.set_xlabel("time (s)")
        axm.set_ylabel("static margin (calibers)")
        axm.set_title("Margin rises as propellant burns (CG moves forward)\n"
                      "and with speed (fins lift more near Mach 0.7), then settles", fontsize=10)
        axm.legend(fontsize=8)
        save(fig2, "stability_margin.png")

    # ---- 3. weathercocking
    if "weathercock" in only:
        results["weathercock"] = {}
        fig, axes = plt.subplots(1, 3, figsize=(14, 3.9))
        for j, name in enumerate(("H128W", "I284W")):
            rows = [o["events"] for o in res["weathercock"][j * len(WINDS):(j + 1) * len(WINDS)]]
            ap = np.array([e["apogee_m"] for e in rows])
            aoa = np.array([e["rail_exit_aoa_deg"] for e in rows])
            drift = np.array([-e["apogee_xy"][0] for e in rows])  # wind from the west, so upwind is -x
            axes[0].plot(WINDS, 100 * (ap / ap[0] - 1), "o-", label=name)
            axes[1].plot(WINDS, aoa, "o-", label=f"{name} (leaves rail at {rows[0]['rail_exit_speed']:.0f} m/s)")
            axes[2].plot(WINDS, drift, "o-", label=name)
            results["weathercock"][name] = {"apogee_loss_pct_at_8ms": round(float(100 * (1 - ap[8] / ap[0])), 1),
                                            "rail_exit_aoa_deg_at_8ms": round(float(aoa[8]), 1),
                                            "upwind_apogee_drift_m_at_8ms": round(float(drift[8]), 0)}
        for ax, lab in zip(axes, ("apogee change (%)", "angle of attack leaving the rail (deg)",
                                  "apogee upwind of the pad (m)")):
            ax.set_xlabel("wind at 10 m (m/s)")
            ax.set_ylabel(lab)
            ax.legend(fontsize=8)
        fig.suptitle("Weathercocking: the rocket turns into the wind, loses altitude and flies upwind (6 ft rail)",
                     fontsize=10)
        save(fig, "weathercock.png")
        results["rail_length_H128W_8ms"] = {
            lab: {"exit_speed": round(o["events"]["rail_exit_speed"], 1),
                  "exit_aoa_deg": round(o["events"]["rail_exit_aoa_deg"], 1)}
            for lab, o in zip(RAILS, res["rails"])}

    # ---- 4. landing dispersion
    if "landing" in only:
        labels = ("dual deploy (drogue, main at 150 m)", "main at apogee")
        dist = {lab: [o["events"]["landing_distance_m"] for o in res["landing"][k::2]] for k, lab in enumerate(labels)}
        times = {lab: [o["events"]["landing_time"] for o in res["landing"][k::2]] for k, lab in enumerate(labels)}
        fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
        for lab in labels:
            axes[0].plot(WINDS, dist[lab], "o-", label=lab)
            axes[1].plot(WINDS, times[lab], "o-", label=lab)
        axes[0].set_ylabel("landing distance from pad (m)")
        axes[1].set_ylabel("time from liftoff to landing (s)")
        for ax in axes:
            ax.set_xlabel("wind at 10 m (m/s)")
            ax.legend(fontsize=8)
        fig.suptitle("I284W flight: why dual deploy exists", fontsize=10)
        save(fig, "landing_dispersion.png")
        results["landing_m_at_8ms"] = {lab: round(float(v[8]), 0) for lab, v in dist.items()}

    # ---- 5. Monte Carlo
    if "montecarlo" in only:
        ys = np.array([o["events"]["apogee_m"] for o in res["mc_raw"]])
        yc = np.array([o["events"]["apogee_m"] for o in res["mc_cal"]])
        src = sensitivity(mc_raw_in, ys)
        share = {k: v**2 / sum(b**2 for b in src.values()) for k, v in src.items()}
        src_c = sensitivity(mc_cal_in, yc)
        share_c = {k: v**2 / sum(b**2 for b in src_c.values()) for k, v in src_c.items()}

        def stats(y):
            return {"p5": round(float(np.percentile(y, 5))), "p50": round(float(np.median(y))),
                    "p95": round(float(np.percentile(y, 95))), "std": round(float(y.std()), 1)}

        results["monte_carlo"] = {
            "n": args.mc, "uncalibrated": stats(ys), "drag_calibrated": stats(yc),
            "variance_share_uncalibrated": {k: round(float(v), 3) for k, v in share.items()},
            "variance_share_calibrated": {k: round(float(v), 3) for k, v in share_c.items()},
            "inputs": "drag x N(1, 0.10) [calibrated: 0.02]; dry mass x N(1, 0.02); impulse x N(1, 0.03); "
                      "wind U(0, 8) m/s from a random direction; temperature ISA + N(0, 8) K"}
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        bins = np.linspace(min(ys.min(), yc.min()), max(ys.max(), yc.max()), 35)
        axes[0].hist(ys, bins=bins, alpha=0.6, label=f"drag model as-is: sigma {ys.std():.0f} m")
        axes[0].hist(yc, bins=bins, alpha=0.6, label=f"drag calibrated from one flight: sigma {yc.std():.0f} m")
        axes[0].set_xlabel("apogee (m)")
        axes[0].set_ylabel(f"count (of {args.mc})")
        axes[0].legend(fontsize=8)
        order = sorted(share, key=share.get)
        yy = np.arange(len(order))
        axes[1].barh(yy - 0.2, [100 * share[k] for k in order], 0.4, label="drag model as-is", color="C0")
        axes[1].barh(yy + 0.2, [100 * share_c[k] for k in order], 0.4, label="drag calibrated", color="C1")
        axes[1].set_yticks(yy, order)
        axes[1].set_xlabel("share of apogee variance (%)")
        axes[1].legend(fontsize=8)
        axes[1].set_title("Where the uncertainty comes from", fontsize=10)
        fig.suptitle(f"I284W, {args.mc} simulated flights per case", fontsize=10)
        save(fig, "monte_carlo.png")

    # ---- 6. thrust gained at altitude
    if "altitude" in only:
        ap = np.array([o["events"]["apogee_m"] for o in res["altitude"]]).reshape(len(SITES), len(EXITS))
        gain = 100 * (ap / ap[:, :1] - 1)  # vs. the as-certified curve at the same site
        fig, ax = plt.subplots(figsize=(6.5, 3.8))
        for j, d in enumerate(EXITS[1:], start=1):
            ax.plot(SITES, gain[:, j], "o-", label=f"nozzle exit {1000 * d:.0f} mm")
        ax.set_xlabel("launch site elevation (m)")
        ax.set_ylabel("apogee gain from pressure thrust (%)")
        ax.set_title("I284W: thrust curves are measured near sea level;\n"
                     "higher up the nozzle pushes against thinner air", fontsize=10)
        ax.legend(fontsize=8)
        save(fig, "altitude_thrust.png")
        results["altitude_thrust"] = {
            "apogee_m_as_certified": {f"{h:.0f}": round(float(a), 1) for h, a in zip(SITES, ap[:, 0])},
            "apogee_gain_pct": {f"{1000 * d:.0f}mm": {f"{h:.0f}": round(float(g), 2) for h, g in zip(SITES, gain[:, j])}
                                for j, d in enumerate(EXITS) if d > 0}}

    # ---- 7. parachute opening loads
    if "recovery" in only:
        r = rocket.load(ROCKET, motor=motor.Motor.from_eng(MOTORS["I284W"]))
        m = r.mass_props(r.motor.burn_time + 1.0)[0]
        speeds = np.arange(10.0, 90.1, 5.0)
        dual, apogee_main = res["recovery"][0]["events"], res["recovery"][1]["events"]
        fig, ax = plt.subplots(figsize=(6.5, 3.9))
        curves = {}
        for n in (4.0, 8.0, 12.0):
            rec = flight.Recovery(drogue_cda=None, main_altitude=None, fill_constant=n)
            loads = []
            for V in speeds:  # main fired at 150 m while falling at V
                ev = {}
                flight.descend(flight.Launch(), rec, m, 0.0, np.array([0.0, 0.0, 150.0, 0.0, 0.0, -V]), ev, t_end=8.0)
                loads.append(ev["main_opening_force_N"])
            curves[n] = np.array(loads)
            ax.plot(speeds, curves[n], "-" if n == 8.0 else "--", label=f"fills in {n:.0f} canopy diameters")
        ax.plot([dual["drogue_descent_speed"]], [dual["main_opening_force_N"]], "ko")
        ax.annotate("dual deploy:\nmain opens under the drogue", (dual["drogue_descent_speed"], dual["main_opening_force_N"]),
                    textcoords="offset points", xytext=(10, 25), fontsize=8, arrowprops={"arrowstyle": "-", "lw": 0.6})
        ax.set_xlabel("falling speed when the main opens (m/s)")
        ax.set_ylabel("peak opening force (N)")
        sec = ax.secondary_yaxis("right", functions=(lambda f: f / (m * 9.80665), lambda g: g * m * 9.80665))
        sec.set_ylabel("(g)")
        ax.set_title(f"60 in main on the {m:.1f} kg rocket: the load grows with the square of the speed", fontsize=10)
        ax.legend(fontsize=8, loc="upper left")
        save(fig, "opening_shock.png")
        results["recovery"] = {
            "mass_kg": round(float(m), 2),
            "dual_deploy": {k: round(float(dual[k]), 1) for k in ("drogue_descent_speed", "drogue_opening_force_N",
                                                                   "main_opening_force_N", "main_opening_g")},
            "main_at_apogee": {k: round(float(apogee_main[k]), 1) for k in ("main_opening_force_N", "main_opening_g")},
            "main_force_N_vs_speed": {f"fill_{n:.0f}D": {f"{v:.0f}": round(float(f), 0) for v, f in zip(speeds, c)}
                                      for n, c in curves.items()}}

    RESULTS.write_text(json.dumps(results, indent=2))
    print(f"updated: {', '.join(sorted(only))}  ->  {RESULTS.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
