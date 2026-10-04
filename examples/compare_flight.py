"""Compare a real flight against the simulator and calibrate drag.

    python examples/compare_flight.py rockets/my_rocket.toml data/motors/AeroTech_I284W.eng \
        data/flights/2026-10-12_I284W.csv --wind 4 --wind-from 200 --elevation 180

Prints the uncalibrated and calibrated apogee and the fitted drag multiplier, and
saves an overlay plot next to the log file.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from hprsim import calibrate, flight, motor, rocket


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rocket", type=Path)
    ap.add_argument("motor", type=Path)
    ap.add_argument("log", type=Path)
    ap.add_argument("--wind", type=float, default=0.0, help="m/s at 10 m")
    ap.add_argument("--wind-from", type=float, default=270.0, help="degrees, where the wind comes from")
    ap.add_argument("--elevation", type=float, default=200.0, help="launch site, m above sea level")
    ap.add_argument("--rail", type=float, default=1.83, help="effective rail length, m")
    ap.add_argument("--temp-offset", type=float, default=0.0, help="K warmer than ISA")
    args = ap.parse_args()

    mot = motor.Motor.from_eng(args.motor)
    launch = flight.Launch(rail_length=args.rail, site_elevation=args.elevation, wind_speed=args.wind,
                           wind_from_deg=args.wind_from, dT=args.temp_offset)
    t_log, alt_log = calibrate.read_altitude_log(args.log)
    measured = float(alt_log.max())

    def make(scale):
        return rocket.load(args.rocket, motor=mot, cd_scale=scale)

    base = flight.simulate(make(1.0), launch)
    scale = calibrate.fit_cd_scale(make, launch, measured)
    tuned = flight.simulate(make(scale), launch)
    print(f"measured apogee   {measured:8.1f} m")
    print(f"simulated (raw)   {base.apogee:8.1f} m  ({100 * (base.apogee - measured) / measured:+.1f}%)")
    print(f"fitted cd_scale   {scale:8.3f}")

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(t_log, alt_log, "k", label=f"altimeter: {measured:.0f} m")
    ax.plot(base.t, base.pos[:, 2], "C0--", label=f"simulated, uncalibrated: {base.apogee:.0f} m")
    ax.plot(tuned.t, tuned.pos[:, 2], "C3", lw=1, label=f"simulated, drag x {scale:.2f}")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("altitude above pad (m)")
    ax.legend()
    fig.tight_layout()
    out = args.log.with_suffix(".png")
    fig.savefig(out, dpi=140)
    print(f"plot -> {out}")


if __name__ == "__main__":
    main()
