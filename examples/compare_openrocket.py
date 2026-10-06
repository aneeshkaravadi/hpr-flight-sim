"""Check hprsim against the results OpenRocket saved inside an .ork file.

    python examples/compare_openrocket.py "Dual parachute deployment.ork" [--figure]

Compares the dry mass, the dry CG, the CP at every flight point OpenRocket
stored (at that point's Mach number and angle of attack), and the drag
coefficient while coasting (at that point's Mach and Reynolds numbers). It only
uses the design and OpenRocket's own output, so no thrust curve is needed.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from hprsim import aero, ork, rocket
from hprsim.motor import Motor

ROOT = Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ork", type=Path)
    ap.add_argument("--figure", action="store_true", help="save docs/figures/openrocket_cp.png and _drag.png")
    args = ap.parse_args()
    cfg = ork.read(args.ork)
    placeholder = Motor("none", 0.038, 0.1, 1e-9, 2e-9, [0.0, 0.1, 1.0, 1.01], [0.0, 1.0, 1.0, 0.0])
    r = rocket.from_config(cfg, args.ork.parent, motor=placeholder)
    m_dry, cg_dry = r.dry_mass_props()
    sims = ork.stored_simulations(args.ork)
    if not sims:
        raise SystemExit("no stored simulation results in this file (run the simulations in OpenRocket and save)")
    print(f"{cfg['name']}: body diameter {1000 * r.diameter:.1f} mm")
    print(f"{'motor':10s} {'dry mass (kg)':>22s} {'dry CG (m)':>20s}")
    mach, diff = [], []
    for s in sims:
        d = s["data"]
        m0, mm0, cg0 = d["Mass"][0], d["Motor mass"][0], d["CG location"][0]
        x_motor = cfg["motor"]["aft_station"] - (s["motor_length"] or 0.0) / 2  # motor CG ~ its middle
        cg_or = (cg0 * m0 - mm0 * x_motor) / (m0 - mm0)
        print(f"{s['motor'] or '?':10s} {m0 - mm0:10.4f} vs {m_dry:8.4f}   {cg_or:8.4f} vs {cg_dry:8.4f}")
        ok = ~np.isnan(d["CP location"]) & ~np.isnan(d["Angle of attack"]) & (d["Angle of attack"] < 0.02) & (d["Mach number"] < 0.8)
        for M, a, cp in zip(d["Mach number"][ok], d["Angle of attack"][ok], d["CP location"][ok]):
            mach.append(M)
            diff.append(1000 * (r.cp_at(max(a, 1e-4), M) - cp))
    mach, diff = np.array(mach), np.array(diff)
    dm, cd_or, cd_ours = [], [], []
    for s in sims:
        d = s["data"]
        coast = (d.get("Thrust", np.zeros_like(d["Time"])) <= 0) & (d["Vertical velocity"] > 0) & (d["Mach number"] > 0.05)
        coast &= ~np.isnan(d["Drag coefficient"]) & ~np.isnan(d["Reynolds number"])
        for M, re, cd in zip(d["Mach number"][coast], d["Reynolds number"][coast], d["Drag coefficient"][coast]):
            dm.append(M)
            cd_or.append(cd)
            cd_ours.append(aero.drag_coefficient(r, M, re, thrusting=False))
    dm, cd_or, cd_ours = np.array(dm), np.array(cd_or), np.array(cd_ours)
    print(f"\nCP, hprsim minus OpenRocket, at {len(diff)} flight points (angle of attack under 1.1 deg):")
    for lo, hi in ((0.0, 0.2), (0.2, 0.4), (0.4, 0.6), (0.6, 0.8)):
        k = (mach >= lo) & (mach < hi)
        if k.any():
            print(f"  Mach {lo:.1f}-{hi:.1f}: {k.sum():4d} points, mean {diff[k].mean():+.2f} mm, largest {np.abs(diff[k]).max():.2f} mm")
    print(f"  overall RMS {np.sqrt(np.mean(diff**2)):.2f} mm")
    if len(dm):
        print(f"\nDrag coefficient while coasting, hprsim / OpenRocket, at {len(dm)} flight points:")
        for lo, hi in ((0.0, 0.3), (0.3, 0.6), (0.6, 0.9), (0.9, 1.0), (1.0, 1.5)):
            k = (dm >= lo) & (dm < hi)
            if k.any():
                print(f"  Mach {lo:.1f}-{hi:.1f}: {k.sum():4d} points, {np.mean(cd_ours[k]):.3f} vs {np.mean(cd_or[k]):.3f} "
                      f"({100 * (np.mean(cd_ours[k]) / np.mean(cd_or[k]) - 1):+.1f}%)")
    if args.figure:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6.5, 3.6))
        ax.plot(mach, diff, ".", ms=3, alpha=0.5)
        ax.axhline(0.0, color="k", lw=0.8)
        ax.set_xlabel("Mach number")
        ax.set_ylabel("CP difference, hprsim - OpenRocket (mm)")
        ax.set_title(f"'{cfg['name']}': {len(diff)} flight points from OpenRocket's saved simulations\n"
                     f"(body diameter {1000 * r.diameter:.0f} mm)", fontsize=10)
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(ROOT / "docs" / "figures" / "openrocket_cp.png", dpi=140)
        if len(dm):
            fig, ax = plt.subplots(figsize=(6.5, 3.6))
            order = np.argsort(dm)
            ax.plot(dm[order], cd_or[order], ".", ms=4, label="OpenRocket")
            ax.plot(dm[order], cd_ours[order], ".", ms=4, label="hprsim")
            ax.set_xlabel("Mach number")
            ax.set_ylabel("drag coefficient (coasting)")
            ax.set_title(f"'{cfg['name']}', square-edged fins: {len(dm)} coasting points\n"
                         "hprsim leaves out the nose's wave drag, which shows above Mach 1", fontsize=10)
            ax.legend(fontsize=8)
            ax.grid(alpha=0.3)
            fig.tight_layout()
            fig.savefig(ROOT / "docs" / "figures" / "openrocket_drag.png", dpi=140)
        print("saved docs/figures/openrocket_cp.png and openrocket_drag.png")


if __name__ == "__main__":
    main()
