"""Turn an OpenRocket design into a rocket TOML file you can read and edit.

    python examples/import_ork.py "My rocket.ork" rockets/my_rocket.toml

rocket.load() can also fly an .ork directly; the TOML is just easier to tweak.
"""
from __future__ import annotations

import argparse
from pathlib import Path

from hprsim import ork


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ork", type=Path)
    ap.add_argument("out", type=Path)
    args = ap.parse_args()
    cfg = ork.read(args.ork)
    args.out.write_text(ork.to_toml(cfg))
    nose, body, fins = cfg["nose"], cfg["body"], cfg["fins"]
    dry = nose["mass"] + fins["mass"] + sum(m["mass"] for m in cfg["mass"].values())
    print(f"{cfg['name']}: {1000 * body['diameter']:.1f} mm x {nose['length'] + body['length']:.3f} m, "
          f"{dry:.3f} kg without the motor, {len(cfg['mass'])} internal parts")
    print(f"motor: {cfg['motor'].get('manufacturer', '')} {cfg['motor'].get('designation', '?')} "
          f"(put its .eng file in data/motors/)")
    for w in cfg["_warnings"]:
        print("note:", w)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
