"""OpenRocket import, checked against masses and positions worked out by hand for a small test design."""
import gzip
import zipfile
from pathlib import Path

import numpy as np
import pytest

from hprsim import flight, ork, rocket

ROOT = Path(__file__).resolve().parents[1]
I284W = ROOT / "data" / "motors" / "AeroTech_I284W.eng"

DESIGN = """<?xml version='1.0' encoding='utf-8'?>
<openrocket version="1.10" creator="hprsim test">
  <rocket>
    <name>Test rocket</name>
    <motorconfiguration configid="cfg1" default="true"><stage number="0" active="true"/></motorconfiguration>
    <subcomponents>
      <stage>
        <name>Sustainer</name>
        <subcomponents>
          <nosecone>
            <name>Nose cone</name>
            <material type="bulk" density="1000.0">Test</material>
            <length>0.3</length><thickness>0.002</thickness><shape>conical</shape><aftradius>auto 0.04</aftradius>
          </nosecone>
          <bodytube>
            <name>Upper tube</name>
            <material type="bulk" density="1200.0">Test</material>
            <length>0.5</length><thickness>0.001</thickness><radius>0.04</radius>
            UPPER_OVERRIDE
            <subcomponents>
              <masscomponent>
                <name>Avionics</name><position type="top">0.1</position>
                <mass>0.2</mass><packedlength>0.05</packedlength>
              </masscomponent>
              <parachute>
                <name>Main</name><position type="bottom">0.0</position>
                <material type="surface" density="0.05">Nylon</material>
                <cd>1.5</cd><diameter>1.2</diameter><deployevent>altitude</deployevent><deployaltitude>200</deployaltitude>
                <linecount>8</linecount><linelength>1.0</linelength>
                <linematerial type="line" density="0.001">Line</linematerial>
                <packedlength>0.1</packedlength>
              </parachute>
            </subcomponents>
          </bodytube>
          <bodytube>
            <name>Lower tube</name>
            <material type="bulk" density="1200.0">Test</material>
            <length>0.7</length><thickness>0.001</thickness><radius>auto 0.04</radius>
            <subcomponents>
              <innertube>
                <name>Motor tube</name><position type="bottom">0.0</position>
                <material type="bulk" density="1000.0">Test</material>
                <length>0.4</length><outerradius>0.02</outerradius><thickness>0.001</thickness>
                <motormount><overhang>0.01</overhang>
                  <motor configid="cfg1"><designation>I284W</designation><manufacturer>AeroTech</manufacturer></motor>
                </motormount>
              </innertube>
              <trapezoidfinset>
                <name>Fins</name><position type="bottom">0.0</position>
                <material type="bulk" density="1500.0">Test</material>
                <fincount>3</fincount><rootchord>0.12</rootchord><tipchord>0.06</tipchord>
                <sweeplength>0.06</sweeplength><height>0.08</height><thickness>0.003</thickness>
              </trapezoidfinset>
              <parachute>
                <name>Drogue</name><position type="top">0.05</position>
                <material type="surface" density="0.05">Nylon</material>
                <cd>auto</cd><diameter>0.3</diameter><deployevent>apogee</deployevent>
                <packedlength>0.04</packedlength>
              </parachute>
            </subcomponents>
          </bodytube>
          EXTRA
        </subcomponents>
      </stage>
    </subcomponents>
  </rocket>
</openrocket>
"""


def write(tmp_path, fmt="zip", upper_override="", extra=""):
    xml = DESIGN.replace("UPPER_OVERRIDE", upper_override).replace("EXTRA", extra).encode()
    f = tmp_path / f"test_{fmt}.ork"
    if fmt == "zip":
        with zipfile.ZipFile(f, "w") as z:
            z.writestr("rocket.ork", xml)
    elif fmt == "gzip":
        f.write_bytes(gzip.compress(xml))
    else:
        f.write_bytes(xml)
    return f


def annulus(rho, ro, t, L):
    return rho * np.pi * (ro**2 - (ro - t) ** 2) * L


def test_import_matches_hand_calculation(tmp_path):
    cfg = ork.read(write(tmp_path))
    assert cfg["_warnings"] == []
    assert cfg["body"] == {"diameter": 0.08, "length": pytest.approx(1.2), "mass": 0.0}
    # conical shell: slant area pi R sqrt(R^2 + L^2), 2 mm wall, 1000 kg/m^3
    assert cfg["nose"]["mass"] == pytest.approx(1000 * np.pi * 0.04 * np.hypot(0.04, 0.3) * 0.002, rel=1e-4)
    m = cfg["mass"]
    assert m["upper_tube"]["mass"] == pytest.approx(annulus(1200, 0.04, 0.001, 0.5))
    assert m["upper_tube"]["station"] == pytest.approx(0.3 + 0.25)
    assert m["lower_tube"]["station"] == pytest.approx(0.8 + 0.35)
    assert m["avionics"] == pytest.approx({"mass": 0.2, "station": 0.3 + 0.1 + 0.025, "inertia": 0.2 * 0.05**2 / 12,
                                           "roll_inertia": 0.0})
    canopy = 0.05 * np.pi * 1.2**2 / 4 + 0.001 * 8 * 1.0  # cloth + shroud lines
    assert m["main"]["mass"] == pytest.approx(canopy)
    assert m["main"]["station"] == pytest.approx(0.8 - 0.1 + 0.05)  # packed at the bottom of the upper tube
    assert m["motor_tube"]["mass"] == pytest.approx(annulus(1000, 0.02, 0.001, 0.4))
    assert m["motor_tube"]["station"] == pytest.approx(1.5 - 0.2)
    f = cfg["fins"]
    assert f["mass"] == pytest.approx(1500 * 3 * 0.003 * 0.08 * (0.12 + 0.06) / 2)
    assert f["aft_offset"] == pytest.approx(0.0, abs=1e-12)  # flush with the aft end
    assert f["cross_section"] == "square"  # OpenRocket's default when the file doesn't say
    assert cfg["motor"] == {"aft_station": pytest.approx(1.51), "designation": "I284W", "manufacturer": "AeroTech"}
    rec = cfg["recovery"]
    assert rec["main_cda"] == pytest.approx(1.5 * np.pi * 1.2**2 / 4)
    assert rec["drogue_cda"] == pytest.approx(ork.AUTO_CHUTE_CD * np.pi * 0.3**2 / 4)  # Cd left on auto
    assert rec["main_altitude"] == 200.0


@pytest.mark.parametrize("fmt", ["gzip", "plain"])
def test_older_file_formats_read_the_same(tmp_path, fmt):
    assert ork.read(write(tmp_path, fmt)) == ork.read(write(tmp_path, "zip"))


def test_subtree_mass_override_keeps_the_cg(tmp_path):
    plain = ork.read(write(tmp_path))
    over = ork.read(write(tmp_path, upper_override="<overridemass>1.0</overridemass>"
                                                 "<overridesubcomponentsmass>true</overridesubcomponentsmass>"))
    names = ("upper_tube", "avionics", "main")

    def total_and_cg(cfg):
        ms = [cfg["mass"][k]["mass"] for k in names]
        return sum(ms), sum(cfg["mass"][k]["mass"] * cfg["mass"][k]["station"] for k in names) / sum(ms)
    assert total_and_cg(over)[0] == pytest.approx(1.0)
    assert total_and_cg(over)[1] == pytest.approx(total_and_cg(plain)[1])
    assert over["mass"]["motor_tube"] == plain["mass"]["motor_tube"]  # outside the override


def test_imported_rocket_flies_and_round_trips_through_toml(tmp_path):
    path = write(tmp_path)
    r = rocket.load(path, motor_file=I284W)
    rec = flight.Recovery.from_file(path)
    fl = flight.simulate(r, flight.Launch(), rec)
    assert fl.apogee > 100 and fl.events["landing_time"] > fl.events["apogee_time"]
    toml = tmp_path / "imported.toml"
    toml.write_text(ork.to_toml(ork.read(path)))
    r2 = rocket.load(toml, motor_file=I284W)
    assert r2.mass_props(0.0) == pytest.approx(r.mass_props(0.0), rel=1e-5)
    assert r2.stability(0.3) == pytest.approx(r.stability(0.3), rel=1e-5)
    assert vars(flight.Recovery.from_file(toml)) == pytest.approx(vars(rec), rel=1e-9)


def test_single_parachute_at_apogee_survives_the_toml_round_trip(tmp_path):
    single = DESIGN.replace("<deployevent>altitude</deployevent><deployaltitude>200</deployaltitude>",
                            "<deployevent>apogee</deployevent>")
    single = single[:single.index("<parachute>\n                <name>Drogue")] + single[single.index("</parachute>", single.index("<name>Drogue")) + len("</parachute>"):]
    f = tmp_path / "single.ork"
    f.write_text(single.replace("UPPER_OVERRIDE", "").replace("EXTRA", ""))
    rec = flight.Recovery.from_file(f)
    assert rec.drogue_cda is None and rec.main_altitude is None
    toml = tmp_path / "single.toml"
    toml.write_text(ork.to_toml(ork.read(f)))
    back = flight.Recovery.from_file(toml)
    assert not back.drogue_cda and back.main_altitude is None and back.main_cda == pytest.approx(rec.main_cda)


def test_finds_the_motor_by_designation():
    assert rocket.find_motor("I284W", ROOT / "data" / "motors") == I284W
    with pytest.raises(FileNotFoundError, match="thrustcurve.org"):
        rocket.find_motor("Z9000", ROOT / "data" / "motors")


def test_transitions_are_rejected_clearly(tmp_path):
    with pytest.raises(ValueError, match="transitions"):
        ork.read(write(tmp_path, extra="<transition><name>Boattail</name><length>0.05</length></transition>"))
