# Flight data

Altimeter logs for calibrating the drag model. `examples/compare_flight.py` reads any CSV with a time column and an altitude column (it finds headers containing "time" and "alt"). It also accepts OpenRocket's CSV export, and converts feet to meters if the header says "ft".

```bash
python examples/compare_flight.py rockets/my_rocket.toml data/motors/AeroTech_I284W.eng data/flights/my_flight.csv \
    --wind 4 --wind-from 200 --elevation 180
```

It prints the simulated vs measured apogee and the drag multiplier that makes them agree, and saves an overlay plot next to the log.

The honest test of the calibration is to fit the drag on one flight, then predict a different flight of the same rocket on another motor.
