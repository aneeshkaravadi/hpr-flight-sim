# Motor data

| File | Motor | Source |
|---|---|---|
| `AeroTech_H128W.eng` | AeroTech H128W (29 mm) | TMT certification data 1997, via [thrustcurve.org](https://www.thrustcurve.org) |
| `AeroTech_I284W.eng` | AeroTech I284W (38 mm) | TMT certification data 1999, via thrustcurve.org |
| `AeroTech_J420R.eng` | AeroTech J420R (38 mm) | certification data, via thrustcurve.org |

All three are listed on thrustcurve.org as public domain. `thrustcurve_metadata.json` has the published totals (impulse, burn time, masses) for the tests.

The `.eng` files are 25-30 point samples of the full test-stand curves, so integrating them gives a total impulse 0.3-2.7% different from the published value. The tests allow 3%. I originally wanted the J350W, but its curve on thrustcurve.org isn't marked public domain, so I used the J420R (same 38 mm case, similar impulse) instead.

To add a motor, download its RASP `.eng` file from thrustcurve.org into this folder.
