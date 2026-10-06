# Derivations

These are the equations the simulator uses, where each one comes from, and the test that checks it.

## 1. Atmosphere (`atmosphere.py`)

The International Standard Atmosphere. Below 11 km, temperature falls linearly, $T = T_0 - Lh$ with $L = 6.5$ K/km. Combining hydrostatic balance ($dP/dh = -\rho g$) with the ideal gas law gives

$$P = P_0\left(\frac{T}{T_0}\right)^{g_0/(R L)}$$

Above 11 km the temperature is constant, so pressure decays exponentially instead. Density comes from $\rho = P/(RT)$, the speed of sound from $a = \sqrt{\gamma R T}$, and viscosity from Sutherland's law.

Checked against the US Standard Atmosphere tables at 0, 5 and 11 km (`test_isa_table_values`).

## 2. Motors (`motor.py`)

Thrust comes from the certified thrust curve (a RASP `.eng` file). The mass needs one assumption: propellant burns in proportion to the impulse delivered so far,

$$m_p(t) = m_{p,0}\left(1 - \frac{I(t)}{I_\text{total}}\right) \quad\Rightarrow\quad \dot m = \frac{m_{p,0}}{I_\text{total}}F(t) = \frac{F}{c}, \qquad c = \frac{I_\text{total}}{m_{p,0}}$$

That makes the effective exhaust velocity $c$ constant, which is the same assumption OpenRocket and RockSim make.

Checked by: the mass flow integrates to exactly the propellant loaded, and the burnout speed with no gravity or drag matches the rocket equation, $\Delta v = c\ln(m_0/m_f)$ (`test_burnout_speed_matches_rocket_equation`).

**Thrust at altitude.** The curve is measured on a static stand near sea level. Thrust is $F = \dot m v_e + (p_e - p_a)A_e$, and only the last term depends on where the rocket is. Relative to the test stand,

$$F(p_a) = F_\text{curve} + (p_\text{ref} - p_a)A_e$$

with $p_\text{ref}$ = 101.325 kPa. The nozzle exit area $A_e$ is an input, because a .eng file doesn't include it. The mass flow is set by the chamber, so it doesn't change.

Checked by: in a drag-free flight at constant mass, the extra burnout speed equals $(A_e/m)\int(p_\text{ref} - p_a)\,dt$ along the trajectory (`test_pressure_thrust_adds_the_expected_speed`).

## 3. Mass properties (`rocket.py`)

Stations are measured from the nose tip, increasing toward the tail. Each part is modeled as follows:

| Part | Model |
|---|---|
| nose cone | a thin shell of uniform thickness, so its mass follows the surface area $dA = 2\pi r\,ds$. The CG and pitch inertia are integrated numerically along the real profile, treating each slice as a ring with inertia $\tfrac12 r^2\,dm$ about its own diameter. |
| body tube | a thin-walled tube, $I = m\left(\tfrac{R^2}{2} + \tfrac{L^2}{12}\right)$ |
| fins | the centroid of the trapezoidal planform, $\bar x = \dfrac{C_r^2 + C_rC_t + C_t^2 + X_t(C_r + 2C_t)}{3(C_r + C_t)}$ |
| propellant | a solid cylinder; the motor case is a thin tube |

The total pitch inertia about the moving CG is $I = \sum I_{cm,i} + \sum m_i x_i^2 - M x_{cg}^2$ (the parallel-axis theorem, summed over all parts).

## 4. Stability: Barrowman's equations (`rocket.py`)

**Nose cone.** Slender-body theory gives a normal-force slope of exactly $C_{N\alpha} = 2$ per radian, whatever the shape. The center of pressure (CP) also has a general formula,

$$X_{CP} = L - \frac{V}{A_\text{base}}$$

where $V$ is the nose volume.
- **Cone:** $V = AL/3$, so $X_{CP} = \tfrac23 L$, which matches Barrowman's table.
- **OpenRocket's parabolic series** ($r = R(2u - u^2)$): $V = \tfrac{8}{15}\pi R^2 L$, so $X_{CP} = \tfrac{7}{15}L$.

**Fins**, with $n$ fins, semi-span $s$, root chord $C_r$, tip chord $C_t$, and $L_F$ the length of the mid-chord line:

$$C_{N\alpha,f} = \left(1 + \frac{R}{s + R}\right)\frac{4n(s/d)^2}{1 + \sqrt{1 + \left(\dfrac{2L_F}{C_r + C_t}\right)^2}}$$

$$X_f = X_b + \frac{X_t}{3}\frac{C_r + 2C_t}{C_r + C_t} + \frac16\left[(C_r + C_t) - \frac{C_rC_t}{C_r + C_t}\right]$$

The first factor in $C_{N\alpha,f}$ accounts for the body increasing the flow over the fins. For a rectangular, unswept fin, $X_f$ lands exactly at the quarter chord, which is the classic thin-airfoil result (`test_rectangular_fin_cp_is_quarter_chord...`).

**Body lift at larger angles.** The terms above are linear in $\alpha$ and leave out the body tube's own lift, which matters once the angle of attack is more than a few degrees. Galejs' correction (the one OpenRocket uses) adds

$$C_{N,\text{body}} = K\,\frac{A_\text{plan}}{A_\text{ref}}\sin^2\alpha, \qquad K = 1.1$$

acting at the centroid of the side-view (planform) area of the nose and body. The nose planform is $2\int r\,dx$, with its centroid from $\int 2rx\,dx$. The CP at a finite angle is the average of the linear CP and the planform centroid, weighted by their normal forces, $C_{N\alpha}\alpha$ and $C_{N,\text{body}}$ (`Rocket.cp_at`).

Checked by: a cone's planform is a triangle, with area $RL$ and centroid at $\tfrac23 L$ (`test_cone_planform_area_and_centroid`). The force and moment the flight code applies match the formula (`test_body_lift_force_and_moment_follow_galejs`).

**Compressibility.** Speed enters the fin slope through the effective aspect ratio (Diederich's semi-empirical form, which OpenRocket also uses):

$$C_{N\alpha,f} = \left(1 + \frac{R}{s + R}\right)\frac{4n(s/d)^2}{1 + \sqrt{1 + \left(\beta\dfrac{2L_F}{C_r + C_t}\right)^2}}, \qquad \beta = \sqrt{1 - M^2}$$

This is capped at $M = 0.8$. It grows much less with speed than dividing the whole slope by $\beta$ (Prandtl–Glauert), which is what I did at first. That overstated the fins' lift above Mach 0.3, and the OpenRocket comparison in section 10 is what showed it (`test_fin_slope_compressibility_acts_through_the_aspect_ratio`).

**Total CP and static margin:**

$$X_{CP} = \frac{\sum C_{N\alpha,i}X_i}{\sum C_{N\alpha,i}}, \qquad \text{static margin} = \frac{X_{CP} - X_{CG}}{d}\ \text{(calibers)}$$

## 5. Drag (`aero.py`)

This is the component buildup from Niskanen's thesis (the basis of OpenRocket).

**Skin friction.** The turbulent friction coefficient $C_f = 1/(1.50\ln Re - 5.6)^2$ is compared with a roughness-limited value $C_{f,r} = 0.032\,(R_s/L)^{0.2}$, and the larger one is used. A compressibility factor $(1 - 0.1M^2)$ is applied.

**Friction drag.** Each surface's wetted area gets a form factor:

$$C_{D,f} = C_f\,\frac{\left(1 + \frac{1}{2f_B}\right)S_\text{body} + \left(1 + \frac{2t}{\bar c}\right)S_\text{fins}}{A_\text{ref}}$$

where $f_B$ is the body length-to-diameter ratio, $t$ the fin thickness and $\bar c$ the mean fin chord.

**Base drag** is $C_{D,b} = 0.12 + 0.13M^2$ below Mach 1 and $0.25/M$ above (they meet at 0.25). It acts on the base area. While the motor is firing, the motor's own cross-section is subtracted, because the exhaust fills that part of the base. OpenRocket doesn't do this, so the two differ while the motor burns.

**Skin friction through Mach 1.** Above Mach 1.1 the compressibility correction is $C_f/(1 + 0.15M^2)^{0.58}$ for turbulent flow and $C_{f,r}/(1 + 0.18M^2)$ when roughness-limited, instead of $C_f(1 - 0.1M^2)$. The subsonic and supersonic forms don't meet at Mach 1 (switching there dropped my drag 3% in one step), so from Mach 0.9 to 1.1 the correction is a straight-line blend of the two, as in OpenRocket. For the roughness-limited value it blends their values at Mach 0.9 and 1.1 (`test_skin_friction_blends_through_mach_one`).

**Fin edges** depend on the cross-section. A square leading edge feels the stagnation pressure, $0.85\,q_\text{stag}/q$ with $q_\text{stag}/q = 1 + M^2/4 + M^4/40$ below Mach 1 and $1.84 - 0.76/M^2 + 0.166/M^4 + 0.035/M^6$ above. A rounded one gets $(1 - M^2)^{-0.417} - 1$ below Mach 0.9, $1 - 1.785(M - 0.9)$ up to Mach 1, and $1.214 - 0.502/M^2 + 0.1095/M^4$ above. Either is scaled by $\cos^2\Lambda_\text{LE}$ times the edge area. The trailing edge is a small base: all of it for a square edge, half for a rounded one (my assumption), and none for an airfoil's sharp edge. The pieces of each formula meet where they change (`test_drag_pieces_join_up_across_mach_one`), and the square-edge case is checked by hand (`test_square_fin_edges_by_hand`).

**Not included:** the nose cone's own pressure (wave) drag. It is small below about Mach 0.8 for smooth noses, but on OpenRocket's example, leaving it out puts my drag 1.4% low at Mach 1.00 and 6.8% low by Mach 1.05 (README, "Against OpenRocket").

The example rocket comes out at $C_D \approx 0.5$ around Mach 0.3, which is in the normal range for high-power rockets (`test_example_drag_is_in_the_usual_range`). For a real rocket you then calibrate it against a flight (section 9).

## 6. Six-degree-of-freedom ascent (`flight.py`)

**State and frames.** The state is position and velocity in the world frame (east, north, up), an attitude quaternion $q$ (body to world), and the body angular rate $\omega$. The body $x$ axis points out the nose.

**Forces.** With $\vec v_a = \vec v - \vec w(z)$ the velocity relative to the air:
- **thrust:** $F\,\hat x_b$
- **gravity:** $-mg\,\hat z$
- **drag:** $-\tfrac12\rho V^2 A C_D\,\hat v_a$
- **normal force:** magnitude $\tfrac12\rho V^2 A C_{N\alpha}\,\alpha$, where the angle of attack is $\alpha = \arctan(|u_\perp|/u_x)$ and $u = R^\top\vec v_a$ is the airspeed in the body frame. It acts at the CP and points against the sideways component of the airspeed, $-\hat u_\perp$. Body lift, $\tfrac12\rho V^2 A\,C_{N,\text{body}}$, acts in the same direction at the planform centroid.

**Why it weathercocks.** The CP sits behind the CG, so the normal force at the CP creates a moment $(x_{cg} - x_{cp})\hat x_b \times \vec N$. That moment turns the nose toward the direction the rocket is moving relative to the air. In a crosswind, that means turning upwind (`test_rocket_weathercocks_into_the_wind`).

**Pitch damping.** Each part of the rocket that swings sideways meets the air, which resists the rotation. The jet leaving the nozzle carries angular momentum away as well. Together they give a damping moment $-C\,\omega_\perp$, where

$$C = \frac{\rho V A}{2}\sum_i C_{N\alpha,i}(x_i - x_{cg})^2 + \dot m\,(x_\text{nozzle} - x_{cg})^2$$

**Equations of motion:**

$$m\dot{\vec v} = \sum\vec F, \qquad I\dot\omega = M - \omega\times I\omega, \qquad \dot q = \tfrac12\,q\otimes(0, \omega)$$

**Roll.** Canted fins drive it, and the fins themselves resist it. With fin cant $\delta$, each fin makes normal force $C_{N\alpha 1}\delta$ at its mean aerodynamic chord, a distance $y_\text{MAC} + r_t$ from the axis (Barrowman). Rolling at rate $p$, a strip of fin at radius $r$ meets the air at an extra angle $pr/V$ and pushes back with the 2-D slope $2\pi/\beta$ (strip theory). Together,

$$M_\text{roll} = qA\,N C_{N\alpha 1}\delta\,(y_\text{MAC} + r_t) - q\frac{p}{V}N\frac{2\pi}{\beta}\int_{r_t}^{r_t+s} c(r)\,r^2\,dr$$

with $C_{N\alpha 1}$ the slope of one fin (the fin formula above, without the body interference factor) and $y_\text{MAC} = \tfrac{s}{3}\tfrac{C_r + 2C_t}{C_r + C_t}$. The roll rate settles where the two terms balance, which is proportional to speed, and it gets there with time constant $I_\text{roll}V/(qK_d)$. Here $K_d$ is everything multiplying $qp/V$ in the damping term (`test_canted_fins_spin_up_to_the_steady_roll_rate`).

**Thrust misalignment.** A thrust line tilted by $\varepsilon$ at the nozzle adds a side force $F\sin\varepsilon$ and a moment $(x_{cg} - x_\text{nozzle})F\sin\varepsilon$ about the CG (`test_thrust_misalignment_torques_the_rocket_about_the_cg`).

Integration is classic 4th-order Runge–Kutta, with a 2 ms step during the burn and 10 ms during the coast. Halving both steps changes the apogee by under 0.1% (`test_time_step_is_converged`).

**Launch rail.** Until the rocket has traveled the rail length, it can only slide along the rail and cannot rotate. The exit moment is interpolated within the time step (`test_vertical_flight_without_drag_matches_kinematics` checks $v = \sqrt{2aL}$ to 0.2%).

## 7. Checking the rotational dynamics

Linearize the pitch motion while coasting at speed $V$:
- $\theta$ is the body pitch angle, $\gamma$ the flight-path angle, and $\alpha = \theta - \gamma$.
- $K = \tfrac12\rho V^2 A C_{N\alpha}(x_{cp} - x_{cg})$ is the restoring stiffness.
- $Z = \tfrac12\rho V^2 A C_{N\alpha}/(mV)$ is how fast the normal force turns the flight path.

The two equations are

$$I\ddot\theta = -K\alpha - C\dot\theta, \qquad \dot\gamma = Z\alpha$$

Eliminating $\theta$ gives a damped oscillator in $\alpha$:

$$\ddot\alpha + \left(\frac{C}{I} + Z\right)\dot\alpha + \left(\frac{K}{I} + \frac{CZ}{I}\right)\alpha = 0$$

So the angle of attack oscillates and decays at the rate

$$\sigma = \frac{C}{2I} + \frac{Z}{2}$$

The test starts the rocket coasting with a small pitch rate and compares its oscillation period and decay rate against these formulas (`test_pitch_oscillation_frequency_and_damping`).

## 8. Descent

After apogee the rocket is treated as a point mass hanging under its parachute's drag area $C_dA$: a drogue first, then the main at the set altitude. At steady descent the speed is

$$v_t = \sqrt{\frac{2mg}{\rho\,C_dA}}$$

(`test_descent_reaches_terminal_velocity_under_main`).

**Opening shock.** A canopy fills over a distance, not instantly. The model assumes its diameter grows linearly with the distance $x$ traveled since deployment, so the drag area grows as $x^2$ and is full after $L = n D_0$, where $n$ is the fill constant. Without gravity, $v\,dv/dx = dv/dt$ turns $m\,dv/dt = -\tfrac12\rho v^2 C_DA\,(x/L)^2$ into

$$\frac{dv}{v} = -\frac{\rho C_DA}{2mL^2}x^2\,dx \quad\Rightarrow\quad v = v_0 e^{-kx^3},\qquad k = \frac{\rho C_DA}{6mL^2}$$

The force $F = \tfrac12\rho v^2 C_DA(x/L)^2 \propto x^2 e^{-2kx^3}$ peaks at $x^* = (3k)^{-1/3}$. If $x^* > L$ (a heavy payload), it peaks at full inflation instead. The test fills a canopy at 30 m/s for a light and a very heavy payload and checks the peak force against this formula (`test_opening_shock_matches_the_closed_form`).

$n$ depends on the canopy type. It is an assumption here (8 by default), which is why the README shows 4 and 12 too.

## 9. Calibration and Monte Carlo

**Calibration.** `calibrate.fit_cd_scale` finds the drag multiplier that makes the simulated apogee match the altimeter, using a root solve (Brent's method). The test flies a rocket with 15% extra drag, logs its altitude in feet, and checks that the fit recovers 1.15.

**Monte Carlo.** Each input is sampled from a distribution. The sensitivity numbers are standardized regression coefficients:
1. Standardize every input and the apogee.
2. Fit a linear regression.
3. Each coefficient squared is roughly that input's share of the apogee variance.

## 10. Importing OpenRocket designs (`ork.py`)

An .ork file is XML (zipped by newer OpenRocket versions, gzipped by older ones). Each part's position is given relative to its parent: from the parent's top, its bottom, its middle, right after the previous part, or absolute. Masses are computed the way OpenRocket does:

| Part | Mass |
|---|---|
| nose cone | shell area × wall thickness × density |
| tubes, couplers, inner tubes | annulus area × length × density |
| bulkheads, centering rings | disc or ring area × length × density |
| fins | (planform + tab area) × thickness × density × count |
| parachutes | cloth area × surface density, plus the shroud lines |

Mass overrides in the file are then applied. An override that covers subcomponents scales the whole subtree, which keeps its CG. Each internal part becomes a point mass at its own CG, carrying its own inertia.

Checked against a design whose masses and positions are worked out by hand (`tests/test_ork.py`). On one of OpenRocket's own example rockets, the result agrees with the dry mass, CG and CP that OpenRocket saved (README, "Against OpenRocket").

## Limitations

- No nose wave drag, so the drag falls behind above Mach 1 (7% low by Mach 1.05 on OpenRocket's example). Normal-force slopes stop changing above Mach 0.8.
- Fin cant adds no drag in the model, and roll damping comes from the fins only.
- Thrust changes with altitude only if you give the nozzle exit diameter.
- Parachute filling is a simple model: drag area grows with distance squared, and the fill constant is an assumption.
- Wind is a steady power-law profile with no gusts.
- The example rocket's dimensions are representative of 3-inch kits. They aren't a specific kit.
- The .ork import handles single-stage rockets with one fin set and no transitions.
