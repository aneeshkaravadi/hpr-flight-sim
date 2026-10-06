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

**Compressibility.** The fin slope is multiplied by $1/\sqrt{1 - M^2}$ (Prandtl–Glauert), capped at $M = 0.8$. That's why the static margin rises near max speed.

**Total CP and static margin:**

$$X_{CP} = \frac{\sum C_{N\alpha,i}X_i}{\sum C_{N\alpha,i}}, \qquad \text{static margin} = \frac{X_{CP} - X_{CG}}{d}\ \text{(calibers)}$$

## 5. Drag (`aero.py`)

This is the component buildup from Niskanen's thesis (the basis of OpenRocket).

**Skin friction.** The turbulent friction coefficient $C_f = 1/(1.50\ln Re - 5.6)^2$ is compared with a roughness-limited value $C_{f,r} = 0.032\,(R_s/L)^{0.2}$, and the larger one is used. A compressibility factor $(1 - 0.1M^2)$ is applied.

**Friction drag.** Each surface's wetted area gets a form factor:

$$C_{D,f} = C_f\,\frac{\left(1 + \frac{1}{2f_B}\right)S_\text{body} + \left(1 + \frac{2t}{\bar c}\right)S_\text{fins}}{A_\text{ref}}$$

where $f_B$ is the body length-to-diameter ratio, $t$ the fin thickness and $\bar c$ the mean fin chord.

**Base drag** is $C_{D,b} = 0.12 + 0.13M^2$ on the base area. While the motor is firing, the motor's own cross-section is subtracted, because the exhaust fills that part of the base.

**Fins also carry** a rounded-leading-edge term and a blunt-trailing-edge term (the trailing edge acts like a small base).

The example rocket comes out at $C_D \approx 0.5$ around Mach 0.3, which is in the normal range for high-power rockets (`test_example_drag_is_in_the_usual_range`). For a real rocket you then calibrate it against a flight (section 9).

## 6. Six-degree-of-freedom ascent (`flight.py`)

**State and frames.** The state is position and velocity in the world frame (east, north, up), an attitude quaternion $q$ (body to world), and the body angular rate $\omega$. The body $x$ axis points out the nose.

**Forces.** With $\vec v_a = \vec v - \vec w(z)$ the velocity relative to the air:
- **thrust:** $F\,\hat x_b$
- **gravity:** $-mg\,\hat z$
- **drag:** $-\tfrac12\rho V^2 A C_D\,\hat v_a$
- **normal force:** magnitude $\tfrac12\rho V^2 A C_{N\alpha}\,\alpha$, where the angle of attack is $\alpha = \arctan(|u_\perp|/u_x)$ and $u = R^\top\vec v_a$ is the airspeed in the body frame. It acts at the CP and points against the sideways component of the airspeed, $-\hat u_\perp$.

**Why it weathercocks.** The CP sits behind the CG, so the normal force at the CP creates a moment $(x_{cg} - x_{cp})\hat x_b \times \vec N$. That moment turns the nose toward the direction the rocket is moving relative to the air. In a crosswind, that means turning upwind (`test_rocket_weathercocks_into_the_wind`).

**Pitch damping.** Each part of the rocket that swings sideways meets the air, which resists the rotation. The jet leaving the nozzle carries angular momentum away as well. Together they give a damping moment $-C\,\omega_\perp$, where

$$C = \frac{\rho V A}{2}\sum_i C_{N\alpha,i}(x_i - x_{cg})^2 + \dot m\,(x_\text{nozzle} - x_{cg})^2$$

**Equations of motion:**

$$m\dot{\vec v} = \sum\vec F, \qquad I\dot\omega = M - \omega\times I\omega, \qquad \dot q = \tfrac12\,q\otimes(0, \omega)$$

Roll is not modeled. Integration is classic 4th-order Runge–Kutta, with a 2 ms step during the burn and 10 ms during the coast. Halving both steps changes the apogee by under 0.1% (`test_time_step_is_converged`).

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

(`test_descent_reaches_terminal_velocity_under_main`). Parachute inflation time is ignored.

## 9. Calibration and Monte Carlo

**Calibration.** `calibrate.fit_cd_scale` finds the drag multiplier that makes the simulated apogee match the altimeter, using a root solve (Brent's method). The test flies a rocket with 15% extra drag, logs its altitude in feet, and checks that the fit recovers 1.15.

**Monte Carlo.** Each input is sampled from a distribution. The sensitivity numbers are standardized regression coefficients:
1. Standardize every input and the apogee.
2. Fit a linear regression.
3. Each coefficient squared is roughly that input's share of the apogee variance.

## Limitations

- Subsonic only. There is no transonic drag rise, and no body lift at higher angles of attack (Galejs).
- No roll dynamics, fin cant or thrust misalignment.
- Thrust changes with altitude only if you give the nozzle exit diameter.
- Parachutes open instantly.
- Wind is a steady power-law profile with no gusts.
- The example rocket's dimensions are representative of 3-inch kits. They aren't a specific kit.
