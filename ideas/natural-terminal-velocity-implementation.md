# Natural Terminal Velocity Formulation and Dynamic Solver Architecture

## 1. Overview and Purpose

The original solver suite in `freeFallingFoam` uses a prescribed reference frame acceleration:

```text
V_frame(t) = u0 + acceleration * t
```

In the prescribed formulation, the `terminalVelocity` parameter acts as an artificial numerical clamp. The body acceleration remains fixed independently of the actual aerodynamic drag forces calculated by OpenFOAM.

The new `trueFreeFallingPimpleFoam` solver removes prescribed accelerations and velocity clamping. The motion of the falling body emerges naturally from the dynamic force balance between gravity and CFD-calculated aerodynamic stresses (pressure and viscous forces):

```text
m * a = m * g + F_aero

a = g + F_aero / m
```

Velocity is integrated dynamically at each time step:

```text
V_frame^(n+1) = V_frame^n + a_frame^n * dt
```

The resulting frame velocity is used directly by the reference frame and the `fallingInletVelocity` boundary condition.

---

## 2. Terminology and Naming Reorganization

To prevent confusion between numerical integration schemes and physics formulations, the repository solvers are organized into two distinct categories:

1. **`legacyPrescribedMotion/`**:
   Solvers where frame acceleration is directly supplied by the user configuration dictionary (`u0` and constant `acceleration`).

   - `solvers/legacyPrescribedMotion/freeFallingPimpleFoam`
   - `solvers/legacyPrescribedMotion/freeFallingRhoPimpleFoam`
   - `solvers/legacyPrescribedMotion/freeFallingSonicFoam`
2. **`trueFreeFallingPimpleFoam/`**:
   Solvers where frame acceleration and velocity are dynamically coupled to the CFD aerodynamic forces at every time step.

---

## 3. Natural Terminal Velocity Concept

Terminal velocity is an emergent physical steady state, not an input parameter.

Terminal velocity occurs naturally when the aerodynamic drag balances gravity:

```text
F_aero + m * g = 0
```

When this balance is reached, net acceleration drops to zero (`a = 0`) while frame velocity `V_frame` remains constant. The solver never imposes or clamps terminal velocity.

---

## 4. Ground Interrupt and Operational Modes

The solver supports ground proximity tracking independent of the atmospheric density model via the `groundInterrupt` configuration option.

### Mode 1: Height-Constrained Drop (`groundInterrupt on`)

- **Configuration**: `initialAltitude` set (e.g. 500 m), `groundAltitude` set (e.g. 0 m), `groundInterrupt on`.
- **Behavior**: Altitude is dynamically integrated as `z^(n+1) = z^n + V_z^n * dt`. When `altitude <= groundAltitude`, the simulation cleanly logs ground impact velocity and terminates (`runTime.write()`, exit time loop).
- **Use Case**: Used to measure maximum speed reached before ground impact, or to determine whether a falling object hits the ground before reaching terminal velocity.

### Mode 2: Unconstrained Terminal Velocity Search (`groundInterrupt off`)

- **Configuration**: `initialAltitude 0`, `groundAltitude -1e30`, `groundInterrupt off`.
- **Behavior**: The simulation falls unconstrained until `endTime`.
- **Use Case**: Simple pure CFD free-fall simulation to computationally find the physical terminal velocity equilibrium point (`a = 0`).

---

## 5. Configuration Dictionary (`constant/fallingFrameDict`)

Below is the complete `fallingFrameDict` with all options enabled and data types specified:

```cpp
/* constant/fallingFrameDict */
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      fallingFrameDict;
}

// Initial velocity of the body at t = 0 [m/s] (type: vector)
u0                  ( 0 0 0 );

// Gravitational / initial frame acceleration vector [m/s^2] (type: vector)
acceleration        ( 0 -9.81 0 );

// Terminal velocity clamping (optional) (type: scalar)
terminalVelocity    85;

// Boundary patch(es) representing the falling body for force integration (type: wordList)
patches             ( sphere );

// Mass of the falling body [kg] (type: scalar)
mass                0.033;

// Atmospheric fluid density model: "rhoInf" (constant) or "altitudeBased" (type: word)
densityModel        rhoInf;

// Constant reference fluid density [kg/m^3] for incompressible force calculation (type: scalar)
rhoInf              1.204;

// Gravity model: "constant" or "altitudeBased" (type: word)
gravityModel        constant;

// Enable ground impact termination check (type: bool)
groundInterrupt     on;

// Starting altitude above ground level [m] (type: scalar)
initialAltitude     500;

// Ground level altitude threshold [m] (type: scalar)
groundAltitude      0;
```

### Parameter Reference Table

| Parameter            | Type         | Required / Default   | Description                                                          |
| :------------------- | :----------- | :------------------- | :------------------------------------------------------------------- |
| `u0`               | `vector`   | Mandatory            | Initial frame velocity at t = 0                                      |
| `acceleration`     | `vector`   | Mandatory            | Base gravitational acceleration vector (g)                           |
| `patches`          | `wordList` | Mandatory            | Target wall patch name(s) for pressure & viscous force integration   |
| `mass`             | `scalar`   | Mandatory            | Physical mass of the falling body in kg                              |
| `rhoInf`           | `scalar`   | Default: 1.204       | Fluid density for kinematic pressure conversion                      |
| `densityModel`     | `word`     | Default:`rhoInf`   | Fluid density formulation (`rhoInf` or `altitudeBased`)          |
| `gravityModel`     | `word`     | Default:`constant` | Gravitational acceleration model (`constant` or `altitudeBased`) |
| `groundInterrupt`  | `bool`     | Default:`off`      | Terminate simulation when`altitude <= groundAltitude`              |
| `initialAltitude`  | `scalar`   | Default: 0           | Initial physical altitude z0 in meters                               |
| `groundAltitude`   | `scalar`   | Default: 0           | Ground reference altitude zg in meters                               |
| `terminalVelocity` | `scalar`   | Optional             | Deprecated parameter; retained for legacy dictionary compatibility   |

---

## 6. Implementation Workflow in Solver Code

### Timestep Control Flow (`trueFreeFallingPimpleFoam`)

```text
At each time step n:

1. Evaluate current body altitude:
   z^(n) = z_initial + integral(V_z dt)

2. Check ground termination:
   if (groundInterrupt && altitude <= groundAltitude)
   {
       Info << "Ground impact reached at altitude " << altitude << " m" << endl;
       Info << "Impact velocity: " << mag(frameVelocity) << " m/s" << endl;
       runTime.write();
       break;
   }

3. Evaluate density & gravity models:
   rho_fluid = evaluateDensityModel(altitude)
   g_vec     = evaluateGravityModel(altitude)

4. Solve CFD momentum and pressure (p_rgh formulation):
   - Solve UEqn with current fallingFrameForce = -a_frame^(n)
   - Solve p_rghEqn, reconstruct total pressure p = p_rgh + gh
   - Correct velocity field U and face flux phi

5. Calculate aerodynamic forces on body patch(es):
   F_pressure = sum( p_rgh * n * A_face ) * rho_fluid
   F_viscous  = sum( tau * n * A_face ) * rho_fluid
   F_aero     = F_pressure + F_viscous

6. Compute net acceleration and update frame velocity:
   F_net = mass * g_vec + F_aero
   a_frame^(n+1) = F_net / mass
   V_frame^(n+1) = V_frame^(n) + a_frame^(n+1) * dt

7. Update object registry:
   Store V_frame^(n+1) and a_frame^(n+1) for fallingInletVelocity BC lookup
   fallingFrameForce = -a_frame^(n+1)
```

### Aerodynamic Force Calculation (p_rgh vs p)

In an accelerating reference frame, total pressure p contains a fictitious hydrostatic gradient gh = (fallingFrameForce . r).

Integrating total pressure over the body surface introduces an artificial buoyancy force proportional to the frame acceleration. To ensure physical accuracy:

F_aero is integrated EXCLUSIVELY from p_rgh (hydrostatic-reduced pressure) and viscous shear stresses.

This isolates pure aerodynamic drag and lift forces from reference frame artifacts.

---

## 7. Boundary Condition (fallingInletVelocity) Integration

The fallingInletVelocity patch condition is updated to query the solver object registry:

```cpp
if (db().foundObject<uniformDimensionedVectorField>("frameVelocity"))
{
    // Dynamic True Free-Fall Mode:
    const vector& Vframe = db().lookupObject<uniformDimensionedVectorField>("frameVelocity").value();
    Uinlet = -Vframe;
}
else
{
    // Legacy Prescribed Acceleration Fallback:
    Uinlet = -(u0 + acceleration * t);
}
```

This guarantees full backward compatibility with legacy prescribed-motion cases while seamlessly functioning in dynamic true free-fall mode.
