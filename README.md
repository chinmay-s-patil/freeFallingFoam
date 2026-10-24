# freeFallingFoam

Custom OpenFOAM solver suite for simulating body-fixed accelerating reference frames, free-falling projectiles, atmospheric re-entry bodies, and 6DoF tumbling dynamics.

> [!WARNING]
> **Work in Progress: RFGTumble & 6DoF Tumbling Solvers are under active development!**
> 
> - **Working & Tested**: `freeFallingPimpleFoam`, `freeFallingRhoPimpleFoam`, `freeFallingSonicFoam`, and static-geometry tutorials (`steelBallRun-pimple`, `rodsFromGod-pimple`, `rodsFromGodRefined-pimple`).
> - **In Progress / Not Entirely Ready**: `RFGTumble` / 6DoF tumbling dynamics (`freeFalling6DoFPimpleFoam`, `RFGTumble2D`, `rodsFromGod-tumble`). Angular stability, aerodynamic damping tuning, and overset grid interpolation for high-aspect-ratio tumbling bodies are currently undergoing active testing.

---

## Technical & Physics Formulation

Simulating high-speed atmospheric fall across long trajectories using traditional domain translation requires massive computational meshes or complex dynamic domain bounds. `freeFallingFoam` formulates the CFD equations in the **body-fixed accelerating reference frame**:

### 1. Frame Pseudo-Force (d'Alembert Force)
In the non-inertial accelerating frame of a body falling with frame acceleration `a_frame(t)`, an apparent inertial body force appears on the right-hand side of the momentum equation:

```
S_U = -a_frame
```

For incompressible flow, the momentum equation becomes:

```
dU/dt + div(phi, U) - div(nu * grad(U)) = -grad(p)/rho - a_frame
```

### 2. Boundary Condition (`fallingInletVelocity`)
The custom boundary condition [`fallingInletVelocity`](file:///home/lavender/OpenFoamUbu/solvers/freeFallingFoam/src/BCs/fallingInletVelocity/fallingInletVelocityFvPatchVectorField.H) prescribes the relative incoming freestream velocity on inflow faces (`(Sf & testU) < 0` where `testU` is the frame-relative velocity direction):

```
U_inlet(t) = -(u0 + a_frame * t)
```

On outflow faces (`(Sf & testU) >= 0`), a Neumann (zero-gradient) condition is applied to prevent numerical reflections of wake vortices. Optional velocity clamping via `terminalVelocity` / `maxFrameSpeed` is also supported.

### 3. Hydrostatic-Reduced Pressure Formulation (`p_rgh`)
For 6DoF tumbling bodies ([`freeFalling6DoFPimpleFoam`](file:///home/lavender/OpenFoamUbu/freeFallingFoam/solvers/freeFalling6DoFPimpleFoam/freeFalling6DoFPimpleFoam.C)), uniform frame acceleration produces a pseudo-hydrostatic pressure gradient. Integrating total pressure `p` over an asymmetric body surface introduces spurious hydrostatic forces and moments, causing rotational blowup. 

`freeFallingFoam` uses the `p_rgh` formulation:

```
p_rgh = p - rho * (f_frame . r)
```

This isolates pure **aerodynamic forces and moments** for the 6DoF rigid body motion solver.

---

## Repository Structure

```
freeFallingFoam/
├── Allwmake                      # Master compilation script
├── Allwclean                     # Master clean script
├── src/                          # Core solver source code & models
│   ├── BCs/                      # Custom boundary conditions
│   │   └── fallingInletVelocity/ # Body-fixed accelerating frame BC
│   ├── models/                   # Shared atmospheric & wind models (libfreeFallingModels.so)
│   │   └── meanFlowModel/        # Runtime selection tables (constant, powerLaw, tabulated)
│   ├── functionObjects/          # Motion tracking & logging function objects
│   │   └── fallingMotion/        # Console & motion.dat logger
│   ├── solvers/                  # True free-falling dynamic solvers
│   │   ├── shared/               # Shared predictor-corrector kinematics & logging headers
│   │   ├── tFFPimpleFoam         # Incompressible transient solver
│   │   ├── tFFRhoPimpleFoam      # Pressure-based compressible solver
│   │   └── tFFRhoCentralFoam     # Density-based central-upwind solver (Kurganov-Tadmor)
│   └── utilities/                # Dictionary generation utility (writeFFDict)
├── tutorials/                    # Active freeFallingFoam core tutorials
│   ├── steelBallRun-tFF/         # Steel ball drop with tFFPimpleFoam & motion tracking
│   └── steelBallRun-Denser/      # Steel ball drop with dense grid resolution
├── legacy/                       # Legacy & prescribed-motion solver suite
│   └── predescribedMotion/       # Prescribed-motion solvers and cases
└── visualizations/               # Simulation animations & recordings
```

---

## Solvers Overview

| Solver | Flow Regime | Frame Motion | Key Features |
| :--- | :--- | :--- | :--- |
| **`tFFPimpleFoam`** | Incompressible | Dynamic (Coupled) | True free-falling PIMPLE solver with dynamic fluid-body force coupling |
| **`tFFRhoPimpleFoam`** | Compressible | Dynamic (Coupled) | True free-falling pressure-based compressible solver with Mach auto-regime switch |
| **`tFFRhoCentralFoam`** | Transonic/Supersonic | Dynamic (Coupled) | True free-falling density-based central-upwind solver (Kurganov-Tadmor) |
| **`freeFallingPimpleFoam`** (Legacy) | Incompressible | Prescribed | Prescribed motion PIMPLE solver |
| **`freeFallingRhoPimpleFoam`** (Legacy) | Compressible | Prescribed | Prescribed motion compressible solver |
| **`freeFallingSonicFoam`** (Legacy) | Transonic/Supersonic | Prescribed | Prescribed motion sonic solver |

---

## Building & Installation

### Prerequisites
- **OpenFOAM**: Compatible with OpenFOAM v2006+ / v2412 (Foundation & ESI/OpenCFD versions).

### Compilation
Build the boundary condition library and all solvers using [`Allwmake`](file:///home/lavender/OpenFoamUbu/freeFallingFoam/Allwmake):

```bash
cd freeFallingFoam
./Allwmake
```

### Cleaning Binaries
Clean built object files and executables using [`Allwclean`](file:///home/lavender/OpenFoamUbu/freeFallingFoam/Allwclean):

```bash
./Allwclean
```

---

## Usage & Configuration

To set up a case, add a `fallingFrameDict` inside `constant/`:

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

// Initial velocity of the body at t = 0
u0                  ( 0  0  0 );

// Frame acceleration vector (e.g., gravity in -y direction)
acceleration        ( 0 -9.81 0 );

// Body mass [kg] (required)
mass                98018;

// Wall patch(es) for aerodynamic force integration (required)
patches             ( sphere );

// Motion degrees of freedom mask: (0 1 0) for 1D drop, (1 1 1) for 3D
freedom             ( 0 1 0 );

// Reference density [kg/m^3]
rhoInf              1.204;

// Altitude and ground collision controls (optional)
initialAltitude     10000;
groundAltitude      0;
groundInterrupt     on;

// Terminal velocity clamping (optional)
terminalVelocity    85;

// Output logging to postProcessing/fallingMotion/motion.dat (default: true)
// Note: Set to false if you are using the fallingMotion function object
writeMotionLog      true;

// External body forces & propulsion (optional)
externalForces
{
    enabled         yes;

    forces
    {
        // Global-frame constant load
        constantLoad
        {
            type        constantForce;
            active      yes;
            force       ( 0 0 -50 );       // [N]
        }

        // Global-frame rocket booster
        booster
        {
            type        thrust;
            active      yes;
            direction   ( 0 1 0 );         // Global unit vector
            magnitude   1200;              // [N]
            startTime   0.0;               // [s]
            endTime     4.0;               // [s]
        }

        // Prescribed time-varying force profile
        prescribedLoad
        {
            type        tabulatedForce;
            active      no;
            interpolation  linear;
            outOfRange    zero;            // zero or hold
            table
            (
                (0.0  (0   0 0))
                (0.5  (0 200 0))
                (1.0  (0   0 0))
            );
        }

        // Body-fixed rotating thrust for 6DoF bodies
        rocketThrust
        {
            type              6DoFThrust;
            active            yes;
            direction         ( 1 0 0 );    // Longitudinal axis in body frame
            magnitude         1200;         // [N]
            startTime         0.0;          // [s]
            endTime           4.0;          // [s]
            applicationPoint  (-0.45 0 0);  // Relative to CM in body frame [m]
        }
    }
}
```

In `0/U`, configure the `fallingInletVelocity` boundary condition on far-field patches:

```cpp
inlet
{
    type            fallingInletVelocity;
    // Optional: inline parameters override constant/fallingFrameDict
    // u0           ( 0 0 0 );
    // acceleration ( 0 -9.81 0 );
    // terminalVelocity 85;
    // value        uniform (0 0 0);
}
```

---

## Tutorials & Status

### 1. Steel Ball Drop ([`tutorials/steelBallRun-tFF`](file:///home/lavender/OpenFoamUbu/solvers/freeFallingFoam/tutorials/steelBallRun-tFF))
- **Status**: Working & Verified
- **Physics**: 2 cm steel sphere falling from rest in air under standard gravity (`g = 9.81 m/s²`).
- **Features**: Axisymmetric wedge mesh, dynamic force coupling, drag calculation, live motion tracking script (`plot_motion.py`).
- **Run**:
  ```bash
  cd tutorials/steelBallRun-tFF
  ./Allrun
  ```

### 2. Legacy Prescribed-Motion Tutorials ([`legacy/predescribedMotion/tutorials`](file:///home/lavender/OpenFoamUbu/solvers/freeFallingFoam/legacy/predescribedMotion/tutorials))
- **Status**: Legacy Prescribed Motion
- **Cases**: `steelBallRun-pimple`, `rodsFromGod-pimple`, `rodsFromGodRefined-pimple`, `rodsFromGod-rhoPimple`, `rodsFromGod-rhoPimple+sonic`, `RFGTumble2D`, `rodsFromGod-tumble`.

---

## License
Built upon the OpenFOAM CFD framework (GNU General Public License v3).
