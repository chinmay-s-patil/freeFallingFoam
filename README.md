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
The custom boundary condition [`fallingInletVelocity`](file:///home/lavender/OpenFoamUbu/freeFallingFoam/BCs/fallingInletVelocity/fallingInletVelocityFvPatchVectorField.H) prescribes the relative incoming freestream velocity on inflow faces (`phi < 0`):

```
U_inlet(t) = -(u0 + a_frame * t)
```

On outflow faces (`phi >= 0`), a Neumann (zero-gradient) condition is applied to prevent numerical reflections of wake vortices. Optional velocity clamping via `terminalVelocity` / `maxFrameSpeed` is also supported.

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
├── BCs/                          # Custom boundary conditions
│   └── fallingInletVelocity/    # Body-fixed accelerating frame BC
├── solvers/                      # Solver suite
│   ├── freeFallingPimpleFoam     # Incompressible transient solver
│   ├── freeFallingRhoPimpleFoam  # Compressible transient solver
│   ├── freeFallingSonicFoam      # Transonic/supersonic compressible solver
│   ├── freeFalling6DoFPimpleFoam # 6DoF tumbling solver (p_rgh formulation)
│   ├── overDyMFreeFallingPimpleFoam    # Incompressible overset dynamic mesh solver
│   └── overDyMFreeFallingRhoPimpleFoam # Compressible overset dynamic mesh solver
├── tutorials/                    # Case tutorials
│   ├── steelBallRun-pimple/      # 2D axisymmetric sphere drop (Working)
│   ├── rodsFromGod-pimple/       # High-speed rod free fall (Working)
│   ├── rodsFromGodRefined-pimple/# Mesh-refined rod simulation (Working)
│   ├── rodsFromGod-rhoPimple/    # Compressible rod drop (Working)
│   ├── rodsFromGod-rhoPimple+sonic/ # Supersonic transition test (Working)
│   ├── RFGTumble2D/              # 2D 6DoF tumbling rod (WIP - Under development)
│   └── rodsFromGod-tumble/       # 3D 6DoF overset tumbling rod (WIP)
└── visualizations/               # Pre-rendered simulation videos
```

---

## Solvers Overview

| Solver | Flow Regime | Mesh / Motion | Key Features |
| :--- | :--- | :--- | :--- |
| **`freeFallingPimpleFoam`** | Incompressible | Fixed mesh | PIMPLE algorithm, `-a_frame` body force |
| **`freeFallingRhoPimpleFoam`** | Compressible | Fixed mesh | Thermophysical energy eq., ideal gas |
| **`freeFallingSonicFoam`** | Transonic/Supersonic | Fixed mesh | Density/pressure-based sonic solver |
| **`freeFalling6DoFPimpleFoam`** | Incompressible | 6DoF Motion | `p_rgh` formulation, aerodynamic 6DoF coupling |
| **`overDyMFreeFallingPimpleFoam`** | Incompressible | Overset Mesh | Overset grid motion with frame acceleration |
| **`overDyMFreeFallingRhoPimpleFoam`** | Compressible | Overset Mesh | Compressible overset dynamic grid solver |

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

// Terminal velocity clamping (optional)
terminalVelocity    85;
```

In `0/U`, configure the `fallingInletVelocity` boundary condition on far-field patches:

```cpp
inlet
{
    type            fallingInletVelocity;
}
```

---

## Tutorials & Status

### 1. Steel Ball Drop ([`steelBallRun-pimple`](file:///home/lavender/OpenFoamUbu/freeFallingFoam/tutorials/steelBallRun-pimple))
- **Status**: Working & Verified
- **Physics**: 2 cm steel sphere falling from rest in air under standard gravity (`g = 9.81 m/s²`).
- **Features**: Axisymmetric wedge mesh, boundary layer development, drag calculation, and normalized velocity scale `magUStar = mag(U) / mag(U_frame)`.
- **Run**:
  ```bash
  cd tutorials/steelBallRun-pimple
  ./Allrun
  ```

### 2. Rods From God ([`rodsFromGod-pimple`](file:///home/lavender/OpenFoamUbu/freeFallingFoam/tutorials/rodsFromGod-pimple) & [`rodsFromGodRefined-pimple`](file:///home/lavender/OpenFoamUbu/freeFallingFoam/tutorials/rodsFromGodRefined-pimple))
- **Status**: Working & Verified
- **Physics**: High-aspect-ratio kinetic energy penetrator falling at high velocity.
- **Features**: Fine wake resolution, shock/drag dynamics, ParaView visualization scripts.

### 3. Compressible Rod Drop ([`rodsFromGod-rhoPimple`](file:///home/lavender/OpenFoamUbu/freeFallingFoam/tutorials/rodsFromGod-rhoPimple) & [`rodsFromGod-rhoPimple+sonic`](file:///home/lavender/OpenFoamUbu/freeFallingFoam/tutorials/rodsFromGod-rhoPimple+sonic))
- **Status**: Working & Verified
- **Physics**: Transonic/compressible atmosphere free fall simulation.

### 4. RFGTumble & 6DoF Dynamics ([`RFGTumble2D`](file:///home/lavender/OpenFoamUbu/freeFallingFoam/tutorials/RFGTumble2D) & [`rodsFromGod-tumble`](file:///home/lavender/OpenFoamUbu/freeFallingFoam/tutorials/rodsFromGod-tumble))
- **Status**: Work in Progress (Not Entirely Ready)
- **Details**: Tests 6DoF angular pitch oscillation, overset grid interpolation, and aerodynamic damping during free fall. Current development focuses on stabilizing angular velocity, tuning linear/angular dampers, and optimizing overset cellZone coupling.

---

## License
Built upon the OpenFOAM CFD framework (GNU General Public License v3).
