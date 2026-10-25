# OpenFOAM 6DOF Propeller Tumbling Simulation Setup Guide

This guide details the complete numerical, geometric, and dynamic setup required to simulate full 3D tumbling of the 10-inch Nylon propeller in OpenFOAM.

---

### Executive Summary: 1DOF Axial Rotation vs. 6DOF Tumbling

| Feature | Current Case (`dummySim`) | Full Tumbling Case (`tumbleSim`) |
| :--------------------------- | :------------------------------------------------------ | :------------------------------------------------------------------------- |
| **Degrees of Freedom** | 1 DOF (Rotation around Y axis only: ω_y) | Full 6 DOFs (3 translations: x, y, z; 3 rotations: roll, pitch, yaw) |
| **Rotation Axis** | Fixed vertical axis (0 1 0) | Freely moving 3D rotation vector ω(t) = (ω_x, ω_y, ω_z) |
| **Mesh Strategy** | Cylindrical sliding AMI interface (`solidBodyMotion`) | **Overset Mesh (Chimera Grid)** or Dynamic Mesh Morphing |
| **Max Rotation Angle** | Infinite continuous rotation (smooth sliding) | Infinite with Overset; limited to ~35° with standard mesh morphing |
| **Inertia Input** | Single scalar: I_yy = 3.95 × 10⁻⁵ kg·m² | Full tensor: (I_xx, I_yy, I_zz) = (3.92e-5, 3.95e-5, 3.77e-7) kg·m² |
| **Solver** | `pimpleFoam` | **`overPimpleDyMFoam`** (Recommended) or `pimpleFoam` (morphing) |

---

### 1. The Fundamental Meshing Challenge of 3D Tumbling

In `dummySim`, we used a cylindrical sliding Arbitrary Mesh Interface (AMI). Because the cylinder was rotating strictly around its own centerline ($Y$), its geometric envelope never changed shape or collided with the outer domain.

In a **tumbling** scenario, the propeller pitches end-over-end (around $X$) and rolls (around $Z$). A cylinder pitching end-over-end cannot slide inside a stationary cylinder without deforming.

There are three ways to model this in OpenFOAM:

#### Option A: Spherical Domain & Cavity with AMI (Sliding Interface) — [HIGHLY RECOMMENDED FOR MOVING FRAME]
- **Architecture:** 
  1. The rotating domain is a **solid sphere** ($D \approx 0.290\text{ m}$) enclosing the propeller.
  2. The stationary outer domain has a matching **spherical cavity** ($D \approx 0.290\text{ m}$).
  3. When combined with a moving reference frame (`freeFallingFoam` / `tFF6DoFPimpleFoam`) where free-fall translation is absorbed into frame acceleration (`-a_frame`), the sphere remains concentric at $(0, 0, 0)$.
- **Why it is ideal:** A sphere is geometrically invariant under *any* 3D rotation ($R \cdot S^2 = S^2$ for all $R \in \text{SO}(3)$). The inner sphere and outer spherical cavity stay concentric at all pitch, yaw, and roll angles, enabling continuous 360° tumbling without mesh deformation or overset hole cutting!

#### Option B: Overset Mesh (Chimera Grid) — [RECOMMENDED FOR UNCONSTRAINED DOMAINS]
- **Architecture:** 
  1. A body-fitted component mesh (ideally spherical) tightly surrounds the propeller.
  2. A large, stationary Cartesian background mesh covers the entire fall domain.
  3. The component mesh moves as a rigid body with all 6 DOFs using `overPimpleDyMFoam`.
- **Why it is robust:** Permits unconstrained 3D translation and multi-turn tumbling without mesh inversion.

#### Option C: Mesh Morphing / Smoothing (`sixDoFRigidBodyMotion`)
- **Limitation:** As soon as the propeller flips beyond **~30° to 45°**, boundary cells invert and the solver terminates with `negative cell volume`. Not recommended.

---

### 2. Propeller Mass Properties & Inertia Tensor

From CAD integration of the 3D STL mesh ([patch_2_1.stl](file:///home/lavender/OpenFoamUbu/solvers/freeFallingFoam/tutorials/propFall/meshing/rotatingMesh/constant/triSurface/patch_2_1.stl)) for solid Nylon PA6 (density = 1.14 g/cm³):

```
Total Mass (m):               0.01151 kg  (11.51 grams)
Center of Mass (CofM):        (0.0  0.0  0.0) meters

Principle Moments of Inertia:
  I_xx (Lateral / Pitch axis):     3.92 × 10⁻⁵ kg·m²  (39.2 g·cm²)
  I_yy (Axial Spin / Yaw axis):    3.95 × 10⁻⁵ kg·m²  (39.5 g·cm²)
  I_zz (Longitudinal Blade Span):  3.77 × 10⁻⁷ kg·m²  ( 0.38 g·cm²)
```

> [!IMPORTANT]
> **The Intermediate Axis / Dzhanibekov Phenomenon:**
> Notice that I_xx and I_yy are approximately **100 times larger** than I_zz. Because the mass is distributed along the blade span ($Z$), rotating about the span axis has almost zero radius of gyration.
> Under Euler's rotational equations, having two large moments of inertia and one very small moment of inertia creates intense dynamic coupling. Any small aerodynamic disturbance will cause the propeller to transition rapidly into a chaotic 3D tumbling flip.

---

### 3. Step-by-Step 6DOF Model Setup

#### Step 1: `constant/dynamicMeshDict` (6DOF Rigid Body Motion)

Below is the production-ready dictionary setup for the 6DOF rigid body solver in OpenFOAM:

```cpp
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2412                                  |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      dynamicMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dynamicFvMesh   dynamicOversetFvMesh; // Use dynamicMotionSolverFvMesh if morphing

motionSolverLibs ("libsixDoFRigidBodyMotion.so");

solver          sixDoFRigidBodyMotion;

sixDoFRigidBodyMotionCoeffs
{
    patches         (propeller);
    innerDistance   0.05;
    outerDistance   0.35;

    mass            0.01151; // kg (Nylon PA6)
    centreOfMass    (0 0 0);

    // Full principle moment of inertia vector (I_xx I_yy I_zz)
    momentOfInertia (3.919575e-05 3.945599e-05 3.769299e-07);

    // Initial orientation: unit matrix [e1x e1y e1z e2x e2y e2z e3x e3y e3z]
    // Or tilted by an initial angle to trigger tumble (see Section 4)
    orientation
    (
        1 0 0
        0 1 0
        0 0 1
    );

    // Velocity / Rotation damping factors (1.0 = undamped natural free fall)
    velocityDamping     1.0;
    angularDamping      1.0;

    // Numerical integration scheme for 6DOF ODE
    // Options: Newmark, symplectic, CrankNicolson
    accelerationRelaxation 0.9;
    accelerationDamping    0.95;

    // Center of rotation for force/moment calculations
    rho             rhoInf;
    rhoInf          1.225;

    // Active Degrees of Freedom
    // If you want full 6DOF (3 translations + 3 rotations):
    // Leave restraints empty!
    restraints
    {
        // No restraints = fully unconstrained tumbling free fall
    }
}

// ************************************************************************* //
```

---

#### Step 2: The Tumble Trigger (Initial Perturbation)

In a pure numerical simulation, if a perfectly symmetric propeller is placed normal to a uniform vertical flow with $0^\circ$ tilt and initial angular velocity = 0, the net lateral aerodynamic forces are theoretically zero. It can take tens of thousands of time steps for round-off noise to trigger a flip.

To model realistic aerodynamic tumbling, apply one of the following initial conditions:

1. **Initial Pitch Angle (Recommended):**
   Tilt the propeller by $10^\circ$ to $15^\circ$ around the $X$-axis in the `orientation` tensor:

   ```cpp
   // Tilted 15 degrees around X-axis (pitch-up):
   // cos(15°) = 0.9659, sin(15°) = 0.2588
   orientation
   (
       1.0000   0.0000   0.0000
       0.0000   0.9659  -0.2588
       0.0000   0.2588   0.9659
   );
   ```

2. **Initial Angular Rate:**
   Give the body an initial rotational impulse in `initialAngularVelocity`:

   ```cpp
   initialAngularVelocity (10.0 0.0 0.0); // 10 rad/s pitch tumble
   ```

3. **Center of Gravity (CG) Offset:**
   Slightly offset the center of mass from the geometric hub (e.g., $1\text{ mm}$ forward):

   ```cpp
   centreOfMass (0.001 0 0);
   ```

   This creates an immediate aerodynamic moment arm M = F_drag × Δx, instigating continuous tumbling.

---

#### Step 3: Governing Dynamic Equations Solved by OpenFOAM

At each CFD time step, OpenFOAM integrates the surface pressure and viscous shear stresses over the `propeller` boundary patch:

```
Aerodynamic Force:   F_aero = ∫ (-p · n + τ · n) dA
Aerodynamic Moment:  M_aero = ∫ (r - r_cm) × (-p · n + τ · n) dA
```

The 6DOF ODE solver then updates body translation and rotation using Euler's equations of motion:

1. **Translational Acceleration:**

   ```
   m · (d v_cm / dt) = F_aero + m · g
   ```

2. **Rotational Acceleration (Body-fixed frame):**

   ```
   I_xx · (dω_x/dt) + (I_zz - I_yy) · ω_y · ω_z = M_x
   I_yy · (dω_y/dt) + (I_xx - I_zz) · ω_z · ω_x = M_y
   I_zz · (dω_z/dt) + (I_yy - I_xx) · ω_x · ω_y = M_z
   ```

3. **Orientation Update:**
   OpenFOAM uses **quaternions** (q = [q0, q1, q2, q3]) to update the orientation matrix without encountering gimbal lock during full $360^\circ$ tumbling.

---

#### Step 4: Overset Mesh Directory Structure & Workflow

To build the tumbling case using `overPimpleDyMFoam`:

```
propFall/cases/tumbleOverset/
├── 0/
│   ├── U
│   ├── p
│   ├── pointDisplacement
│   └── zoneID                 <-- Generated by setFields (0=background, 1=component)
├── constant/
│   ├── dynamicMeshDict        <-- 6DOF motion specification
│   ├── polyMesh/              <-- Merged mesh (background + propeller component)
│   ├── transportProperties
│   └── turbulenceProperties
└── system/
    ├── controlDict            <-- application: overPimpleDyMFoam
    ├── fvSchemes              <-- overset schemes included
    ├── fvSolution             <-- solvers for overset and PIMPLE
    └── setFieldsDict          <-- marks cellZones for overset hole cutting
```

#### Workflow Steps

1. **Background Mesh:** Generate a simple Cartesian background box (`blockMesh`).
2. **Propeller Component Mesh:** Generate a refined cylindrical/box mesh tightly enclosing [patch_2_1.stl](file:///home/lavender/OpenFoamUbu/solvers/freeFallingFoam/tutorials/propFall/meshing/rotatingMesh/constant/triSurface/patch_2_1.stl) (`snappyHexMesh`), with the outer boundary defined as `type overset`.
3. **Merge Meshes:** Run `mergeMeshes . . -overwrite` to combine both meshes into a single case.
4. **Mark Overset Zone:** Run `setFields` to initialize `zoneID` (0 for background cells, 1 for component cells).
5. **Solve:** Run `overPimpleDyMFoam -parallel`.

---

### 5. Summary Recommendation

- **For full, realistic, continuous tumbling:** Use **`overPimpleDyMFoam`**. It is the only standard method in OpenFOAM capable of handling multi-turn tumbling without mesh inversion.
- **For immediate testing / early tumble onset (< 30°):** You can run standard `pimpleFoam` with `sixDoFRigidBodyMotion` on your current mesh by changing `dynamicMeshDict` to the dictionary provided in Section 3.
