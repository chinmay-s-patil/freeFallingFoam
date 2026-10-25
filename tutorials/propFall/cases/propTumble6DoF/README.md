# propTumble6DoF: Full 6DoF Propeller Free-Fall & Tumbling Simulation

## Overview
`propTumble6DoF` models the unconstrained 6DoF motion of a 10-inch propeller (Nylon PA6, 11.51 g) in free fall under gravity:
- **3 Translational DOFs:** Freely falling along $Y$ plus full lateral drift along $X$ and $Z$ driven by aerodynamic side forces ($\vec{a}_{\text{frame}} = (\vec{F}_{\text{net}})/m$).
- **3 Rotational DOFs:** Unconstrained 3D attitude motion (pitch, spin, roll) integrated via Euler's equations with RK4 and quaternion attitude propagation across a spherical AMI interface.

---

## Key Differences from `propTumble4DoF`

| Feature | `propTumble4DoF` | `propTumble6DoF` |
| :--- | :--- | :--- |
| **Kinematic Freedom** | `freedom (0 1 0);` (vertical drop only) | `freedom (1 1 1);` (full 3D translation) |
| **Active DOFs** | 4 (1 Translation: $Y$ + 3 Rotations: $\phi, \theta, \psi$) | 6 (3 Translations: $X, Y, Z$ + 3 Rotations: $\phi, \theta, \psi$) |
| **Inertia Tensor Formulation** | Unrotated diagonal principal inertia | Full 3D inertia tensor $\mathbf{J}_b = \mathbf{R}_0 \mathbf{I}_{\text{principal}} \mathbf{R}_0^T$ accounting for the 15° mesh pre-tilt |
| **Kinematics Tracked** | $V_y, A_y, \omega_x, \omega_y, \omega_z$ ($V_x=V_z=0$) | $V_x, V_y, V_z, A_x, A_y, A_z, \omega_x, \omega_y, \omega_z$ |

---

## Physical Properties & Dynamics Configuration

In `constant/fallingFrameDict`:
```foam
// Enable full 3D translation
freedom             ( 1 1 1 );

rotor
{
    hubOrigin         (0 0 0);

    // Principal moments of inertia [kg*m^2]
    inertiaPrincipal  (3.919575e-05 3.945599e-05 3.769299e-07);

    // Mesh is physically pre-tilted 15 deg about +X in polyMesh
    meshTiltAngleDeg  15.0;
    meshTiltAxis      (1 0 0);

    // Full 3D rotational freedom
    rotationFreedom   (1 1 1);
    initialOmega      (0 0 0);
    initialAngleDeg   0.0;
    damping           (0 0 0);
    report            true;
}
```

---

## Running the Simulation

### Parallel Execution (Default: 8 Cores)
```bash
./Allrun
```
This automatically decomposes the domain using Scotch, runs `tFF6DoFPimpleFoam -parallel`, and reconstructs time steps.

### Serial Execution
```bash
tFF6DoFPimpleFoam
```

### Clean Case
```bash
./Allclean
```

### Monitoring & Visualization
- **Live Terminal Dashboard:**
  ```bash
  python3 track_motion.py --watch
  ```
- **Generate 4-Panel Analysis Plots:**
  ```bash
  python3 track_motion.py --plot
  ```
- **Export Combined 6DoF Kinematics to CSV:**
  ```bash
  python3 track_motion.py --export tumble_6dof.csv
  ```
- **3D Animation Rendering:**
  ```bash
  python3 animate_motion.py
  ```
