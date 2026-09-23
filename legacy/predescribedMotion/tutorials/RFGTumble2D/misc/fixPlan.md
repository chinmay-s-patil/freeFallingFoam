# Fix Runaway Angular Velocity in freeFalling6DoFPimpleFoam

## Problem

The velocity field in the RFGTumble2D simulation goes haywire — the angular velocity ω_z grows without bound from 0 to 11.5+ rad/s (a 6m rod spinning at ~1.8 rev/s). The rod should reach a physical equilibrium tumble rate, but instead it accelerates continuously. The time step shrinks to ~10μs as the solver struggles to keep up.

## Root Cause Analysis

I've identified **two bugs** in the solver and **one critical missing feature**:

### Bug 1: `correctPhi.H` uses `p` instead of `p_rgh`

In [correctPhi.H](file:///home/lavender/OpenFoamUbu/freeFallingFoam/solvers/freeFalling6DoFPimpleFoam/correctPhi.H), the flux correction after mesh motion uses the *total pressure* `p`:

```cpp
CorrectPhi(U, phi, p, ...);
```

But this solver uses `p_rgh` as the primary pressure variable. The `CorrectPhi` function uses the pressure to drive a flux correction solve — using `p` (which includes the hydrostatic head `gh`) means the flux correction is contaminated by the body force potential. This creates a small but persistent error in the face fluxes that biases the force/torque integration on the rod.

**Fix:** Change to `p_rgh`, consistent with the pressure equation.

### Bug 2: Body force flux inconsistency with PIMPLE consistent mode

In [pEqn.H](file:///home/lavender/OpenFoamUbu/freeFallingFoam/solvers/freeFalling6DoFPimpleFoam/pEqn.H) line 36, the body force flux is:

```cpp
phiHbyA += fvc::interpolate(rAU)*(fallingFrameForce & mesh.Sf());
```

But when `pimple.consistent()` is true, the solver switches to `rAtU` (line 42-46). The body force flux contribution should also use `rAtU` for consistency, but it was added before the consistent-mode correction block. While PIMPLE consistent mode may not be active currently, this creates a silent inconsistency.

More critically, the face-interpolated body force flux `(fallingFrameForce & mesh.Sf())` should use `ghf` (`fallingFrameForce & mesh.Cf()` dotted with the face area normal) as the standard buoyant solvers do. Our implementation directly dots the body force with face area vectors, which is algebraically equivalent for a uniform force, but should match the standard pattern for robustness with non-uniform meshes.

### Missing Feature: No angular velocity damping restraint

The 6DoF configuration in [dynamicMeshDict](file:///home/lavender/OpenFoamUbu/freeFallingFoam/tutorials/RFGTumble2D/combiRun/constant/dynamicMeshDict) has constraints (fixed point, fixed plane, fixed axis) but the `restraints` block is empty. There is **no physical mechanism** to damp the angular velocity growth.

In reality, a falling rod would have aerodynamic drag torque proportional to ω² that limits the tumble rate. But in this simulation:

1. The `accelerationRelaxation 0.4` under-relaxes the 6DoF response, which slows but doesn't prevent divergence
2. The p_rgh formulation correctly removes hydrostatic forces from the 6DoF computation, but the aerodynamic torque itself grows without natural bound because the mesh motion feedback loop amplifies small asymmetries

**Fix:** Add a `linearDamper` restraint to the 6DoF system. This provides an `M_damping = -coeff * ω` restoring torque that represents the missing aerodynamic rotational drag.

> [!IMPORTANT]
> The damping coefficient needs tuning. A reasonable first estimate based on the rod geometry: for a 6m rod at 85 m/s in air, the aerodynamic torque at equilibrium should balance the driving torque. I'll compute a physically-motivated value.

## Proposed Changes

### Solver Bug Fixes

#### [MODIFY] [correctPhi.H](file:///home/lavender/OpenFoamUbu/freeFallingFoam/solvers/freeFalling6DoFPimpleFoam/correctPhi.H)

Change `p` → `p_rgh` in the `CorrectPhi` call to be consistent with the p_rgh pressure formulation.

#### [MODIFY] [pEqn.H](file:///home/lavender/OpenFoamUbu/freeFallingFoam/solvers/freeFalling6DoFPimpleFoam/pEqn.H)

Move the body force flux addition after the consistent-mode block so it uses `rAtU` instead of `rAU` when PIMPLE consistent mode is active.

---

### Case Configuration Fix

#### [MODIFY] [dynamicMeshDict](file:///home/lavender/OpenFoamUbu/freeFallingFoam/tutorials/RFGTumble2D/combiRun/constant/dynamicMeshDict)

Add a `linearDamper` restraint in the `restraints` block:

```
restraints
{
    angularDamper
    {
        sixDoFRigidBodyMotionRestraint linearDamper;
        coeff   5000;    // N·m·s/rad — needs tuning
    }
}
```

## Open Questions

> [!IMPORTANT]
> **Damping coefficient tuning**: The `linearDamper` coefficient determines the equilibrium tumble rate. Too high → the rod barely tumbles. Too low → still diverges. Would you like me to:
>
> 1. Start with a physically-estimated value (~5000 N·m·s/rad based on the rod inertia and expected equilibrium ω of ~2 rad/s)?
> 2. Or do you have a target tumble rate in mind?

> [!NOTE]
> **Alternative approach**: Instead of `linearDamper`, we could use a `sphericalAngularDamper` which applies damping equally in all rotational directions, or implement a custom quadratic-drag restraint (`M = -C_d * ω|ω|`) which is more physically correct for aerodynamic drag. The `linearDamper` is simplest to start with but may over-damp at low ω and under-damp at high ω compared to quadratic drag.

> [!WARNING]
> **LES on a 2D case**: The solver log shows warnings about LES not being strictly applicable for 2D cases. Since the WALE model is active, it may be adding significant unphysical dissipation (or too little). Consider switching to `laminar` or a RANS model. This is tangential to the velocity blowup but affects the overall physics fidelity.

## Verification Plan

### Automated Tests

1. Rebuild the solver: `wmake` in the solver directory
2. Re-run a short test (5-10 seconds of simulation time) and check:
   - Angular velocity should stabilize or grow much more slowly
   - Time step should remain manageable (≥1e-4 s)
   - Continuity errors should remain small (~1e-6)

### Manual Verification

- Repost-process with ParaView to check velocity field smoothness
- Plot ω_z vs time using the existing `plot_angle.py` script
