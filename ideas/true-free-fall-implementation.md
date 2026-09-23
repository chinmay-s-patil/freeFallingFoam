
Yep. And **“explicit formulation” vs “implicit formulation” is not quite the right naming** here.

Your old solver is better described as:

* `prescribedAcceleration` / `prescribedMotion` — acceleration is supplied by the user.
* `forceCoupled` / `dynamicallyCoupled` — acceleration is obtained from the forces.

“Explicit” and “implicit” have a specific numerical meaning, and you don't want someone reading the repo to assume the distinction refers to that.

I'd call the old code:

```text
legacyPrescribedMotion/
```

or, if you want it very obvious:

```text
legacyFixedAcceleration/
```

I prefer **`legacyPrescribedMotion`**, because the old solver isn't necessarily always fixed acceleration; it prescribes the frame motion.

Here's the Markdown I'd put in the repo:

````md
# Natural Terminal Velocity Formulation

## Purpose

The original `freeFallingPimpleFoam` formulation prescribes the
acceleration of the falling reference frame:

\[
\mathbf{V}_{frame}(t)
=
\mathbf{V}_0+\mathbf{a}_{frame}t
\]

A `terminalVelocity` parameter can limit the imposed frame velocity,
but this does **not** represent a physical terminal velocity. The
acceleration remains prescribed independently of the aerodynamic forces.

The new `trueFreeFallingPimpleFoam` formulation removes this limitation.

The body's acceleration and velocity should emerge naturally from the
balance between gravity and the aerodynamic forces calculated by the CFD
solution.

---

## Governing Concept

For a body of mass \(m\):

\[
m\mathbf{a}
=
m\mathbf{g}
+
\mathbf{F}_{aero}
\]

Therefore,

\[
\boxed{
\mathbf{a}
=
\mathbf{g}
+
\frac{\mathbf{F}_{aero}}{m}
}
\]

The velocity is then integrated in time:

\[
\boxed{
\mathbf{V}^{n+1}
=
\mathbf{V}^{n}
+
\mathbf{a}^{n}\Delta t
}
\]

The resulting velocity is used to update the falling reference frame.

No terminal velocity is prescribed.

---

## Natural Terminal Velocity

Terminal velocity is not an input.

It emerges when the aerodynamic force balances gravity:

\[
\mathbf{F}_{aero}
+
m\mathbf{g}
\approx
0
\]

which gives

\[
\mathbf{a}\approx0
\]

while the velocity continues to remain approximately constant.

Thus:

\[
\boxed{
\text{Terminal velocity}
=
\text{steady-state velocity resulting from force balance}
}
\]

The solver must therefore never directly impose a terminal velocity.

---

## Required Changes

### 1. Remove prescribed frame acceleration

The current formulation reads a constant acceleration from the
falling-frame configuration.

This must no longer be the source of the body's acceleration.

Instead, acceleration becomes a runtime state variable:

```text
frameAcceleration
frameVelocity
````

initialized from the user's initial conditions.

---

### 2. Calculate aerodynamic force

At every timestep, obtain the aerodynamic force acting on the falling
body.

For the initial implementation, only the translational force in the
falling direction is required.

The force should be obtained from the actual CFD pressure and viscous
stresses rather than from an assumed drag coefficient.

---

### 3. Calculate net acceleration

After obtaining the aerodynamic force:

$$
\mathbf{F}_{net}
=
m\mathbf{g}
+
\mathbf{F}_{aero}
$$

then:

$$
\mathbf{a}
=
\frac{\mathbf{F}_{net}}{m}
$$

This acceleration becomes the acceleration of the falling reference frame.

---

### 4. Update velocity

Integrate the acceleration:

$$
\mathbf{V}^{n+1}
=
\mathbf{V}^{n}
+
\mathbf{a}^{n}\Delta t
$$

The resulting velocity becomes the frame velocity used by the next
fluid timestep.

---

### 5. Update the inlet velocity

The inlet boundary condition must no longer reconstruct the frame velocity
from:

```text
initialVelocity + acceleration * time
```

Instead it should obtain the current runtime frame velocity.

Conceptually:

```text
U_inlet = -frameVelocity;
```

---

### 6. Update the accelerating-frame force

The momentum equation must use the dynamically calculated frame
acceleration.

Conceptually:

```text
frameAcceleration = acceleration;
```

and the corresponding pseudo-force must be updated consistently in the
momentum equation.

---

### 7. Handle accelerating-frame pressure correctly

Changing frame acceleration changes the effective hydrostatic pressure
field.

The pressure formulation must therefore remain consistent with the
accelerating reference frame.

The existing `p_rgh` treatment used by the 6DoF solver should be used as
the basis for this implementation.

The aerodynamic force calculation must not include the artificial
pressure contribution associated with the accelerating frame.

---

### 8. Remove `terminalVelocity` as a physical input

The old parameter should either be removed or retained only as a
deprecated compatibility option.

It must **not** clamp the physical velocity.

The solver should instead determine the terminal velocity from:

```text
gravity
+
pressure forces
+
viscous forces
```

---

## Initial Coupling Strategy

The first implementation should use explicit fluid-body coupling.

At timestep `n`:

```text
1. Known frame velocity V[n]
2. Known frame acceleration a[n]
3. Solve CFD
4. Integrate aerodynamic force
5. Calculate net force
6. Calculate acceleration a[n+1]
7. Update velocity V[n+1]
8. Advance to timestep n+1
```

Do not introduce strong coupling or sub-iterations initially.

Once the basic formulation is validated, relaxation/sub-iterations can be
introduced if required for stability.

## Design Principle

The fundamental distinction is:

### Old formulation

```text
User
 ↓
acceleration
 ↓
velocity
 ↓
CFD
```

### New formulation

```text
             ┌──────────────┐
             │     CFD      │
             └──────┬───────┘
                    │
             aerodynamic force
                    │
                    ▼
             ┌──────────────┐
             │ Force balance│
             └──────┬───────┘
                    │
               acceleration
                    │
                    ▼
             ┌──────────────┐
             │   velocity   │
             └──────┬───────┘
                    │
                    ▼
             falling frame
                    │
                    └──────► CFD
```

The body's motion is therefore an output of the coupled
fluid-body system rather than an input imposed independently of the
aerodynamic forces.
