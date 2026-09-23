
Yes — **I like this addition**, with one important correction to the wording:

> The atmosphere should be a function of **altitude**, not altitude *and falling speed*.

The instantaneous falling speed matters to the **relative flow/Mach/Reynolds number**, which then affects the aerodynamic forces. But the ambient atmospheric density itself is primarily determined by the atmospheric model at the current altitude.

So the chain should be:

$$
z(t) \rightarrow T(z),p(z),\rho(z),g(z)
$$

and independently:

$$
V(t),\rho(t),T(t) \rightarrow Re,\;Ma,\;F_\mathrm{aero}
$$

For a standard atmosphere, NASA's models explicitly provide density, pressure and temperature as functions of altitude; the U.S. Standard Atmosphere 1976 is piecewise-defined through atmospheric layers. ([NASA Technical Reports Server][1])

## One other important distinction

I'd make **three independent concepts** in the solver:

```text
densityModel
gravityModel
atmosphereModel
```

Rather than making `altitudeBased` one giant option.

Because eventually you might want:

```text
densityModel  atmosphere;
gravityModel  altitudeBased;
```

or:

```text
densityModel  rhoInf;
gravityModel  constant;
```

or eventually:

```text
atmosphereModel US1976;
```

That gives you much more flexibility.

For your intended atmospheric fall, I'd implement:

```text
densityModel  altitudeBased;
gravityModel  altitudeBased;
```

with an atmosphere dictionary providing the required constants/model.

Also, **altitude should be an actual solver state**:

$$
z_{n+1}=z_n+V_{z,n}\Delta t
$$

rather than something reconstructed from time. This becomes particularly important once the body accelerates dynamically.

And because you want a fall all the way to the ground, I'd define a clear **reference altitude / ground altitude** rather than assuming computational `y=0` automatically means sea level.

---

# Design document

I'd put this beside the true-free-fall design as something like:

`docs/altitudeAtmosphere.md`

````md
# Altitude-Dependent Atmosphere and Gravity

## Purpose

`trueFreeFallingPimpleFoam` should support both constant-property and
altitude-dependent atmospheric simulations.

The current solver can use a constant reference density:

```text
densityModel  rhoInf;
rhoInf        1.225;
````

The new formulation additionally supports:

```text
densityModel  altitudeBased;
```

When `altitudeBased` is selected, atmospheric properties are obtained
from the current physical altitude rather than from a single constant
density.

The same framework should support altitude-dependent gravitational
acceleration.

---

# 1. Density Models

The solver shall support at least two density models.

## Constant Density

```text
densityModel  rhoInf;

rhoInf        1.225;
```

The density remains constant throughout the simulation:

$$
\rho=\rho_\infty
$$

This is intended for incompressible or simplified test cases.

---

## Altitude-Based Density

```text
densityModel  altitudeBased;
```

The atmospheric density is calculated from the current altitude.

Conceptually:

$$
z=z(t)
$$

$$
\rho=\rho(z)
$$

The instantaneous falling velocity is NOT used directly to determine
atmospheric density.

Velocity instead affects the aerodynamic solution through quantities
such as Reynolds number and Mach number.

---

# 2. Atmospheric Model

Altitude-dependent density should be implemented through an atmospheric
model rather than through a hard-coded density equation.

Example:

```text
atmosphere
{
    model       USStandard1976;

    reference
    {
        altitude        0;
        temperature     288.15;
        pressure        101325;
    }

    gas
    {
        R               287.05287;
    }
}
```

The atmospheric model may provide:

```text
temperature(altitude)
pressure(altitude)
density(altitude)
speedOfSound(altitude)
```

For an ideal gas:

$$
\rho=\frac{p}{RT}
$$

The U.S. Standard Atmosphere 1976 provides piecewise atmospheric
properties as functions of geopotential altitude. The standard atmosphere
is an idealized reference atmosphere rather than a prediction of the
actual weather at a particular location.

---

# 3. Geometric vs Geopotential Altitude

The atmosphere model should distinguish between geometric altitude and
geopotential altitude where required.

For sufficiently low-altitude drone simulations, the difference is
small, but the atmosphere implementation should retain the distinction
so that higher-altitude simulations remain physically consistent.

The U.S. Standard Atmosphere 1976 uses geopotential altitude in its
definition.

The atmosphere interface should therefore internally handle the
conversion rather than requiring the user to manually convert altitude.

---

# 4. Initial Altitude

The solver must define the physical altitude of the body at the start of
the simulation.

Example:

```text
initialAltitude   500;
groundAltitude    0;
```

This means the body begins at 500 m above the selected ground/reference
level.

The altitude is then dynamically updated.

For purely vertical motion:

$$
z^{n+1}
=
z^n+V_z^n\Delta t
$$

For a general 3D body:

$$
\mathbf{x}^{n+1}
=
\mathbf{x}^n+\mathbf{V}^n\Delta t
$$

and altitude is obtained from the appropriate vertical coordinate.

---

# 5. Altitude-Dependent Gravity

The solver should support:

```text
gravityModel constant;
```

and:

```text
gravityModel altitudeBased;
```

For constant gravity:

```text
g
(
    0
    -9.80665
    0
);
```

For altitude-dependent gravity:

$$
g(r)=\frac{GM}{r^2}
$$

where

$$
r=R_\mathrm{Earth}+z
$$

and \(z\) is the geometric altitude above the selected reference
radius.

For low-altitude simulations this correction is small, but it provides a
physically consistent model for higher-altitude falls.

The Earth gravitational constant and WGS 84 reference parameters can be
used as the basis for the implementation. WGS 84 defines the Earth
semi-major axis and geocentric gravitational constant among its defining
parameters.

---

# 6. Optional Local Gravity

For higher precision near the Earth's surface, the solver may eventually
support a local gravity model depending on latitude and altitude.

Example:

```text
gravityModel local;

latitude       48.137;
longitude      11.575;
```

The initial implementation should not require this.

Recommended implementation order:

1. constant gravity
2. altitude-dependent spherical gravity
3. latitude + altitude dependent normal gravity

The final option is useful when simulations are tied to a specific
geographic location.

---

# 7. Atmospheric Reference Location

The atmosphere model should distinguish between:

```text
referenceAltitude
groundAltitude
initialAltitude
```

These are not necessarily the same quantity.

Example:

```text
referenceAltitude  0;
groundAltitude     500;
initialAltitude    1500;
```

This would represent a simulation beginning 1000 m above the ground,
while the atmospheric model is referenced to sea level.

---

# 8. Runtime Update

For every timestep:

```text
1. Determine current body altitude.
2. Evaluate atmospheric model.
3. Update rho, p, T and other atmospheric properties.
4. Evaluate gravitational acceleration.
5. Solve the CFD equations.
6. Calculate aerodynamic forces.
7. Calculate net body force.
8. Update body acceleration.
9. Update body velocity.
10. Advance body position/altitude.
11. Repeat.
```

The important distinction is:

```text
Altitude
   │
   ├──► atmospheric density
   ├──► atmospheric pressure
   ├──► atmospheric temperature
   └──► gravitational acceleration

Velocity
   │
   └──► relative flow / aerodynamic forces
```

---

# 9. Interaction With Natural Terminal Velocity

The terminal velocity is no longer a prescribed quantity.

At every altitude:

$$
m\mathbf{a}
=
m\mathbf{g}(z)
+
\mathbf{F}_{aero}
$$

The aerodynamic force depends on the local atmospheric state and
relative velocity:

$$
\mathbf{F}_{aero}
=
\mathbf{F}_{aero}
(
\rho(z),
T(z),
p(z),
\mathbf{V}
)
$$

Therefore the terminal velocity can itself change with altitude.

A local terminal velocity satisfies:

$$
\mathbf{F}_{aero}
+
m\mathbf{g}(z)
=
0
$$

Thus:

$$
V_t=V_t(z)
$$

rather than necessarily being a single constant value throughout the
fall.

This is particularly important for high-altitude simulations.

---

# 10. Compressible vs Incompressible

The atmospheric model must remain independent of whether the CFD solver
is compressible or incompressible.

For example:

```text
densityModel  rhoInf;
rhoInf        1.225;
```

may be used by an incompressible solver.

The same solver family may instead use:

```text
densityModel  altitudeBased;

atmosphere
{
    model USStandard1976;
}
```

For compressible simulations, the atmospheric model should provide
consistent initial/reference thermodynamic conditions.

For incompressible simulations, the altitude-dependent density can be
treated as an externally prescribed atmospheric property.

The implementation must avoid creating an inconsistent thermodynamic
state when used by a compressible solver.

---

# 11. Ground Termination

The solver should detect when the body reaches the specified ground
altitude.

Conceptually:

```text
if (altitude <= groundAltitude)
{
    terminateSimulation();
}
```

This should be implemented independently of the atmospheric model.

Ground altitude is a physical simulation parameter, not necessarily the
same as the atmosphere reference altitude.

---

# 12. Suggested Dictionary Structure

A complete configuration could eventually look like:

```text
atmosphere
{
    densityModel       altitudeBased;
    gravityModel       altitudeBased;

    initialAltitude    1000;
    groundAltitude     0;

    atmosphereModel
    {
        type            USStandard1976;

        referenceAltitude 0;
    }

    gravity
    {
        model           sphericalEarth;

        earthRadius     6.371e6;
        gravitationalParameter 3.986004418e14;
    }
}
```

The simple constant-density case remains:

```text
atmosphere
{
    densityModel       rhoInf;
    rhoInf             1.225;

    gravityModel       constant;
    g                  (0 -9.80665 0);

    initialAltitude    0;
    groundAltitude     -1e30;
}
```

---

# 13. Design Principle

Atmospheric properties should be determined by the physical location of
the body.

The solver should therefore follow:

```text
              Position / Altitude
                      │
          ┌───────────┴───────────┐
          ▼                       ▼
   Atmosphere model          Gravity model
          │                       │
     rho, p, T, a                 g
          │                       │
          └───────────┬───────────┘
                      ▼
                 CFD solution
                      │
                      ▼
               Aerodynamic force
                      │
                      ▼
              Body acceleration
                      │
                      ▼
               Body velocity
                      │
                      ▼
                 Position
                      │
                      └──────────────►
```

This allows atmospheric properties and gravitational acceleration to
change naturally during a falling-body simulation without introducing
them into the terminal-velocity calculation as separate assumptions.

````

## One thing I'd add to the eventual implementation

I'd make the atmosphere system **more general than `altitudeBased`**.

Something like:

```text
atmosphereModel
{
    type USStandard1976;
}
````

Then later you can have:

```text
type USStandard1976;
type ISA;
type tabulated;
type exponential;
type custom;
```

That way `densityModel altitudeBased` means *"density comes from an atmospheric model evaluated at altitude"* rather than *"the solver contains one particular density-vs-altitude formula."*

That's a much cleaner architecture.

And for your drone scenario, **500–1000 m → ground**, this is probably a relatively small correction compared with the aerodynamic uncertainties, but it's still an excellent feature to build into `trueFreeFallingPimpleFoam` now. The U.S. Standard Atmosphere has sea-level \(\rho=1.225\,\mathrm{kg/m^3}\), \(T=288.15\,K\), and \(g_0=9.80665\,m/s^2\), so your existing `rhoInf 1.225` case can remain a clean baseline against which the altitude-dependent version is validated. ([NASA Technical Reports Server][1])

One terminology change I'd strongly recommend: **call the input `initialAltitude`, not `startingAltitude`**, and call the dynamic variable simply `altitude`. It will make the code/readme much clearer once 6DoF enters the picture.

[1]: https://ntrs.nasa.gov/api/citations/20060053240/downloads/20060053240.pdf?utm_source=chatgpt.com
