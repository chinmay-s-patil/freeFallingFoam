// ============================================================
//  predictFallingFrameMotion.H — Predictor-Corrector Predictor Step
//
//  Semi-implicit predictor-corrector scheme:
//    Predict:  V_pred   = V^n + a^n * dt  (using old acceleration)
//    Solve:    PIMPLE / density equations with V_pred BCs
//    Correct:  a^{n+1}  from converged aero forces
//              V^{n+1}  = V^n + 0.5*(a^n + a^{n+1})*dt (trapezoidal)
//
//  This is O(dt²) accurate. It is NOT a full explicit Heun scheme
//  because the predictor does not recompute aerodynamic forces — it
//  reuses a^n from the previous timestep's corrector.
//
//  NOTE: Variables declared here (dt, gDir, V_up_prev, etc.) remain
//  in scope for correctFallingFrameMotion.H because both headers
//  are #include'd at the same brace level inside the time loop body.
// ============================================================

// Store state at beginning of timestep (time n)
const scalar dt = runTime.deltaTValue();
a_prev = frameAcceleration.value();
v_prev = frameVelocity.value();
alt_prev = altitude.value();

// Determine vertical direction from gravity vector (-g / |g|)
vector gDir = gravity.value();
const scalar magG = mag(gDir);
if (magG > SMALL)
{
    gDir /= magG;
}
else
{
    gDir = vector(0, -1, 0);
}

// ------------------------------------------------------------
//  1. Predict velocity: V_pred = V^n + a^n * dt
// ------------------------------------------------------------
const vector V_pred = v_prev + a_prev * dt;

// ------------------------------------------------------------
//  2. Predict altitude: alt_pred = alt^n + V_up^n * dt
// ------------------------------------------------------------
const scalar V_up_prev = -(v_prev & gDir);
const scalar alt_pred = alt_prev + V_up_prev * dt;

frameVelocity.value() = V_pred;
altitude.value() = alt_pred;

// ------------------------------------------------------------
//  3. Update atmospheric density and gravity at predicted altitude
// ------------------------------------------------------------
if (densityModel == "altitudeBased")
{
    const scalar effectiveAlt = max(alt_pred - groundAltitude, 0.0);
    rhoAtmosphere.value() = rhoRef * std::exp(-effectiveAlt / scaleHeight);
}
else
{
    rhoAtmosphere.value() = rhoRef;
}

if (gravityModel == "altitudeBased")
{
    const scalar magG0 = max(mag(initialGravity), SMALL);
    const scalar effectiveAlt = max(alt_pred - groundAltitude, 0.0);
    const scalar gMag = magG0 * sqr(earthRadius / (earthRadius + effectiveAlt));
    gravity.value() = (initialGravity / magG0) * gMag;
}
else
{
    gravity.value() = initialGravity;
}

// ------------------------------------------------------------
//  4. Update mean flow and relative freestream velocity
// ------------------------------------------------------------
vector Umean_pred = vector::zero;
if (meanFlowPtr.valid())
{
    Umean_pred = meanFlowPtr->velocity(altitude.value());
}
meanFlowVelocity.value() = Umean_pred;
relativeVelocity.value() = meanFlowVelocity.value() - frameVelocity.value();

// ------------------------------------------------------------
//  5. Update apparent body force in momentum equations: g - a_frame
// ------------------------------------------------------------
fallingFrameForce.value() = gravity.value() - a_prev;

#if defined(FF_INCOMPRESSIBLE)
{
    const volScalarField gh_old("gh_old", gh);
    gh = (fallingFrameForce & mesh.C());
    ghf = (fallingFrameForce & mesh.Cf());
    p_rgh += gh_old - gh;
}
#endif

// ------------------------------------------------------------
//  6. Update boundary conditions with predicted frame velocity
// ------------------------------------------------------------
U.correctBoundaryConditions();
