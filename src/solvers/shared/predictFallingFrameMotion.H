// ============================================================
//  predictFallingFrameMotion.H — 2nd-order Heun Predictor Step
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
    const scalar magG0 = max(mag(initialAcc), SMALL);
    const scalar effectiveAlt = max(alt_pred - groundAltitude, 0.0);
    const scalar gMag = magG0 * sqr(earthRadius / (earthRadius + effectiveAlt));
    gravity.value() = (initialAcc / magG0) * gMag;
}
else
{
    gravity.value() = initialAcc;
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
//  5. Update pseudo-force in momentum equations
// ------------------------------------------------------------
fallingFrameForce.value() = -a_prev;

#if defined(FF_INCOMPRESSIBLE)
gh = (fallingFrameForce & mesh.C());
ghf = (fallingFrameForce & mesh.Cf());
#endif

// ------------------------------------------------------------
//  6. Update boundary conditions with predicted frame velocity
// ------------------------------------------------------------
U.correctBoundaryConditions();
