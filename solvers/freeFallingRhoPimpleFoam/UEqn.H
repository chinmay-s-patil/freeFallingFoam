// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

// Compute current frame velocity for reference / BC update
const scalar t = runTime.value();
const vector frameVelocity = u0 + frameAcc * t;

// Clamp if maxFrameSpeed is set
const scalar frameSpeed = mag(frameVelocity);
const vector frameVelocityClamped =
    (frameSpeed > maxFrameSpeed && frameSpeed > SMALL)
    ? frameVelocity * (maxFrameSpeed / frameSpeed)
    : frameVelocity;

if (runTime.writeTime())
{
    Info<< "  Frame velocity U_frame = " << frameVelocityClamped
        << "  |U_frame| = " << mag(frameVelocityClamped) << " m/s" << nl;
}

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(rho, U) + fvm::div(phi, U)
  + MRF.DDt(rho, U)
  + turbulence->divDevRhoReff(U)
 ==
    fvOptions(rho, U)
  + rho*fallingFrameForce
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
    K = 0.5*magSqr(U);
}

// ============================================================
//  Update magUStar = mag(U) / max(mag(acc*t), SMALL)
//
//  Normalises by the *current* free-fall speed (a*t), so values
//  approach 1 once the fluid has accelerated to match the frame.
//  Uses the clamped speed so the field saturates gracefully past
//  terminal velocity.
// ============================================================
{
    const scalar Uref = max(mag(frameVelocityClamped), SMALL);
    magUStar = mag(U) / dimensionedScalar("Uref", dimVelocity, Uref);
}
