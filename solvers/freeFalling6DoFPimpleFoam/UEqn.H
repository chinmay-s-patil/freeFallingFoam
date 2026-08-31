// ============================================================
//  UEqn.H  (p_rgh formulation for falling reference frame)
//
//  The body force is NO LONGER a source term on the RHS.
//  Instead it appears in the pressure gradient:
//    grad(p) = grad(p_rgh) + fallingFrameForce
//
//  This means the 6DoF solver only sees p_rgh (aerodynamic
//  pressure) when computing forces on the body.
// ============================================================

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

// Build UEqn WITHOUT the body force (it's in the pressure gradient now)
tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
 ==
    fvOptions(U)
    // NOTE: fallingFrameForce is NOT here.
    // It enters via: solve(UEqn == -grad(p_rgh) - fallingFrameForce)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
     ==
        -fvc::grad(p_rgh) - fallingFrameForce
    );

    fvOptions.correct(U);
}

{
    const scalar Uref = max(mag(frameVelocityClamped), SMALL);
    magUStar = mag(U) / dimensionedScalar("Uref", dimVelocity, Uref);
}
