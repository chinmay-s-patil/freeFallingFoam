// ============================================================
//  UEqn.H  (modified for falling/accelerating reference frame)
//
//  Extra term: fallingFrameForce = -a_frame
//  This is the d'Alembert pseudo-force that appears when the
//  coordinate system accelerates with the body.
//
//  Inlet BC must also be updated each time-step (see pimpleFoam.C)
//  to prescribe U_inlet(t) = -(u0 + a*t) in the lab-frame sense,
//  which looks like an upward flow in the body-fixed frame.
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

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
 ==
    fvOptions(U)
  + fallingFrameForce          // <-- inertial pseudo-force [m/s2]
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
}

{
    const scalar Uref = max(mag(frameVelocityClamped), SMALL);
    magUStar = mag(U) / dimensionedScalar("Uref", dimVelocity, Uref);
}
