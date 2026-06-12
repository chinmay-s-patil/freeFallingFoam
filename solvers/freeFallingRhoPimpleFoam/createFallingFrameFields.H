// ============================================================
//  createFallingFrameFields.H
//
//  Reads the falling-frame parameters from transportProperties
//  (or a dedicated fallingFrameDict) and exposes:
//
//    frameVelocity(t)    = u0 + a*t          [m/s]  (body velocity)
//    frameAcceleration   = a                 [m/s2] (constant for now)
//
//  The momentum equation gains an inertial pseudo-force:
//    S_U = -rho * a(t)   (per unit volume, rho=1 for incompressible)
//
//  For a time-varying acceleration you can replace `frameAcc` with
//  a runtime-evaluated Function1.
// ============================================================

Info<< "Reading fallingFrameDict\n" << endl;

IOdictionary fallingFrameDict
(
    IOobject
    (
        "fallingFrameDict",          // file: constant/fallingFrameDict
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Initial velocity of the body at t=0  (e.g. (0 0 0))
const vector u0
(
    fallingFrameDict.get<vector>("u0")
);

// Constant acceleration of the body  (e.g. (0 -9.81 0) for free-fall)
const vector frameAcc
(
    fallingFrameDict.get<vector>("acceleration")
);

// Optional: clamp maximum frame speed (set to GREAT to disable)
const scalar maxFrameSpeed
(
    fallingFrameDict.getOrDefault<scalar>("terminalVelocity", GREAT)
);

Info<< "  Falling frame: u0 = " << u0
    << "  acceleration = " << frameAcc
    << "  maxFrameSpeed = " << maxFrameSpeed
    << nl << endl;

// Pseudo-force field: -acceleration (uniform, changes sign because
// we are in the body frame, so fluid sees +bodyAcc as apparent gravity)
// Stored as a dimensioned vector for use in UEqn source term.
//
//   DU/Dt|_frame = DU/Dt|_inertial - a_frame
//   => extra RHS term: -a_frame   (per unit mass)
//
dimensionedVector fallingFrameForce
(
    "fallingFrameForce",
    dimAcceleration,
    -frameAcc          // sign: source on RHS of momentum eq.
);

// ============================================================
//  Dimensionless velocity magnitude field
//
//   magUStar = mag(U) / max(mag(U_frame(t)), SMALL)
//
//  U_frame(t) = u0 + frameAcc * t  (clamped to maxFrameSpeed)
//
//  Provides a stable colour map while the body accelerates
// ============================================================

volScalarField magUStar
(
    IOobject
    (
        "magUStar",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("magUStar", dimless, 0.0)
);
