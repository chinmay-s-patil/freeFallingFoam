// ============================================================
//  createFallingFrameFields.H
//
//  Reads the falling-frame parameters from fallingFrameDict
//  and creates the hydrostatic-reduced pressure helpers for
//  the p_rgh formulation.
//
//  Fields created:
//    fallingFrameForce  — inertial pseudo-force = -acceleration
//    gh                 — hydrostatic head at cell centres
//    ghf                — hydrostatic head at face centres
//    magUStar           — normalised velocity magnitude
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

// Pseudo-force = -acceleration (sign: source on RHS of momentum eq.)
//   DU/Dt|_frame = DU/Dt|_inertial - a_frame
//   => extra RHS term: -a_frame (per unit mass)
dimensionedVector fallingFrameForce
(
    "fallingFrameForce",
    dimAcceleration,
    -frameAcc          // sign: upward in body frame for downward falling
);

// ============================================================
//  Hydrostatic head fields for p_rgh formulation
//
//  p = p_rgh + gh    where  gh = fallingFrameForce & position
//
//  This splits the total pressure into:
//    p_rgh — aerodynamic pressure only (what the 6DoF solver sees)
//    gh    — hydrostatic component from the body force
//
//  By solving for p_rgh, the 6DoF force computation is free of
//  the spurious hydrostatic/buoyancy contribution that caused
//  the angular velocity blowup in freeFallingPimpleFoam.
// ============================================================

// Hydrostatic head at cell centres
volScalarField gh
(
    "gh",
    (fallingFrameForce & mesh.C())
);

// Hydrostatic head at face centres
surfaceScalarField ghf
(
    "ghf",
    (fallingFrameForce & mesh.Cf())
);

// ============================================================
//  Normalised velocity magnitude for post-processing
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
