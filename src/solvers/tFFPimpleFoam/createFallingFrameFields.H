// ============================================================
//  createFallingFrameFields.H — tFFPimpleFoam
// ============================================================

#include "uniformDimensionedFields.H"

Info<< "Reading fallingFrameDict\n" << endl;

IOdictionary fallingFrameDict
(
    IOobject
    (
        "fallingFrameDict",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const vector u0
(
    fallingFrameDict.get<vector>("u0")
);

const vector initialAcc
(
    fallingFrameDict.get<vector>("acceleration")
);

const scalar bodyMass
(
    fallingFrameDict.get<scalar>("mass")
);

const wordList patchNames
(
    fallingFrameDict.get<wordList>("patches")
);

const scalar rhoRef
(
    fallingFrameDict.getOrDefault<scalar>("rhoInf", 1.204)
);

const word densityModel
(
    fallingFrameDict.getOrDefault<word>("densityModel", "rhoInf")
);

const word gravityModel
(
    fallingFrameDict.getOrDefault<word>("gravityModel", "constant")
);

const bool groundInterrupt
(
    fallingFrameDict.getOrDefault<bool>("groundInterrupt", false)
);

const scalar initialAltitude
(
    fallingFrameDict.getOrDefault<scalar>("initialAltitude", 0.0)
);

const scalar groundAltitude
(
    fallingFrameDict.getOrDefault<scalar>("groundAltitude", 0.0)
);

const scalar terminalVelocityTolerance
(
    fallingFrameDict.getOrDefault<scalar>("terminalVelocityTolerance", 1e-3)
);

const scalar terminalAccelerationTolerance
(
    fallingFrameDict.getOrDefault<scalar>("terminalAccelerationTolerance", 1e-3)
);

const vector freedomMask
(
    fallingFrameDict.getOrDefault<vector>("freedom", vector(0, 1, 0))
);

const label terminalVelocitySamples
(
    fallingFrameDict.getOrDefault<label>("terminalVelocitySamples", 20)
);

Info<< "  Falling body properties:" << nl
    << "    u0             = " << u0 << " m/s" << nl
    << "    acceleration   = " << initialAcc << " m/s2" << nl
    << "    freedom        = " << freedomMask << nl
    << "    mass           = " << bodyMass << " kg" << nl
    << "    patches        = " << patchNames << nl
    << "    rhoInf         = " << rhoRef << " kg/m3" << nl
    << "    densityModel   = " << densityModel << nl
    << "    gravityModel   = " << gravityModel << nl
    << "    groundInterrupt= " << (groundInterrupt ? "on" : "off") << nl
    << "    initialAltitude= " << initialAltitude << " m" << nl
    << "    groundAltitude = " << groundAltitude << " m" << nl
    << nl;

// ============================================================
//  Registered Uniform Fields (stored in db() and written to uniform/)
// ============================================================

uniformDimensionedVectorField frameVelocity
(
    IOobject
    (
        "frameVelocity",
        runTime.timeName(),
        "uniform",
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    dimensionedVector("frameVelocity", dimVelocity, u0)
);

uniformDimensionedVectorField frameAcceleration
(
    IOobject
    (
        "frameAcceleration",
        runTime.timeName(),
        "uniform",
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    dimensionedVector("frameAcceleration", dimAcceleration, initialAcc)
);

uniformDimensionedScalarField altitude
(
    IOobject
    (
        "altitude",
        runTime.timeName(),
        "uniform",
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    dimensionedScalar("altitude", dimLength, initialAltitude)
);

uniformDimensionedVectorField gravity
(
    IOobject
    (
        "gravity",
        runTime.timeName(),
        "uniform",
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    dimensionedVector("gravity", dimAcceleration, initialAcc)
);

uniformDimensionedScalarField rhoAtmosphere
(
    IOobject
    (
        "rhoAtmosphere",
        runTime.timeName(),
        "uniform",
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    dimensionedScalar("rhoAtmosphere", dimDensity, rhoRef)
);

uniformDimensionedVectorField pressureForce
(
    IOobject
    (
        "pressureForce",
        runTime.timeName(),
        "uniform",
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    dimensionedVector("pressureForce", dimForce, vector::zero)
);

uniformDimensionedVectorField viscousForce
(
    IOobject
    (
        "viscousForce",
        runTime.timeName(),
        "uniform",
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    dimensionedVector("viscousForce", dimForce, vector::zero)
);

uniformDimensionedVectorField aeroForce
(
    IOobject
    (
        "aeroForce",
        runTime.timeName(),
        "uniform",
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    dimensionedVector("aeroForce", dimForce, vector::zero)
);

uniformDimensionedVectorField gravityForce
(
    IOobject
    (
        "gravityForce",
        runTime.timeName(),
        "uniform",
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    dimensionedVector("gravityForce", dimForce, bodyMass * initialAcc)
);

uniformDimensionedVectorField netForce
(
    IOobject
    (
        "netForce",
        runTime.timeName(),
        "uniform",
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    dimensionedVector("netForce", dimForce, bodyMass * initialAcc)
);

uniformDimensionedVectorField meanFlowVelocity
(
    IOobject
    (
        "meanFlowVelocity",
        runTime.timeName(),
        "uniform",
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    dimensionedVector("meanFlowVelocity", dimVelocity, vector::zero)
);

uniformDimensionedVectorField relativeVelocity
(
    IOobject
    (
        "relativeVelocity",
        runTime.timeName(),
        "uniform",
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    dimensionedVector("relativeVelocity", dimVelocity, -u0)
);


// ============================================================
//  Spatial Fields
// ============================================================

dimensionedVector fallingFrameForce
(
    "fallingFrameForce",
    dimAcceleration,
    -frameAcceleration.value()
);

volScalarField gh
(
    "gh",
    (fallingFrameForce & mesh.C())
);

surfaceScalarField ghf
(
    "ghf",
    (fallingFrameForce & mesh.Cf())
);

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
