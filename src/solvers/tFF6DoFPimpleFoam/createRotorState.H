// ============================================================
//  createRotorState.H — 6DoF Rotational State Initialization
// ============================================================

#include "uniformDimensionedFields.H"
#include "quaternion.H"
#include "mathematicalConstants.H"
#include "OFstream.H"

Info<< "Reading rotor configuration from fallingFrameDict\n" << endl;

const dictionary& rotorDict =
    fallingFrameDict.found("rotor")
  ? fallingFrameDict.subDict("rotor")
  : fallingFrameDict;

const point hubOrigin
(
    rotorDict.getOrDefault<point>("hubOrigin", point::zero)
);

const vector Ib
(
    rotorDict.found("inertiaPrincipal")
  ? rotorDict.get<vector>("inertiaPrincipal")
  : rotorDict.getOrDefault<vector>("momentOfInertia", vector(3.919575e-05, 3.945599e-05, 3.769299e-07))
);

if (Ib.x() <= 0.0 || Ib.y() <= 0.0 || Ib.z() <= 0.0)
{
    FatalErrorInFunction
        << "Principal moments of inertia must be strictly positive: "
        << Ib << exit(FatalError);
}

// Body inertia tensor in solver/mesh frame
tensor JbCalc(Ib.x(), 0, 0,  0, Ib.y(), 0,  0, 0, Ib.z());

if (rotorDict.found("meshTiltAngleDeg"))
{
    const scalar tiltAngleDeg = rotorDict.get<scalar>("meshTiltAngleDeg");
    const vector tiltAxis = rotorDict.getOrDefault<vector>("meshTiltAxis", vector(1, 0, 0));
    const tensor R0 = quaternion(normalised(tiltAxis), tiltAngleDeg * constant::mathematical::pi / 180.0).R();
    JbCalc = R0 & JbCalc & R0.T();
    Info<< "  Mesh pre-tilt angle: " << tiltAngleDeg << " deg about axis " << tiltAxis << endl;
}
else if (rotorDict.found("meshTiltOrientation"))
{
    const tensor R0(rotorDict.get<tensor>("meshTiltOrientation"));
    JbCalc = R0 & JbCalc & R0.T();
}
else if (rotorDict.found("inertiaTensor"))
{
    JbCalc = rotorDict.get<tensor>("inertiaTensor");
}

const tensor Jb = JbCalc;
Info<< "  Body inertia tensor in solver/mesh frame: " << Jb << endl;
const tensor JbInv = inv(Jb);

const vector rotFreedom
(
    rotorDict.getOrDefault<vector>("rotationFreedom", vector::one)
);

const vector rotorDamping
(
    rotorDict.getOrDefault<vector>("damping", vector::zero)
);

const bool reportRotor
(
    rotorDict.getOrDefault<bool>("report", true)
);

// Initial angular velocity in body frame [rad/s]
vector omegaB
(
    rotorDict.getOrDefault<vector>("initialOmega", vector::zero)
);

// Initial attitude quaternion (body -> world)
quaternion qBody = quaternion::I;

if (rotorDict.found("initialOrientation"))
{
    tensor rotTensor(rotorDict.get<tensor>("initialOrientation"));
    qBody = quaternion(rotTensor);
    qBody.normalise();
}
else if (rotorDict.found("initialAngleDeg"))
{
    const scalar angleDeg = rotorDict.get<scalar>("initialAngleDeg");
    const vector axis = rotorDict.getOrDefault<vector>("initialAxis", vector(1, 0, 0));
    const scalar angleRad = angleDeg * constant::mathematical::pi / 180.0;
    qBody = quaternion(normalised(axis), angleRad);
    qBody.normalise();
}

// Check for restart fields in the current time directory
uniformDimensionedScalarField attitudeW
(
    IOobject
    (
        "attitudeW",
        runTime.timeName(),
        "uniform",
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    dimensionedScalar("attitudeW", dimless, qBody.w())
);

uniformDimensionedVectorField attitudeVec
(
    IOobject
    (
        "attitudeVec",
        runTime.timeName(),
        "uniform",
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    dimensionedVector("attitudeVec", dimless, qBody.v())
);

// If restart file was read, reconstruct qBody
if (attitudeW.headerOk() && attitudeVec.headerOk())
{
    qBody = quaternion(attitudeW.value(), attitudeVec.value());
    qBody.normalise();
    Info<< "  [Restart] Reconstructed attitude quaternion from time "
        << runTime.timeName() << ": w=" << qBody.w() << " v=" << qBody.v() << endl;
}

uniformDimensionedVectorField angularVelocity
(
    IOobject
    (
        "angularVelocity",
        runTime.timeName(),
        "uniform",
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    dimensionedVector("angularVelocity", inv(dimTime), omegaB)
);

if (angularVelocity.headerOk())
{
    omegaB = angularVelocity.value();
    Info<< "  [Restart] Restored angular velocity from time "
        << runTime.timeName() << ": " << omegaB << " rad/s" << endl;
}

uniformDimensionedVectorField aeroMoment
(
    IOobject
    (
        "aeroMoment",
        runTime.timeName(),
        "uniform",
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    dimensionedVector("aeroMoment", dimForce*dimLength, vector::zero)
);

uniformDimensionedVectorField eulerAngleDeg
(
    IOobject
    (
        "eulerAngleDeg",
        runTime.timeName(),
        "uniform",
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    dimensionedVector
    (
        "eulerAngleDeg",
        dimless,
        qBody.eulerAngles(quaternion::XYZ) * (180.0 / constant::mathematical::pi)
    )
);

// Predictor state holders
quaternion qBody_old = qBody;
vector omegaB_old = omegaB;

// Initialize log file for tumble tracking
autoPtr<OFstream> tumbleLogPtr;
if (Pstream::master())
{
    fileName tumbleDir = runTime.globalPath() / "postProcessing" / "tumbleMotion";
    mkDir(tumbleDir);
    fileName tumbleFile = tumbleDir / "motion.dat";

    bool exists = isFile(tumbleFile);
    tumbleLogPtr.reset(new OFstream(tumbleFile, IOstreamOption(), IOstreamOption::APPEND));

    if (!exists)
    {
        tumbleLogPtr()
            << "# Time(s)\t"
            << "Omega_x(rad/s)\tOmega_y(rad/s)\tOmega_z(rad/s)\t"
            << "RPM_x\tRPM_y\tRPM_z\t"
            << "Moment_x(N.m)\tMoment_y(N.m)\tMoment_z(N.m)\t"
            << "Angle_x(deg)\tAngle_y(deg)\tAngle_z(deg)\t"
            << "Quat_w\tQuat_x\tQuat_y\tQuat_z" << nl;
        tumbleLogPtr().flush();
    }
}
