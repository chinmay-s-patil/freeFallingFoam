// ============================================================
//  createFallingFrameFields.H — Shared freeFallingFoam field init
// ============================================================

#include "uniformDimensionedFields.H"
#include "meanFlowModel.H"

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
    fallingFrameDict.getOrDefault<vector>("u0", vector::zero)
);

// Support both 'gravity' and 'acceleration' dictionary keys for planetary gravity
const vector initialGravity
(
    fallingFrameDict.found("gravity")
  ? fallingFrameDict.get<vector>("gravity")
  : fallingFrameDict.getOrDefault<vector>("acceleration", vector(0, -9.81, 0))
);

// Motion degrees of freedom mask (vector): 1 = free, 0 = constrained
const vector freedomMask
(
    fallingFrameDict.getOrDefault<vector>("freedom", vector(0, 1, 0))
);

if (freedomMask.x() < 0.0 || freedomMask.x() > 1.0 ||
    freedomMask.y() < 0.0 || freedomMask.y() > 1.0 ||
    freedomMask.z() < 0.0 || freedomMask.z() > 1.0)
{
    Info<< "  [Warning] freedomMask has non-binary entries: " << freedomMask
        << ". Acceleration will be scaled accordingly." << nl;
}

// Initial kinematic frame acceleration (defaults to gravity masked by freedom)
const vector initialFrameAcc
(
    fallingFrameDict.found("initialAcceleration")
  ? cmptMultiply(fallingFrameDict.get<vector>("initialAcceleration"), freedomMask)
  : cmptMultiply(initialGravity, freedomMask)
);

const scalar bodyMass
(
    fallingFrameDict.get<scalar>("mass")
);

if (bodyMass <= 0.0)
{
    FatalErrorInFunction
        << "bodyMass must be strictly positive (> 0), but found: "
        << bodyMass << exit(FatalError);
}

const wordList patchNames
(
    fallingFrameDict.get<wordList>("patches")
);

// Validate that specified patches exist in the mesh
forAll(patchNames, patchI)
{
    const label patchID = mesh.boundaryMesh().findPatchID(patchNames[patchI]);
    if (patchID == -1)
    {
        FatalErrorInFunction
            << "Patch '" << patchNames[patchI] << "' specified in fallingFrameDict "
            << "was not found in boundary mesh! Available patches: "
            << mesh.boundaryMesh().names()
            << exit(FatalError);
    }
}

const scalar rhoRef
(
    fallingFrameDict.getOrDefault<scalar>("rhoInf", 1.204)
);

if (rhoRef <= 0.0)
{
    FatalErrorInFunction
        << "rhoInf must be strictly positive (> 0), but found: "
        << rhoRef << exit(FatalError);
}

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

if (groundInterrupt && (initialAltitude <= groundAltitude))
{
    FatalErrorInFunction
        << "groundInterrupt is enabled, but initialAltitude (" << initialAltitude
        << " m) is <= groundAltitude (" << groundAltitude
        << " m). Case would terminate immediately on startup!"
        << exit(FatalError);
}

const scalar terminalVelocityTolerance
(
    fallingFrameDict.getOrDefault<scalar>("terminalVelocityTolerance", 1e-3)
);

const scalar terminalAccelerationTolerance
(
    fallingFrameDict.getOrDefault<scalar>("terminalAccelerationTolerance", 1e-3)
);

if (terminalVelocityTolerance <= 0.0 || terminalAccelerationTolerance <= 0.0)
{
    FatalErrorInFunction
        << "terminalVelocityTolerance and terminalAccelerationTolerance must be positive (> 0), but found: "
        << terminalVelocityTolerance << ", " << terminalAccelerationTolerance
        << exit(FatalError);
}

#if defined(FF_INCOMPRESSIBLE)
const bool includeBuoyancy
(
    fallingFrameDict.getOrDefault<bool>("includeBuoyancy", true)
);
#endif

label terminalVelocitySamples
(
    fallingFrameDict.getOrDefault<label>("terminalVelocitySamples", 20)
);
if (terminalVelocitySamples < 2)
{
    terminalVelocitySamples = 2;
}

const bool terminalVelocityWrite
(
    fallingFrameDict.getOrDefault<bool>("terminalVelocityWrite", true)
);

const bool terminalVelocityStop
(
    fallingFrameDict.getOrDefault<bool>("terminalVelocityStop", false)
);

// Configurable environmental and physical constants
const scalar scaleHeight
(
    fallingFrameDict.getOrDefault<scalar>("scaleHeight", 8400.0)
);

const scalar earthRadius
(
    fallingFrameDict.getOrDefault<scalar>("earthRadius", 6371000.0)
);

#if defined(FF_INCOMPRESSIBLE)
const scalar gammaCoeff
(
    fallingFrameDict.getOrDefault<scalar>("gamma", 1.4)
);

const scalar pAtmSeaLevel
(
    fallingFrameDict.getOrDefault<scalar>("seaLevelPressure", 101325.0)
);
#endif

const scalar speedOfSoundRef
(
    fallingFrameDict.getOrDefault<scalar>("speedOfSound", 343.0)
);

// Construct mean flow model once during initialization
autoPtr<meanFlowModel> meanFlowPtr;
if (fallingFrameDict.found("meanFlow"))
{
    meanFlowPtr = meanFlowModel::New(fallingFrameDict);
}

bool triggerTerminalVelocityWrite = false;
bool triggerGroundInterrupt = false;
bool triggerMachStopper = false;

Info<< "  Falling body properties:" << nl
    << "    u0             = " << u0 << " m/s" << nl
    << "    initialAcc     = " << initialFrameAcc << " m/s2" << nl
    << "    gravity        = " << initialGravity << " m/s2" << nl
    << "    freedom        = " << freedomMask << nl
    << "    mass           = " << bodyMass << " kg" << nl
    << "    patches        = " << patchNames << nl
    << "    rhoInf         = " << rhoRef << " kg/m3" << nl
    << "    densityModel   = " << densityModel << nl
    << "    gravityModel   = " << gravityModel << nl
    << "    groundInterrupt= " << (groundInterrupt ? "on" : "off") << nl
    << "    initialAltitude= " << initialAltitude << " m" << nl
    << "    groundAltitude = " << groundAltitude << " m" << nl
    << "    termVelSamples = " << terminalVelocitySamples << nl
    << "    termVelWrite   = " << (terminalVelocityWrite ? "on" : "off") << nl
    << "    termVelStop    = " << (terminalVelocityStop ? "on" : "off") << nl;

#if defined(FF_COMPRESSIBLE)
if (densityModel == "altitudeBased")
{
    Info<< "    Note: In compressible solvers, flow density is governed by the" << nl
        << "          thermodynamic equation of state. rhoAtmosphere tracks" << nl
        << "          the barometric reference density for diagnostics." << nl;
}
#endif

#if defined(FF_PIMPLE_COMPRESSIBLE)
const word machRegime(fallingFrameDict.getOrDefault<word>("machRegime", "subsonic"));
const scalar machRegimeSwitchCutOff
(
    fallingFrameDict.found("machCutOff")
  ? fallingFrameDict.get<scalar>("machCutOff")
  : fallingFrameDict.getOrDefault<scalar>("machRegimeSwitchCutOff", 0.7)
);
const word machRegimeAction
(
    fallingFrameDict.found("machAction")
  ? fallingFrameDict.get<word>("machAction")
  : fallingFrameDict.getOrDefault<word>("machRegimeAction", "switch")
);
bool machRegimeSwitched = false;

Info<< "    machRegime     = " << machRegime << nl
    << "    machCutOff     = " << machRegimeSwitchCutOff << nl
    << "    machAction     = " << machRegimeAction << nl;

if (machRegime != "subsonic" && machRegime != "transonic" &&
    machRegime != "auto" && machRegime != "autoSwitch" &&
    machRegime != "on" && machRegime != "off")
{
    FatalErrorInFunction
        << "Unknown machRegime: '" << machRegime
        << "'. Valid options are: subsonic, transonic, auto, autoSwitch, on, off."
        << exit(FatalError);
}

if ((machRegime == "auto" || machRegime == "autoSwitch") &&
    machRegimeAction != "switch" && machRegimeAction != "stop")
{
    FatalErrorInFunction
        << "Unknown machAction: '" << machRegimeAction
        << "'. Valid options are: switch, stop."
        << exit(FatalError);
}

if (mesh.foundObject<IOdictionary>("fvSolution"))
{
    IOdictionary& solDict = mesh.lookupObjectRef<IOdictionary>("fvSolution");
    if (solDict.found("PIMPLE"))
    {
        dictionary& pimpleDict = solDict.subDict("PIMPLE");
        if (machRegime == "transonic" || machRegime == "on")
        {
            pimpleDict.set("transonic", true);
            machRegimeSwitched = true;
        }
        else if (machRegime == "subsonic" || machRegime == "off")
        {
            pimpleDict.set("transonic", false);
        }
        else if (machRegime == "auto" || machRegime == "autoSwitch")
        {
            // On restart, preserve existing transonic setting from fvSolution
            if (runTime.timeIndex() > runTime.startTimeIndex() || runTime.value() > 0)
            {
                machRegimeSwitched = pimpleDict.getOrDefault<bool>("transonic", false);
            }
            else
            {
                pimpleDict.set("transonic", false);
            }
        }
    }
}
#endif

Info<< nl;


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
    dimensionedVector("frameAcceleration", dimAcceleration, initialFrameAcc)
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
    dimensionedVector("gravity", dimAcceleration, initialGravity)
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
    dimensionedVector("gravityForce", dimForce, bodyMass * initialGravity)
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
    dimensionedVector("netForce", dimForce, bodyMass * initialFrameAcc)
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

uniformDimensionedScalarField speedOfSound
(
    IOobject
    (
        "speedOfSound",
        runTime.timeName(),
        "uniform",
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    dimensionedScalar("speedOfSound", dimVelocity, speedOfSoundRef)
);

uniformDimensionedScalarField machNumber
(
    IOobject
    (
        "machNumber",
        runTime.timeName(),
        "uniform",
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    dimensionedScalar("machNumber", dimless, 0.0)
);

#if defined(FF_PIMPLE_COMPRESSIBLE)
uniformDimensionedScalarField maxLocalMachNumber
(
    IOobject
    (
        "maxLocalMachNumber",
        runTime.timeName(),
        "uniform",
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    dimensionedScalar("maxLocalMachNumber", dimless, 0.0)
);
#endif


// ============================================================
//  Spatial Fields
// ============================================================

// Apparent body force in the non-inertial body frame: g - a_frame
dimensionedVector fallingFrameForce
(
    "fallingFrameForce",
    dimAcceleration,
    gravity.value() - frameAcceleration.value()
);

#if defined(FF_INCOMPRESSIBLE)
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
#endif

// Dimensionless normalized velocity magnitude |U| / max(|V_frame|, 1.0)
// used for non-dimensional flow field scaling and post-processing visualization
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


// ============================================================
//  State & History for Predictor-Corrector and Diagnostics
// ============================================================

// O(1) circular ring buffer for natural terminal velocity detection
List<scalar> velHistory(terminalVelocitySamples, 0.0);
List<vector> accHistory(terminalVelocitySamples, vector::zero);
label velHistoryCount = 0;
label velHistoryIndex = 0;
bool terminalVelocityTriggered = false;

// Store state at time n for 2nd-order Heun predictor-corrector
vector a_prev = frameAcceleration.value();
vector v_prev = frameVelocity.value();
scalar alt_prev = altitude.value();
