/*---------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::functionObjects::fallingMotion

Implementation
    Foam::functionObjects::fallingMotion

\*---------------------------------------------------------------------------*/

#include "fallingMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "Pstream.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fallingMotion, 0);
    addToRunTimeSelectionTable(functionObject, fallingMotion, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::fallingMotion::openFile()
{
    if (Pstream::master() && writeCSV_ && !motionFilePtr_.valid())
    {
        fileName outputDir = time_.path()/"postProcessing"/"fallingMotion";
        mkDir(outputDir);

        const fileName motionPath = outputDir/"motion.dat";
        const bool append = (time_.timeIndex() > time_.startTimeIndex()) || (isFile(motionPath) && time_.value() > 0);

        motionFilePtr_.reset
        (
            new OFstream
            (
                motionPath,
                IOstreamOption(),
                (append ? IOstreamOption::APPEND : IOstreamOption::NO_APPEND)
            )
        );

        if (!append)
        {
            motionFilePtr_()
                << "# Time Altitude Vx Vy Vz Ax Ay Az g_x g_y g_z rho a_sound Mach "
                << "Umean_x Umean_y Umean_z Urel_x Urel_y Urel_z "
                << "Fpress_x Fpress_y Fpress_z Fvisc_x Fvisc_y Fvisc_z "
                << "Faero_x Faero_y Faero_z Fgrav_x Fgrav_y Fgrav_z "
                << "Fnet_x Fnet_y Fnet_z Fext_x Fext_y Fext_z" << endl;
        }
    }
}

void Foam::functionObjects::fallingMotion::printConsole(const fvMesh& mesh)
{
    if (!mesh.foundObject<uniformDimensionedVectorField>("frameVelocity"))
    {
        return;
    }

    const vector& Vframe = mesh.lookupObject<uniformDimensionedVectorField>("frameVelocity").value();
    const vector& Aframe = mesh.lookupObject<uniformDimensionedVectorField>("frameAcceleration").value();
    const scalar alt     = mesh.lookupObject<uniformDimensionedScalarField>("altitude").value();
    const vector& Fpress = mesh.lookupObject<uniformDimensionedVectorField>("pressureForce").value();
    const vector& Fvisc  = mesh.lookupObject<uniformDimensionedVectorField>("viscousForce").value();
    const vector& Faero  = mesh.lookupObject<uniformDimensionedVectorField>("aeroForce").value();
    const vector& Fgrav  = mesh.lookupObject<uniformDimensionedVectorField>("gravityForce").value();
    const vector& Fnet   = mesh.lookupObject<uniformDimensionedVectorField>("netForce").value();

    Info<< "\n[Falling Motion] t = " << time_.timeName() << " s:" << endl;

    if (printFrameVel_)
    {
        Info<< "  V = " << mag(Vframe) << " m/s (vector: " << Vframe << ")" << endl;
    }
    if (printFrameAcc_)
    {
        Info<< "  a = " << mag(Aframe) << " m/s2 (vector: " << Aframe << ")" << endl;
    }
    if (printAltitude_)
    {
        Info<< "  Altitude = " << alt << " m" << endl;
    }
    if (printMeanFlow_ && mesh.foundObject<uniformDimensionedVectorField>("meanFlowVelocity"))
    {
        const vector& Umean = mesh.lookupObject<uniformDimensionedVectorField>("meanFlowVelocity").value();
        Info<< "  Umean = " << Umean << " m/s" << endl;
    }
    if (printRelativeVel_ && mesh.foundObject<uniformDimensionedVectorField>("relativeVelocity"))
    {
        const vector& Urel = mesh.lookupObject<uniformDimensionedVectorField>("relativeVelocity").value();
        Info<< "  Urelative = " << Urel << " m/s" << endl;
    }
    if (printPressureForce_)
    {
        Info<< "  Fpressure = " << Fpress << " N" << endl;
    }
    if (printViscousForce_)
    {
        Info<< "  Fviscous = " << Fvisc << " N" << endl;
    }
    if (printAeroForce_)
    {
        Info<< "  Faero = " << Faero << " N" << endl;
    }
    if (printGravityForce_)
    {
        Info<< "  Fgravity = " << Fgrav << " N" << endl;
    }
    if (printNetForce_)
    {
        Info<< "  Fnet = " << Fnet << " N" << endl;
    }
    if (printExternalForce_ && mesh.foundObject<uniformDimensionedVectorField>("externalForce"))
    {
        const vector& Fext = mesh.lookupObject<uniformDimensionedVectorField>("externalForce").value();
        Info<< "  Fexternal = " << Fext << " N" << endl;
    }
    if (printSpeedOfSound_ && mesh.foundObject<uniformDimensionedScalarField>("speedOfSound"))
    {
        Info<< "  Speed of Sound = " << mesh.lookupObject<uniformDimensionedScalarField>("speedOfSound").value() << " m/s" << endl;
    }
    if (printMachNumber_ && mesh.foundObject<uniformDimensionedScalarField>("machNumber"))
    {
        Info<< "  Mach Number = " << mesh.lookupObject<uniformDimensionedScalarField>("machNumber").value() << endl;
    }
    if (printMaxLocalMachNumber_ && mesh.foundObject<uniformDimensionedScalarField>("maxLocalMachNumber"))
    {
        Info<< "  Max Local Mach = " << mesh.lookupObject<uniformDimensionedScalarField>("maxLocalMachNumber").value() << endl;
    }
    if (printDynamicPressure_ && mesh.foundObject<uniformDimensionedScalarField>("rhoAtmosphere") && mesh.foundObject<uniformDimensionedVectorField>("relativeVelocity"))
    {
        const scalar rho = mesh.lookupObject<uniformDimensionedScalarField>("rhoAtmosphere").value();
        const vector& Urel = mesh.lookupObject<uniformDimensionedVectorField>("relativeVelocity").value();
        Info<< "  Dynamic Pressure = " << (0.5 * rho * magSqr(Urel)) << " Pa" << endl;
    }
}

void Foam::functionObjects::fallingMotion::writeLog(const fvMesh& mesh)
{
    if (!Pstream::master() || !writeCSV_)
    {
        return;
    }

    if (!mesh.foundObject<uniformDimensionedVectorField>("frameVelocity"))
    {
        return;
    }

    openFile();

    const vector& Vframe = mesh.lookupObject<uniformDimensionedVectorField>("frameVelocity").value();
    const vector& Aframe = mesh.lookupObject<uniformDimensionedVectorField>("frameAcceleration").value();
    const scalar alt     = mesh.lookupObject<uniformDimensionedScalarField>("altitude").value();
    const vector& grav   = mesh.lookupObject<uniformDimensionedVectorField>("gravity").value();
    const scalar rho     = mesh.lookupObject<uniformDimensionedScalarField>("rhoAtmosphere").value();
    const scalar a_sound = mesh.foundObject<uniformDimensionedScalarField>("speedOfSound")
                         ? mesh.lookupObject<uniformDimensionedScalarField>("speedOfSound").value()
                         : 0.0;
    const scalar Mach    = mesh.foundObject<uniformDimensionedScalarField>("machNumber")
                         ? mesh.lookupObject<uniformDimensionedScalarField>("machNumber").value()
                         : 0.0;
    const vector Umean   = mesh.foundObject<uniformDimensionedVectorField>("meanFlowVelocity")
                         ? mesh.lookupObject<uniformDimensionedVectorField>("meanFlowVelocity").value()
                         : vector::zero;
    const vector Urel    = mesh.foundObject<uniformDimensionedVectorField>("relativeVelocity")
                         ? mesh.lookupObject<uniformDimensionedVectorField>("relativeVelocity").value()
                         : -Vframe;
    const vector& Fpress = mesh.lookupObject<uniformDimensionedVectorField>("pressureForce").value();
    const vector& Fvisc  = mesh.lookupObject<uniformDimensionedVectorField>("viscousForce").value();
    const vector& Faero  = mesh.lookupObject<uniformDimensionedVectorField>("aeroForce").value();
    const vector& Fgrav  = mesh.lookupObject<uniformDimensionedVectorField>("gravityForce").value();
    const vector& Fnet   = mesh.lookupObject<uniformDimensionedVectorField>("netForce").value();
    const vector Fext    = mesh.foundObject<uniformDimensionedVectorField>("externalForce")
                         ? mesh.lookupObject<uniformDimensionedVectorField>("externalForce").value()
                         : vector::zero;

    motionFilePtr_()
        << time_.timeName() << " "
        << alt << " "
        << Vframe.x() << " " << Vframe.y() << " " << Vframe.z() << " "
        << Aframe.x() << " " << Aframe.y() << " " << Aframe.z() << " "
        << grav.x() << " " << grav.y() << " " << grav.z() << " "
        << rho << " "
        << a_sound << " "
        << Mach << " "
        << Umean.x() << " " << Umean.y() << " " << Umean.z() << " "
        << Urel.x() << " " << Urel.y() << " " << Urel.z() << " "
        << Fpress.x() << " " << Fpress.y() << " " << Fpress.z() << " "
        << Fvisc.x() << " " << Fvisc.y() << " " << Fvisc.z() << " "
        << Faero.x() << " " << Faero.y() << " " << Faero.z() << " "
        << Fgrav.x() << " " << Fgrav.y() << " " << Fgrav.z() << " "
        << Fnet.x() << " " << Fnet.y() << " " << Fnet.z() << " "
        << Fext.x() << " " << Fext.y() << " " << Fext.z() << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fallingMotion::fallingMotion
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    time_(runTime),
    writeCSV_(true),
    motionFilePtr_(nullptr),
    regionName_(polyMesh::defaultRegion),
    printFrameVel_(true),
    printFrameAcc_(true),
    printAltitude_(true),
    printPressureForce_(false),
    printViscousForce_(false),
    printAeroForce_(true),
    printGravityForce_(true),
    printNetForce_(true),
    printExternalForce_(true),
    printMeanFlow_(true),
    printRelativeVel_(true),
    printSpeedOfSound_(true),
    printMachNumber_(true),
    printMaxLocalMachNumber_(false),
    printDynamicPressure_(false)
{
    read(dict);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fallingMotion::read(const dictionary& dict)
{
    writeCSV_ = dict.getOrDefault<bool>("writeCSV", true);
    regionName_ = dict.getOrDefault<word>("region", polyMesh::defaultRegion);

    if (dict.found("print"))
    {
        const dictionary& printDict = dict.subDict("print");
        printFrameVel_          = printDict.getOrDefault<bool>("frameVel", true);
        printFrameAcc_          = printDict.getOrDefault<bool>("frameAcc", true);
        printAltitude_          = printDict.getOrDefault<bool>("altitude", true);
        printPressureForce_     = printDict.getOrDefault<bool>("pressureForce", false);
        printViscousForce_      = printDict.getOrDefault<bool>("viscousForce", false);
        printAeroForce_         = printDict.getOrDefault<bool>("aeroForce", true);
        printGravityForce_      = printDict.getOrDefault<bool>("gravityForce", true);
        printNetForce_          = printDict.getOrDefault<bool>("netForce", true);
        printExternalForce_     = printDict.getOrDefault<bool>("externalForce", true);
        printMeanFlow_          = printDict.getOrDefault<bool>("meanFlow", true);
        printRelativeVel_       = printDict.getOrDefault<bool>("relativeVel", true);
        printSpeedOfSound_      = printDict.getOrDefault<bool>("speedOfSound", true);
        printMachNumber_        = printDict.getOrDefault<bool>("machNumber", true);
        printMaxLocalMachNumber_= printDict.getOrDefault<bool>("maxLocalMachNumber", false);
        printDynamicPressure_   = printDict.getOrDefault<bool>("dynamicPressure", false);
    }

    return true;
}


bool Foam::functionObjects::fallingMotion::execute()
{
    if (!time_.foundObject<fvMesh>(regionName_))
    {
        return true;
    }
    const fvMesh& mesh = time_.lookupObject<fvMesh>(regionName_);
    printConsole(mesh);
    writeLog(mesh);
    return true;
}

bool Foam::functionObjects::fallingMotion::write()
{
    return true;
}

// ************************************************************************* //
