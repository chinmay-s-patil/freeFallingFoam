/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::solidBodyMotionFunctions::sixDoFRotatingMotion

Description
    Implementation of the torque-driven 6DoF axial rotation model.

\*---------------------------------------------------------------------------*/

#include "sixDoFRotatingMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "forces.H"
#include "mathematicalConstants.H"
#include "Pstream.H"
#include "IOdictionary.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(sixDoFRotatingMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        sixDoFRotatingMotion,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solidBodyMotionFunctions::sixDoFRotatingMotion::readState()
{
    IOobject io
    (
        "sixDoFRotatingMotionState",
        time_.timeName(),
        "uniform",
        time_,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    if (io.typeHeaderOk<IOdictionary>(true))
    {
        IOdictionary dict(io);
        dict.readIfPresent("omega", omega_);
        dict.readIfPresent("theta", theta_);
        dict.readIfPresent("torque", torque_);
        if (Pstream::master())
        {
            Info<< "sixDoFRotatingMotion: Restored state at time " << time_.timeName()
                << " -> omega = " << omega_ << " rad/s, theta = " << theta_ << " rad" << endl;
        }
    }
}


void Foam::solidBodyMotionFunctions::sixDoFRotatingMotion::logMotion
(
    const scalar t,
    const scalar torque,
    const scalar omega,
    const scalar rpm,
    const scalar deg
) const
{
    if (!Pstream::master())
    {
        return;
    }

    if (!logFilePtr_)
    {
        fileName dir = time_.globalPath() / "postProcessing" / "sixDoFRotatingMotion";
        mkDir(dir);
        fileName logFile = dir / "motion.dat";

        bool exists = isFile(logFile);
        logFilePtr_.reset(new OFstream(logFile, IOstreamOption(), IOstreamOption::APPEND));

        if (!exists)
        {
            logFilePtr_()
                << "# Time(s)\tTorque(N.m)\tOmega(rad/s)\tRPM\tAngle(deg)" << nl;
        }
    }

    logFilePtr_()
        << t << tab
        << torque << tab
        << omega << tab
        << rpm << tab
        << deg << nl;
    logFilePtr_().flush();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::sixDoFRotatingMotion::sixDoFRotatingMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime),
    origin_(SBMFCoeffs.getOrDefault<point>("origin", point::zero)),
    axis_(SBMFCoeffs.getOrDefault<vector>("axis", vector(0, 1, 0))),
    patches_(SBMFCoeffs.get<wordRes>("patches")),
    mass_(SBMFCoeffs.getOrDefault<scalar>("mass", 0.012)),
    momentOfInertia_(SBMFCoeffs.getOrDefault<scalar>("momentOfInertia", 6.5e-5)),
    omega0_(SBMFCoeffs.getOrDefault<scalar>("omega0", 0.0)),
    damping_(SBMFCoeffs.getOrDefault<scalar>("damping", 0.0)),
    rhoInf_(SBMFCoeffs.getOrDefault<scalar>("rhoInf", 1.225)),
    rhoName_(SBMFCoeffs.getOrDefault<word>("rho", "rhoInf")),
    report_(SBMFCoeffs.getOrDefault<bool>("report", true)),
    omega_(omega0_),
    theta_(0.0),
    torque_(0.0),
    lastTimeIndex_(-1),
    logFilePtr_(nullptr)
{
    axis_.normalise();
    readState();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::sixDoFRotatingMotion::transformation() const
{
    const scalar dt = time_.deltaTValue();
    const label curTimeIndex = time_.timeIndex();

    if (curTimeIndex != lastTimeIndex_ && curTimeIndex > 0)
    {
        scalar aeroTorque = 0.0;

        try
        {
            dictionary forcesDict;
            forcesDict.add("type", "forces");
            forcesDict.add("patches", patches_);
            forcesDict.add("rhoInf", rhoInf_);
            forcesDict.add("rho", rhoName_);
            forcesDict.add("CofR", origin_);
            forcesDict.add("writeToFile", false);
            forcesDict.add("log", false);

            functionObjects::forces f("sixDoFRotorForces", time_, forcesDict, true);
            f.calcForcesMoments();

            vector totalMoment = f.momentEff();
            aeroTorque = (totalMoment & axis_);
        }
        catch (...)
        {
            aeroTorque = 0.0;
        }

        torque_ = aeroTorque;

        if (dt > 1e-15 && momentOfInertia_ > 1e-15)
        {
            // Angular momentum equation:
            // I * (omega^{n+1} - omega^n)/dt = Torque - damping * omega^{n+1}
            // omega^{n+1} * (1 + damping*dt/I) = omega^n + (Torque/I)*dt
            scalar denom = 1.0 + (damping_ / momentOfInertia_) * dt;
            omega_ = (omega_ + (torque_ / momentOfInertia_) * dt) / denom;

            // Orientation angle increment
            theta_ += omega_ * dt;
        }

        lastTimeIndex_ = curTimeIndex;

        scalar rpm = omega_ * 60.0 / (2.0 * constant::mathematical::pi);
        scalar deg = theta_ * 180.0 / constant::mathematical::pi;

        if (report_ && Pstream::master())
        {
            Info<< "sixDoFRotatingMotion:"
                << "  Time = " << time_.timeName() << " s"
                << ", Torque = " << torque_ << " N.m"
                << ", Omega = " << omega_ << " rad/s (" << rpm << " RPM)"
                << ", Angle = " << deg << " deg"
                << endl;
        }

        logMotion(time_.value(), torque_, omega_, rpm, deg);
    }

    quaternion q(axis_, theta_);
    return septernion(-origin_) * q * septernion(origin_);
}


bool Foam::solidBodyMotionFunctions::sixDoFRotatingMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    if (solidBodyMotionFunction::read(SBMFCoeffs))
    {
        SBMFCoeffs.readIfPresent("origin", origin_);
        SBMFCoeffs.readIfPresent("axis", axis_);
        axis_.normalise();
        SBMFCoeffs.readIfPresent("patches", patches_);
        SBMFCoeffs.readIfPresent("mass", mass_);
        SBMFCoeffs.readIfPresent("momentOfInertia", momentOfInertia_);
        SBMFCoeffs.readIfPresent("damping", damping_);
        SBMFCoeffs.readIfPresent("rhoInf", rhoInf_);
        SBMFCoeffs.readIfPresent("report", report_);
        return true;
    }

    return false;
}


void Foam::solidBodyMotionFunctions::sixDoFRotatingMotion::writeData
(
    Ostream& os
) const
{
    solidBodyMotionFunction::writeData(os);

    os.writeEntry("omega", omega_);
    os.writeEntry("theta", theta_);
    os.writeEntry("torque", torque_);

    if (Pstream::master())
    {
        fileName uniformDir = time_.timePath() / "uniform";
        mkDir(uniformDir);
        IOdictionary dict
        (
            IOobject
            (
                "sixDoFRotatingMotionState",
                time_.timeName(),
                "uniform",
                time_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        );
        dict.add("omega", omega_);
        dict.add("theta", theta_);
        dict.add("torque", torque_);
        OFstream stateFile(uniformDir / "sixDoFRotatingMotionState");
        dict.writeHeader(stateFile);
        dict.writeData(stateFile);
    }
}

// ************************************************************************* //
