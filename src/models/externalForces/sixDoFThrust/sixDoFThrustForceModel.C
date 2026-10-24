/*---------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::externalForceModels::sixDoFThrustForce

Implementation
    Foam::externalForceModels::sixDoFThrustForce
\*---------------------------------------------------------------------------*/

#include "sixDoFThrustForceModel.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace externalForceModels
{
    defineTypeNameAndDebug(sixDoFThrustForce, 0);
    addToRunTimeSelectionTable(externalForceModel, sixDoFThrustForce, dictionary);
    addNamedToRunTimeSelectionTable(externalForceModel, sixDoFThrustForce, dictionary, sixDoFThrust);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::externalForceModels::sixDoFThrustForce::sixDoFThrustForce
(
    const word& name,
    const dictionary& dict
)
:
    externalForceModel(name, dict),
    direction_(dict.get<vector>("direction")),
    magnitude_(dict.get<scalar>("magnitude")),
    startTime_(dict.getOrDefault<scalar>("startTime", 0.0)),
    endTime_(dict.getOrDefault<scalar>("endTime", GREAT)),
    applicationPoint_(dict.getOrDefault<point>("applicationPoint", point::zero))
{
    const scalar magD = mag(direction_);
    if (magD < SMALL)
    {
        FatalErrorInFunction
            << "Zero length direction vector provided for 6DoFThrust '" << name << "'."
            << exit(FatalError);
    }
    direction_ /= magD;

    if (endTime_ < startTime_)
    {
        FatalErrorInFunction
            << "endTime (" << endTime_ << " s) is less than startTime ("
            << startTime_ << " s) for 6DoFThrust '" << name << "'."
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::bodyLoad Foam::externalForceModels::sixDoFThrustForce::evaluate
(
    const scalar t,
    const quaternion& q
) const
{
    if (!active_ || t < startTime_ || t > endTime_)
    {
        return bodyLoad();
    }

    // Force in body coordinates
    const vector F_body = magnitude_ * direction_;

    // Transform force and application point from body to global coordinates
    const vector F_global = q.transform(F_body);
    const vector r_global = q.transform(applicationPoint_);

    const vector moment_global = r_global ^ F_global;
    return bodyLoad(F_global, moment_global);
}

// ************************************************************************* //
