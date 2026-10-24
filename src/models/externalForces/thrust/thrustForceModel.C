/*---------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::externalForceModels::thrustForce

Implementation
    Foam::externalForceModels::thrustForce
\*---------------------------------------------------------------------------*/

#include "thrustForceModel.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace externalForceModels
{
    defineTypeNameAndDebug(thrustForce, 0);
    addToRunTimeSelectionTable(externalForceModel, thrustForce, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::externalForceModels::thrustForce::thrustForce
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
            << "Zero length direction vector provided for thrust '" << name << "'."
            << exit(FatalError);
    }
    direction_ /= magD;

    if (endTime_ < startTime_)
    {
        FatalErrorInFunction
            << "endTime (" << endTime_ << " s) is less than startTime ("
            << startTime_ << " s) for thrust '" << name << "'."
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::bodyLoad Foam::externalForceModels::thrustForce::evaluate
(
    const scalar t,
    const quaternion& /* q */
) const
{
    if (!active_ || t < startTime_ || t > endTime_)
    {
        return bodyLoad();
    }

    const vector F = magnitude_ * direction_;
    const vector moment = applicationPoint_ ^ F;
    return bodyLoad(F, moment);
}

// ************************************************************************* //
