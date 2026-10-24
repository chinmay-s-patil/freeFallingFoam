/*---------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::externalForceModels::constantForce

Implementation
    Foam::externalForceModels::constantForce
\*---------------------------------------------------------------------------*/

#include "constantForceModel.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace externalForceModels
{
    defineTypeNameAndDebug(constantForce, 0);
    addToRunTimeSelectionTable(externalForceModel, constantForce, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::externalForceModels::constantForce::constantForce
(
    const word& name,
    const dictionary& dict
)
:
    externalForceModel(name, dict),
    force_(dict.get<vector>("force")),
    applicationPoint_(dict.getOrDefault<point>("applicationPoint", point::zero))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::bodyLoad Foam::externalForceModels::constantForce::evaluate
(
    const scalar /* t */,
    const quaternion& /* q */
) const
{
    if (!active_)
    {
        return bodyLoad();
    }

    const vector moment = applicationPoint_ ^ force_;
    return bodyLoad(force_, moment);
}

// ************************************************************************* //
