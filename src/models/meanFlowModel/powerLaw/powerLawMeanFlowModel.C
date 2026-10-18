/*---------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::meanFlowModels::powerLaw

Implementation
    Foam::meanFlowModels::powerLaw

\*---------------------------------------------------------------------------*/

#include "powerLawMeanFlowModel.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace meanFlowModels
{
    defineTypeNameAndDebug(powerLaw, 0);
    addToRunTimeSelectionTable(meanFlowModel, powerLaw, dictionary);
}
}

Foam::meanFlowModels::powerLaw::powerLaw(const dictionary& dict)
:
    meanFlowModel(dict),
    groundHeight_(dict.getOrDefault<scalar>("groundHeight", 0.0)),
    referenceHeight_(dict.getOrDefault<scalar>("referenceHeight", 10.0)),
    referenceVelocity_(dict.getOrDefault<vector>("referenceVelocity", vector::zero)),
    exponent_(dict.getOrDefault<scalar>("exponent", 0.14))
{
    if (referenceHeight_ <= groundHeight_)
    {
        FatalErrorInFunction
            << "referenceHeight (" << referenceHeight_
            << " m) must be strictly greater than groundHeight (" << groundHeight_
            << " m) in powerLaw meanFlow model!" << exit(FatalError);
    }
}

Foam::vector Foam::meanFlowModels::powerLaw::velocity(const scalar altitude) const
{
    if (altitude <= groundHeight_)
    {
        return vector::zero;
    }

    const scalar hNumerator   = altitude - groundHeight_;
    const scalar hDenominator = max(referenceHeight_ - groundHeight_, SMALL);
    const scalar ratio        = max(hNumerator / hDenominator, 0.0);

    return referenceVelocity_ * pow(ratio, exponent_);
}
