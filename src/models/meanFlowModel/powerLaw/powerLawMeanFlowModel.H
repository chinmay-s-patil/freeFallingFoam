/*---------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::meanFlowModels::powerLaw

Description
    Power-law atmospheric mean flow model:
    U(z) = U_ref * ((z - z0) / (z_ref - z0))^alpha

\*---------------------------------------------------------------------------*/

#ifndef powerLawMeanFlowModel_H
#define powerLawMeanFlowModel_H

#include "meanFlowModel.H"

namespace Foam
{
namespace meanFlowModels
{

class powerLaw
:
    public meanFlowModel
{
    scalar groundHeight_;
    scalar referenceHeight_;
    vector referenceVelocity_;
    scalar exponent_;

public:

    TypeName("powerLaw");

    powerLaw(const dictionary& dict);
    virtual ~powerLaw() = default;

    virtual vector velocity(const scalar altitude) const override;
};

} // End namespace meanFlowModels
} // End namespace Foam

#endif
