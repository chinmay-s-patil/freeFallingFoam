/*---------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::meanFlowModels::constant

Description
    Constant mean atmospheric flow model U_mean(z) = U_ref.

\*---------------------------------------------------------------------------*/

#ifndef constantMeanFlowModel_H
#define constantMeanFlowModel_H

#include "meanFlowModel.H"

namespace Foam
{
namespace meanFlowModels
{

class constant
:
    public meanFlowModel
{
    vector velocity_;

public:

    TypeName("constant");

    constant(const dictionary& dict);
    virtual ~constant() = default;

    virtual vector velocity(const scalar altitude) const override;
};

} // End namespace meanFlowModels
} // End namespace Foam

#endif
