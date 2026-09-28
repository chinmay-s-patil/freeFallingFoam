/*---------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::meanFlowModels::constant

Implementation
    Foam::meanFlowModels::constant

\*---------------------------------------------------------------------------*/

#include "constantMeanFlowModel.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace meanFlowModels
{
    defineTypeNameAndDebug(constant, 0);
    addToRunTimeSelectionTable(meanFlowModel, constant, dictionary);
}
}

Foam::meanFlowModels::constant::constant(const dictionary& dict)
:
    meanFlowModel(dict),
    velocity_(dict.getOrDefault<vector>("velocity", vector::zero))
{}

Foam::vector Foam::meanFlowModels::constant::velocity(const scalar altitude) const
{
    return velocity_;
}
