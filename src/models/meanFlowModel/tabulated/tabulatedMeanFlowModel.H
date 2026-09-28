/*---------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::meanFlowModels::tabulated

Description
    Altitude-tabulated mean flow model with linear vector interpolation.

\*---------------------------------------------------------------------------*/

#ifndef tabulatedMeanFlowModel_H
#define tabulatedMeanFlowModel_H

#include "meanFlowModel.H"
#include "List.H"
#include "Tuple2.H"

namespace Foam
{
namespace meanFlowModels
{

class tabulated
:
    public meanFlowModel
{
    List<Tuple2<scalar, vector>> table_;

public:

    TypeName("tabulated");

    tabulated(const dictionary& dict);
    virtual ~tabulated() = default;

    virtual vector velocity(const scalar altitude) const override;
};

} // End namespace meanFlowModels
} // End namespace Foam

#endif
