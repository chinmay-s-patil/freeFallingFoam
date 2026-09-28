/*---------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::meanFlowModels::tabulated

Implementation
    Foam::meanFlowModels::tabulated

\*---------------------------------------------------------------------------*/

#include "tabulatedMeanFlowModel.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace meanFlowModels
{
    defineTypeNameAndDebug(tabulated, 0);
    addToRunTimeSelectionTable(meanFlowModel, tabulated, dictionary);
}
}

Foam::meanFlowModels::tabulated::tabulated(const dictionary& dict)
:
    meanFlowModel(dict),
    table_()
{
    dict.readIfPresent("table", table_);
}

Foam::vector Foam::meanFlowModels::tabulated::velocity(const scalar altitude) const
{
    if (table_.empty())
    {
        return vector::zero;
    }

    if (altitude <= table_.first().first())
    {
        return table_.first().second();
    }
    if (altitude >= table_.last().first())
    {
        return table_.last().second();
    }

    for (label i = 0; i < table_.size() - 1; ++i)
    {
        if (altitude >= table_[i].first() && altitude <= table_[i+1].first())
        {
            const scalar z0 = table_[i].first();
            const scalar z1 = table_[i+1].first();
            const vector& u0 = table_[i].second();
            const vector& u1 = table_[i+1].second();

            const scalar frac = (altitude - z0) / max(z1 - z0, SMALL);
            return u0 + frac * (u1 - u0);
        }
    }

    return table_.last().second();
}
