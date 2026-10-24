/*---------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::externalForceModels::tabulatedForce

Implementation
    Foam::externalForceModels::tabulatedForce
\*---------------------------------------------------------------------------*/

#include "tabulatedForceModel.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace externalForceModels
{
    defineTypeNameAndDebug(tabulatedForce, 0);
    addToRunTimeSelectionTable(externalForceModel, tabulatedForce, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::externalForceModels::tabulatedForce::tabulatedForce
(
    const word& name,
    const dictionary& dict
)
:
    externalForceModel(name, dict),
    table_(),
    outOfRange_(dict.getOrDefault<word>("outOfRange", "zero")),
    interpolation_(dict.getOrDefault<word>("interpolation", "linear")),
    applicationPoint_(dict.getOrDefault<point>("applicationPoint", point::zero))
{
    dict.readEntry("table", table_);

    if (table_.empty())
    {
        FatalErrorInFunction
            << "Empty 'table' specified for tabulatedForce '" << name << "'."
            << exit(FatalError);
    }

    if (outOfRange_ != "zero" && outOfRange_ != "hold" && outOfRange_ != "clamp")
    {
        FatalErrorInFunction
            << "Invalid 'outOfRange' setting: '" << outOfRange_
            << "' for tabulatedForce '" << name
            << "'. Valid options are: zero, hold."
            << exit(FatalError);
    }

    if (interpolation_ != "linear")
    {
        FatalErrorInFunction
            << "Invalid 'interpolation' setting: '" << interpolation_
            << "' for tabulatedForce '" << name
            << "'. Currently only 'linear' is supported."
            << exit(FatalError);
    }

    // Verify strictly increasing time entries
    for (label i = 1; i < table_.size(); ++i)
    {
        if (table_[i].first() <= table_[i-1].first())
        {
            FatalErrorInFunction
                << "Time entries in table for tabulatedForce '" << name
                << "' must be strictly increasing. Found duplicate or non-increasing time at index "
                << i << " (t[" << (i-1) << "] = " << table_[i-1].first()
                << ", t[" << i << "] = " << table_[i].first() << ")."
                << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::bodyLoad Foam::externalForceModels::tabulatedForce::evaluate
(
    const scalar t,
    const quaternion& /* q */
) const
{
    if (!active_ || table_.empty())
    {
        return bodyLoad();
    }

    vector F = vector::zero;

    if (t < table_.first().first())
    {
        if (outOfRange_ == "hold" || outOfRange_ == "clamp")
        {
            F = table_.first().second();
        }
        else // "zero"
        {
            return bodyLoad();
        }
    }
    else if (t > table_.last().first())
    {
        if (outOfRange_ == "hold" || outOfRange_ == "clamp")
        {
            F = table_.last().second();
        }
        else // "zero"
        {
            return bodyLoad();
        }
    }
    else
    {
        // Inside table range
        for (label i = 0; i < table_.size() - 1; ++i)
        {
            if (t >= table_[i].first() && t <= table_[i+1].first())
            {
                const scalar t0 = table_[i].first();
                const scalar t1 = table_[i+1].first();
                const vector& F0 = table_[i].second();
                const vector& F1 = table_[i+1].second();

                const scalar frac = (t - t0) / max(t1 - t0, SMALL);
                F = F0 + frac * (F1 - F0);
                break;
            }
        }
    }

    const vector moment = applicationPoint_ ^ F;
    return bodyLoad(F, moment);
}

// ************************************************************************* //
