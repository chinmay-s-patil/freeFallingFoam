/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 AUTHOR,AFFILIATION
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fallingInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "uniformDimensionedFields.H"

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void fallingInletVelocityFvPatchVectorField::readParams
(
    const dictionary& dict
)
{
    if (dict.found("u0") || dict.found("acceleration"))
    {
        // Parameters supplied inline in the patch entry — take them directly.
        // NOTE: Inline parameters operate in prescribed-motion mode, bypassing
        // dynamic aero-force coupling. To use dynamic solver coupling, omit u0
        // and acceleration from 0/U and configure them in constant/fallingFrameDict.
        overrideParams_ = true;
        u0_            = dict.getOrDefault<vector>("u0",           vector::zero);
        acceleration_  = dict.getOrDefault<vector>("acceleration", vector::zero);

        Info<< "  [fallingInletVelocity] Patch '" << this->patch().name()
            << "': using inline parameters (prescribed-motion mode)." << nl;

        if (dict.found("maxFrameSpeed"))
        {
            maxFrameSpeed_ = dict.get<scalar>("maxFrameSpeed");
        }
        else if (dict.found("terminalVelocity"))
        {
            maxFrameSpeed_ = dict.get<scalar>("terminalVelocity");
        }
        else
        {
            maxFrameSpeed_ = GREAT;
        }
    }
    else
    {
        // Fall through to constant/fallingFrameDict (shared with the solver).
        overrideParams_ = false;

        const IOdictionary frameDict
        (
            IOobject
            (
                "fallingFrameDict",
                this->db().time().constant(),
                this->db().time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false   // do not register — temporary read
            )
        );

        u0_ = frameDict.getOrDefault<vector>("u0", vector::zero);

        if (frameDict.found("initialAcceleration"))
        {
            acceleration_ = frameDict.get<vector>("initialAcceleration");
        }
        else if (frameDict.found("acceleration"))
        {
            acceleration_ = frameDict.get<vector>("acceleration");
        }
        else if (frameDict.found("gravity"))
        {
            acceleration_ = frameDict.get<vector>("gravity");
        }
        else
        {
            acceleration_ = vector(0, -9.81, 0);
        }

        if (frameDict.found("maxFrameSpeed"))
        {
            maxFrameSpeed_ = frameDict.get<scalar>("maxFrameSpeed");
        }
        else if (frameDict.found("terminalVelocity"))
        {
            maxFrameSpeed_ = frameDict.get<scalar>("terminalVelocity");
        }
        else
        {
            maxFrameSpeed_ = GREAT;
        }
    }
}


vector fallingInletVelocityFvPatchVectorField::inletVelocity
(
    const scalar t
) const
{
    vector Uframe = u0_ + acceleration_ * t;

    const scalar spd = mag(Uframe);
    if (spd > maxFrameSpeed_ && spd > SMALL)
    {
        Uframe *= maxFrameSpeed_ / spd;
    }

    // Fluid moves opposite to the falling body in the body-fixed frame.
    return -Uframe;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

fallingInletVelocityFvPatchVectorField::
fallingInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(p, iF),
    u0_(vector::zero),
    acceleration_(vector::zero),
    maxFrameSpeed_(GREAT),
    overrideParams_(false)
{
    refValue()      = vector::zero;
    refGrad()       = vector::zero;
    valueFraction() = 1.0;          // pure Dirichlet until first updateCoeffs
}


fallingInletVelocityFvPatchVectorField::
fallingInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchVectorField(p, iF),
    u0_(vector::zero),
    acceleration_(vector::zero),
    maxFrameSpeed_(GREAT),
    overrideParams_(false)
{
    readParams(dict);

    // Initialise mixed-BC internals.
    refGrad()       = vector::zero;
    valueFraction() = 1.0;          // start fully Dirichlet
    refValue()      = inletVelocity(0.0);

    // Initialise patch field value: honour user-specified "value" if present in 0/U,
    // otherwise default to refValue()
    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        fvPatchVectorField::operator=(refValue());
    }
}


fallingInletVelocityFvPatchVectorField::
fallingInletVelocityFvPatchVectorField
(
    const fallingInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchVectorField(ptf, p, iF, mapper),
    u0_(ptf.u0_),
    acceleration_(ptf.acceleration_),
    maxFrameSpeed_(ptf.maxFrameSpeed_),
    overrideParams_(ptf.overrideParams_)
{}


fallingInletVelocityFvPatchVectorField::
fallingInletVelocityFvPatchVectorField
(
    const fallingInletVelocityFvPatchVectorField& ptf
)
:
    mixedFvPatchVectorField(ptf),
    u0_(ptf.u0_),
    acceleration_(ptf.acceleration_),
    maxFrameSpeed_(ptf.maxFrameSpeed_),
    overrideParams_(ptf.overrideParams_)
{}


fallingInletVelocityFvPatchVectorField::
fallingInletVelocityFvPatchVectorField
(
    const fallingInletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(ptf, iF),
    u0_(ptf.u0_),
    acceleration_(ptf.acceleration_),
    maxFrameSpeed_(ptf.maxFrameSpeed_),
    overrideParams_(ptf.overrideParams_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fallingInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    vector Uinlet;
    if (overrideParams_)
    {
        const scalar t = this->db().time().timeOutputValue();
        Uinlet = inletVelocity(t);
    }
    else if (this->db().foundObject<uniformDimensionedVectorField>("relativeVelocity"))
    {
        const uniformDimensionedVectorField& Vrel =
            this->db().lookupObject<uniformDimensionedVectorField>("relativeVelocity");
        Uinlet = Vrel.value();
    }
    else if (this->db().foundObject<uniformDimensionedVectorField>("frameVelocity"))
    {
        const uniformDimensionedVectorField& Vframe =
            this->db().lookupObject<uniformDimensionedVectorField>("frameVelocity");
        Uinlet = -Vframe.value();
    }
    else
    {
        const scalar t = this->db().time().timeOutputValue();
        Uinlet = inletVelocity(t);
    }


    // Determine inflow vs outflow geometrically from test flow direction
    vector testU = Uinlet;
    if (mag(testU) < SMALL)
    {
        // At t=0 when velocity is zero, test direction using -acceleration
        if (this->db().foundObject<uniformDimensionedVectorField>("frameAcceleration"))
        {
            const uniformDimensionedVectorField& Acc =
                this->db().lookupObject<uniformDimensionedVectorField>("frameAcceleration");
            testU = -Acc.value();
        }
        else
        {
            testU = -acceleration_;
        }
    }

    const vectorField& Sf = patch().Sf();

    refGrad() = vector::zero;

    const scalar magTestU = mag(testU);

    forAll(Sf, facei)
    {
        const scalar flux = (Sf[facei] & testU);
        const scalar tol = 1e-6 * mag(Sf[facei]) * magTestU;

        if (flux < -tol)
        {
            // Inflow face: prescribe the frame-relative inlet velocity.
            refValue()[facei]      = Uinlet;
            valueFraction()[facei] = 1.0;
        }
        else
        {
            // Outflow or parallel face: zeroGradient
            refValue()[facei]      = vector::zero;
            valueFraction()[facei] = 0.0;
        }
    }

    mixedFvPatchVectorField::updateCoeffs();
}


void fallingInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    mixedFvPatchVectorField::write(os);

    if (overrideParams_)
    {
        os.writeEntry("u0",           u0_);
        os.writeEntry("acceleration", acceleration_);
        if (maxFrameSpeed_ < GREAT)
        {
            os.writeEntry("maxFrameSpeed", maxFrameSpeed_);
        }
    }
    // When overrideParams_ is false the reader finds parameters in
    // fallingFrameDict, so we do not duplicate them here.
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    fallingInletVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
