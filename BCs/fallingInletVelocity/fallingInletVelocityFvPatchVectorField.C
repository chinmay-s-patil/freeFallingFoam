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
        overrideParams_ = true;
        u0_            = dict.getOrDefault<vector>("u0",           vector::zero);
        acceleration_  = dict.getOrDefault<vector>("acceleration", vector::zero);
        maxFrameSpeed_ = dict.getOrDefault<scalar>("maxFrameSpeed", GREAT);
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

        u0_           = frameDict.get<vector>("u0");
        acceleration_ = frameDict.get<vector>("acceleration");

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

    // Do NOT call the mixedFvPatchVectorField dict constructor — we handle
    // the "value" entry ourselves to avoid requiring it in 0/U.
    fvPatchVectorField::operator=(refValue());
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

    const scalar t = this->db().time().timeOutputValue();
    const vector Uinlet = inletVelocity(t);

    const scalarField& phip =
        patch().lookupPatchField<surfaceScalarField>("phi");

    // Pull the adjacent cell values so outflow faces can mirror them.
    const vectorField Uinternal(patchInternalField());

    // refGrad is never used: we always set valueFraction = 1 (pure Dirichlet).
    refGrad() = vector::zero;

    forAll(phip, facei)
    {
        if (phip[facei] < 0)
        {
            // Inflow face: prescribe the frame-relative inlet velocity.
            refValue()[facei]      = Uinlet;
            valueFraction()[facei] = 1.0;
        }
        else
        {
            // Outflow face: mirror the interior cell value.
            // valueFraction=1 with refValue=internal is mathematically
            // identical to zeroGradient but assembles as a Dirichlet
            // coefficient, which keeps the face flux consistent with
            // the divergence-free constraint and avoids continuity errors.
            refValue()[facei]      = Uinternal[facei];
            valueFraction()[facei] = 1.0;
        }
    }

    mixedFvPatchVectorField::updateCoeffs();
}


void fallingInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);

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

    fvPatchVectorField::writeValueEntry(os);
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
