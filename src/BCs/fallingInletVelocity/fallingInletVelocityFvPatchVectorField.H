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

Class
    Foam::fallingInletVelocityFvPatchVectorField

Group
    grpInletBoundaryConditions

Description
    Inlet/outlet velocity boundary condition for a body-fixed (falling/
    accelerating) reference frame.

    On faces where flow enters the domain (phi < 0) the patch velocity is
    set to the frame-relative inlet value:

        U_inlet(t) = -(u0 + acceleration * t)   [clamped to maxFrameSpeed]

    On faces where flow exits (phi >= 0) a zero-gradient (Neumann) condition
    is applied, allowing wake fluctuations to leave without reflecting.

    This mixed Dirichlet/Neumann behaviour prevents blow-up when unsteady
    structures reach the far-field boundaries.

    Parameters are read either from the patch dictionary (inline) or from
    constant/fallingFrameDict (shared with the solver).

Usage
    Inline in 0/U patch entry (overrides fallingFrameDict):
    \verbatim
    <patchName>
    {
        type            fallingInletVelocity;
        u0              (0 0 0);
        acceleration    (0 -9.81 0);
        maxFrameSpeed   85;          // optional, default: no clamp
    }
    \endverbatim

    Or shared-dict mode (reads constant/fallingFrameDict automatically):
    \verbatim
    <patchName>
    {
        type            fallingInletVelocity;
    }
    \endverbatim

SourceFiles
    fallingInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef fallingInletVelocityFvPatchVectorField_H
#define fallingInletVelocityFvPatchVectorField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class fallingInletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class fallingInletVelocityFvPatchVectorField
:
    public mixedFvPatchVectorField          // <-- was fixedValue
{
    // Private data

        //- Initial frame velocity [m/s]
        vector u0_;

        //- Constant frame acceleration [m/s²]
        vector acceleration_;

        //- Maximum frame speed clamp [m/s]
        scalar maxFrameSpeed_;

        //- True when parameters were read from the patch dict itself
        bool overrideParams_;


    // Private member functions

        //- Read u0, acceleration, maxFrameSpeed from dict or fallingFrameDict
        void readParams(const dictionary& dict);

        //- Return U_inlet = -(u0 + acceleration*t), magnitude-clamped
        vector inletVelocity(const scalar t) const;


public:

    //- Runtime type information
    TypeName("fallingInletVelocity");


    // Constructors

        //- Construct from patch and internal field
        fallingInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fallingInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patch field onto a new patch
        fallingInletVelocityFvPatchVectorField
        (
            const fallingInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy construct
        fallingInletVelocityFvPatchVectorField
        (
            const fallingInletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new fallingInletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fallingInletVelocityFvPatchVectorField
        (
            const fallingInletVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new fallingInletVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
