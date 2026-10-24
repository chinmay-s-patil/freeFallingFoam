/*---------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::externalForceModels::tabulatedForce

Description
    Time-dependent force vector linearly interpolated from a user-prescribed
    time-force table.

    Dictionary syntax:
    \verbatim
    prescribedLoad
    {
        type              tabulatedForce;
        active            yes;
        interpolation     linear;              // Interpolation scheme (default: linear)
        outOfRange        zero;                // Behaviour outside range: zero, hold (default: zero)
        applicationPoint  (0 0 0);             // Optional application point [m]
        table
        (
            (0.0  (0 0 0))
            (0.5  (0 200 0))
            (1.0  (0 0 0))
            (2.0  (0 -50 0))
        );
    }
    \endverbatim
\*---------------------------------------------------------------------------*/

#ifndef tabulatedForceModel_H
#define tabulatedForceModel_H

#include "externalForceModel.H"
#include "List.H"
#include "Tuple2.H"
#include "point.H"

namespace Foam
{
namespace externalForceModels
{

/*---------------------------------------------------------------------------*\
                     Class tabulatedForce Declaration
\*---------------------------------------------------------------------------*/

class tabulatedForce
:
    public externalForceModel
{
private:

    // Private Data

        //- Table of (time, forceVector) pairs
        List<Tuple2<scalar, vector>> table_;

        //- Out of range handling ("zero" or "hold")
        word outOfRange_;

        //- Interpolation scheme ("linear")
        word interpolation_;

        //- Application point in global coordinates [m]
        point applicationPoint_;


public:

    //- Runtime type information
    TypeName("tabulatedForce");


    // Constructors

        //- Construct from name and dictionary
        tabulatedForce(const word& name, const dictionary& dict);


    //- Destructor
    virtual ~tabulatedForce() = default;


    // Member Functions

        //- Evaluate force and moment in global coordinates at time t
        virtual bodyLoad evaluate
        (
            const scalar t,
            const quaternion& q = quaternion::I
        ) const override;
};

} // End namespace externalForceModels
} // End namespace Foam

#endif
