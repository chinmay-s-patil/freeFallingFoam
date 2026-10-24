/*---------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::externalForceModels::constantForce

Description
    Constant external force vector in global coordinates.

    Dictionary syntax:
    \verbatim
    constantLoad
    {
        type              constantForce;
        active            yes;
        force             (0 0 -50);           // Force vector [N]
        applicationPoint  (0 0 0);             // Optional point [m]
    }
    \endverbatim
\*---------------------------------------------------------------------------*/

#ifndef constantForceModel_H
#define constantForceModel_H

#include "externalForceModel.H"
#include "point.H"

namespace Foam
{
namespace externalForceModels
{

/*---------------------------------------------------------------------------*\
                     Class constantForce Declaration
\*---------------------------------------------------------------------------*/

class constantForce
:
    public externalForceModel
{
private:

    // Private Data

        //- Constant force vector in global coordinates [N]
        vector force_;

        //- Application point relative to centre of mass in global coordinates [m]
        point applicationPoint_;


public:

    //- Runtime type information
    TypeName("constantForce");


    // Constructors

        //- Construct from name and dictionary
        constantForce(const word& name, const dictionary& dict);


    //- Destructor
    virtual ~constantForce() = default;


    // Member Functions

        //- Evaluate force and moment in global coordinates
        virtual bodyLoad evaluate
        (
            const scalar t,
            const quaternion& q = quaternion::I
        ) const override;
};

} // End namespace externalForceModels
} // End namespace Foam

#endif
