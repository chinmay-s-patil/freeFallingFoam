/*---------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::externalForceModels::sixDoFThrustForce

Description
    Body-fixed frame thrust rotating dynamically with body attitude.
    Thrust direction and application point are defined in the body-fixed
    coordinate frame, and transformed dynamically to global coordinates:
        F_global = q * F_body * q^-1
        r_global = q * r_body * q^-1
        M_global = r_global ^ F_global

    Dictionary syntax:
    \verbatim
    rocketThrust
    {
        type              6DoFThrust;
        active            yes;
        direction         (1 0 0);             // Direction in body coordinates
        magnitude         1200;                // Thrust magnitude [N]
        startTime         0.0;                 // Ignition time [s] (default: 0.0)
        endTime           4.0;                 // Shutdown time [s] (default: GREAT)
        applicationPoint  (-0.45 0 0);         // Point in body coordinates [m] (default: zero)
    }
    \endverbatim
\*---------------------------------------------------------------------------*/

#ifndef sixDoFThrustForceModel_H
#define sixDoFThrustForceModel_H

#include "externalForceModel.H"
#include "point.H"

namespace Foam
{
namespace externalForceModels
{

/*---------------------------------------------------------------------------*\
                   Class sixDoFThrustForce Declaration
\*---------------------------------------------------------------------------*/

class sixDoFThrustForce
:
    public externalForceModel
{
private:

    // Private Data

        //- Unit direction vector in body-fixed coordinates
        vector direction_;

        //- Thrust magnitude [N]
        scalar magnitude_;

        //- Ignition start time [s]
        scalar startTime_;

        //- Shutdown end time [s]
        scalar endTime_;

        //- Application point in body-fixed coordinates [m]
        point applicationPoint_;


public:

    //- Runtime type information
    TypeName("6DoFThrust");


    // Constructors

        //- Construct from name and dictionary
        sixDoFThrustForce(const word& name, const dictionary& dict);


    //- Destructor
    virtual ~sixDoFThrustForce() = default;


    // Member Functions

        //- Evaluate thrust force and moment transformed to global coordinates
        virtual bodyLoad evaluate
        (
            const scalar t,
            const quaternion& q = quaternion::I
        ) const override;
};

} // End namespace externalForceModels
} // End namespace Foam

#endif
