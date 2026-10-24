/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::solidBodyMotionFunctions::sixDoFRotatingMotion

Description
    6DoF axial rotation model driven dynamically by fluid forces and moments
    integrated over specified body patches (e.g. propeller wall).

    Solves the angular momentum equation about the rotation axis:
        I * d(omega)/dt = Torque_aero - damping * omega
        d(theta)/dt = omega

    Produces continuous, unbounded solid-body rotation of a cellZone across an
    AMI (Arbitrary Mesh Interface) without mesh deformation or distortion.

SourceFiles
    sixDoFRotatingMotion.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_sixDoFRotatingMotion_H
#define Foam_sixDoFRotatingMotion_H

#include "solidBodyMotionFunction.H"
#include "primitiveFields.H"
#include "point.H"
#include "wordRes.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                   Class sixDoFRotatingMotion Declaration
\*---------------------------------------------------------------------------*/

class sixDoFRotatingMotion
:
    public solidBodyMotionFunction
{
    // Private Data

        //- Center of rotation / hub position
        point origin_;

        //- Rotation axis unit vector
        vector axis_;

        //- Patches on which fluid forces/moments are evaluated
        wordRes patches_;

        //- Mass of the body [kg]
        scalar mass_;

        //- Moment of inertia about the spin axis [kg*m^2]
        scalar momentOfInertia_;

        //- Initial angular velocity [rad/s]
        scalar omega0_;

        //- Rotational damping coefficient [N*m*s/rad]
        scalar damping_;

        //- Reference density for kinematic pressure [kg/m^3]
        scalar rhoInf_;

        //- Density field name ("rhoInf" or "rho")
        word rhoName_;

        //- Report progress to stdout
        bool report_;

        //- Current angular velocity [rad/s]
        mutable scalar omega_;

        //- Current angular position [rad]
        mutable scalar theta_;

        //- Current aerodynamic torque [N*m]
        mutable scalar torque_;

        //- Last updated time index
        mutable label lastTimeIndex_;

        //- Optional motion tracking log file
        mutable autoPtr<OFstream> logFilePtr_;


    // Private Member Functions

        //- Read initial/restarted motion state from time directory
        void readState();

        //- Log motion dynamics to postProcessing/sixDoFRotatingMotion/motion.dat
        void logMotion
        (
            const scalar t,
            const scalar torque,
            const scalar omega,
            const scalar rpm,
            const scalar deg
        ) const;

        //- No copy construct
        sixDoFRotatingMotion(const sixDoFRotatingMotion&) = delete;

        //- No copy assignment
        void operator=(const sixDoFRotatingMotion&) = delete;


public:

    //- Runtime type information
    TypeName("sixDoFRotatingMotion");


    // Constructors

        //- Construct from components
        sixDoFRotatingMotion
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );

        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const
        {
            return autoPtr<solidBodyMotionFunction>
            (
                new sixDoFRotatingMotion
                (
                    SBMFCoeffs_,
                    time_
                )
            );
        }


    //- Destructor
    virtual ~sixDoFRotatingMotion() = default;


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs);

        //- Write state in dictionary format
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
