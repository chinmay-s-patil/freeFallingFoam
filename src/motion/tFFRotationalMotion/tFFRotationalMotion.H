/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::solidBodyMotionFunctions::tFFRotationalMotion

Description
    Solid-body rotational motion driven by 3D attitude quaternions registered
    in the mesh registry by tFF6DoFPimpleFoam (attitudeW, attitudeVec).

    Accommodates arbitrary 3D multi-axis tumbling across a spherical AMI
    interface or overset grid without mesh distortion.

SourceFiles
    tFFRotationalMotion.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_tFFRotationalMotion_H
#define Foam_tFFRotationalMotion_H

#include "solidBodyMotionFunction.H"
#include "primitiveFields.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                    Class tFFRotationalMotion Declaration
\*---------------------------------------------------------------------------*/

class tFFRotationalMotion
:
    public solidBodyMotionFunction
{
    // Private Data

        //- Center of rotation / hub position
        point origin_;


public:

    //- Runtime type information
    TypeName("tFFRotationalMotion");


    // Constructors

        //- Construct from components
        tFFRotationalMotion
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );

        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const override
        {
            return autoPtr<solidBodyMotionFunction>
            (
                new tFFRotationalMotion
                (
                    SBMFCoeffs_,
                    time_
                )
            );
        }

        //- Disallow default bitwise copy construct
        tFFRotationalMotion(const tFFRotationalMotion&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const tFFRotationalMotion&) = delete;


    //- Destructor
    virtual ~tFFRotationalMotion() = default;


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const override;

        //- Update properties from a dictionary
        virtual bool read(const dictionary& SBMFCoeffs) override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
