/*---------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::functionObjects::fallingMotion

Description
    Dedicated function object for tracking, displaying, and logging
    dynamic free-fall motion fields (frameVelocity, frameAcceleration, altitude,
    forces, and atmospheric density).

    Allows fine-grained toggling of console stdout output and CSV logging to
    postProcessing/fallingMotion/motion.dat.

\*---------------------------------------------------------------------------*/

#ifndef fallingMotion_H
#define fallingMotion_H

#include "functionObject.H"
#include "IOdictionary.H"
#include "OFstream.H"
#include "fvMesh.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class fallingMotion Declaration
\*---------------------------------------------------------------------------*/

class fallingMotion
:
    public functionObject
{
    // Private Data

        //- Time database reference
        const Time& time_;

        //- File output flag
        bool writeCSV_;

        //- Output file stream pointer
        autoPtr<OFstream> motionFilePtr_;

        //- Console print flags
        bool printFrameVel_;
        bool printFrameAcc_;
        bool printAltitude_;
        bool printPressureForce_;
        bool printViscousForce_;
        bool printAeroForce_;
        bool printGravityForce_;
        bool printNetForce_;
        bool printMeanFlow_;
        bool printRelativeVel_;


    // Private Member Functions

        //- Open log file and write header
        void openFile();

        //- Output console information
        void printConsole(const fvMesh& mesh);

        //- Write data to motion.dat
        void writeLog(const fvMesh& mesh);

public:

    //- Runtime type information
    TypeName("fallingMotion");

    // Constructors

        //- Construct from time and dictionary
        fallingMotion
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

    //- Destructor
    virtual ~fallingMotion() = default;

    // Member Functions

        //- Read configuration dictionary
        virtual bool read(const dictionary& dict) override;

        //- Execute at runtime time steps
        virtual bool execute() override;

        //- Write data at write time steps
        virtual bool write() override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
