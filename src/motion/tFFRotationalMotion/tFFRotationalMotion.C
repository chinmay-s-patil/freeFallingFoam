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
    Implementation of the registry-driven 3D attitude motion function.

\*---------------------------------------------------------------------------*/

#include "tFFRotationalMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"
#include "uniformDimensionedFields.H"
#include "quaternion.H"
#include "septernion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(tFFRotationalMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        tFFRotationalMotion,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::tFFRotationalMotion::tFFRotationalMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime),
    origin_(SBMFCoeffs.getOrDefault<point>("origin", point::zero))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::tFFRotationalMotion::transformation() const
{
    if (!time_.foundObject<polyMesh>(polyMesh::defaultRegion))
    {
        return septernion::I;
    }

    const polyMesh& mesh = time_.lookupObject<polyMesh>(polyMesh::defaultRegion);

    if
    (
        !mesh.foundObject<uniformDimensionedScalarField>("attitudeW")
     || !mesh.foundObject<uniformDimensionedVectorField>("attitudeVec")
    )
    {
        return septernion::I;
    }

    const scalar w =
        mesh.lookupObject<uniformDimensionedScalarField>("attitudeW").value();
    const vector v =
        mesh.lookupObject<uniformDimensionedVectorField>("attitudeVec").value();

    quaternion q(w, v);
    const scalar magQ = mag(q);
    if (magQ > 1e-12)
    {
        q.normalise();
    }
    else
    {
        q = quaternion::I;
    }

    return septernion(-origin_) * q * septernion(origin_);
}


bool Foam::solidBodyMotionFunctions::tFFRotationalMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    if (solidBodyMotionFunction::read(SBMFCoeffs))
    {
        SBMFCoeffs.readIfPresent("origin", origin_);
        return true;
    }

    return false;
}

// ************************************************************************* //
