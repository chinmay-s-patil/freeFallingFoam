/*---------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::externalForceModel

Implementation
    Foam::externalForceModel
\*---------------------------------------------------------------------------*/

#include "externalForceModel.H"

namespace Foam
{
    defineTypeNameAndDebug(externalForceModel, 0);
    defineRunTimeSelectionTable(externalForceModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::externalForceModel::externalForceModel
(
    const word& name,
    const dictionary& dict
)
:
    name_(name),
    dict_(dict),
    active_(dict.getOrDefault<bool>("active", true))
{}


// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::externalForceModel> Foam::externalForceModel::New
(
    const word& name,
    const dictionary& dict
)
{
    const word modelType(dict.get<word>("type"));

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Unknown externalForceModel type '" << modelType
            << "' for force instance '" << name << "'." << nl << nl
            << "Valid externalForceModel types are:" << nl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<externalForceModel>(ctorPtr(name, dict));
}

// ************************************************************************* //
