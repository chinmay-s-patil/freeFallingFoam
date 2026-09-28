/*---------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::meanFlowModel

Implementation
    Foam::meanFlowModel

\*---------------------------------------------------------------------------*/

#include "meanFlowModel.H"

namespace Foam
{
    defineTypeNameAndDebug(meanFlowModel, 0);
    defineRunTimeSelectionTable(meanFlowModel, dictionary);
}

Foam::meanFlowModel::meanFlowModel(const dictionary& dict)
:
    dict_(dict)
{}

Foam::autoPtr<Foam::meanFlowModel> Foam::meanFlowModel::New(const dictionary& dict)
{
    const dictionary& meanFlowDict = dict.found("meanFlow") ? dict.subDict("meanFlow") : dict;
    const word modelType = meanFlowDict.getOrDefault<word>("model", "constant");

    Info<< "Selecting meanFlowModel " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Unknown meanFlowModel type " << modelType << nl << nl
            << "Valid meanFlowModel types are:" << nl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<meanFlowModel>(ctorPtr(meanFlowDict));
}
