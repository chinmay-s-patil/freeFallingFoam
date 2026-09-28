/*---------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::meanFlowModel

Description
    Abstract base class for mean atmospheric flow velocity models U_mean(z).

\*---------------------------------------------------------------------------*/

#ifndef meanFlowModel_H
#define meanFlowModel_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "dictionary.H"
#include "vector.H"
#include "autoPtr.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class meanFlowModel Declaration
\*---------------------------------------------------------------------------*/

class meanFlowModel
{
protected:

    // Private Data
        const dictionary& dict_;

public:

    //- Runtime type information
    TypeName("meanFlowModel");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        meanFlowModel,
        dictionary,
        (const dictionary& dict),
        (dict)
    );

    // Constructors
        meanFlowModel(const dictionary& dict);

    // Selector
        static autoPtr<meanFlowModel> New(const dictionary& dict);

    //- Destructor
    virtual ~meanFlowModel() = default;

    // Member Functions
        virtual vector velocity(const scalar altitude) const = 0;
};

} // End namespace Foam

#endif
