/*---------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::externalForces

Description
    Manager class for configurable collections of external body forces
    specified in constant/fallingFrameDict under the "externalForces" subdictionary.
\*---------------------------------------------------------------------------*/

#ifndef externalForces_H
#define externalForces_H

#include "typeInfo.H"
#include "dictionary.H"
#include "PtrList.H"
#include "fvMesh.H"
#include "externalForceModel.H"
#include "bodyLoad.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class externalForces Declaration
\*---------------------------------------------------------------------------*/

class externalForces
{
private:

    // Private Data

        //- Master enable switch
        bool enabled_;

        //- List of active external force models
        PtrList<externalForceModel> forces_;


public:

    //- Runtime type information
    TypeName("externalForces");


    // Constructors

        //- Construct from fallingFrameDict or parent dictionary
        externalForces(const dictionary& dict);


    //- Destructor
    virtual ~externalForces() = default;


    // Member Functions

        //- Check if external forces are enabled
        bool enabled() const
        {
            return enabled_;
        }

        //- Number of configured force models
        label size() const
        {
            return forces_.size();
        }

        //- Access to configured force models
        const PtrList<externalForceModel>& forces() const
        {
            return forces_;
        }

        //- Evaluate total load at time t with explicit orientation quaternion
        bodyLoad evaluate
        (
            const scalar t,
            const quaternion& q = quaternion::I
        ) const;

        //- Evaluate total load at time t, querying attitude from mesh if 6DoF
        bodyLoad evaluate
        (
            const scalar t,
            const fvMesh& mesh
        ) const;

        //- Print detailed load summary to Info
        void printLoads
        (
            const scalar t,
            const fvMesh& mesh
        ) const;
};

} // End namespace Foam

#endif
