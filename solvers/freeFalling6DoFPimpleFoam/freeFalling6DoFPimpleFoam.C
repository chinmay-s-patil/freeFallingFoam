/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    freeFalling6DoFPimpleFoam

Group
    grpIncompressibleSolvers

Description
    Transient solver for incompressible, turbulent flow of Newtonian fluids
    on a moving mesh, with a falling/accelerating reference frame.

    Uses the p_rgh (hydrostatic-reduced pressure) formulation so that
    the 6DoF rigid-body motion solver sees only aerodynamic forces,
    not the spurious hydrostatic/buoyancy contribution from the
    falling-frame body force.

    p_rgh = p - (fallingFrameForce & position)

    This prevents the angular velocity blowup that occurs when the
    6DoF solver integrates the total pressure (which includes the
    hydrostatic gradient) over the body surface.

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p_rgh   | Hydrostatic-reduced kinematic pressure [m2/s2]
        \<turbulence fields\> | As required by user selection
    \endplaintable

    \heading Required dictionaries
    \plaintable
        constant/fallingFrameDict | u0, acceleration, maxFrameSpeed
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, turbulent flow"
        " with a falling reference frame and p_rgh formulation"
        " for stable 6DoF rigid-body coupling."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"          // reads p_rgh, U, phi
    #include "createUfIfPresent.H"

    // ---------------------------------------------------------------
    //  Read falling-frame parameters, create gh, ghf, magUStar
    // ---------------------------------------------------------------
    #include "createFallingFrameFields.H"

    // ---------------------------------------------------------------
    //  Create total pressure p = p_rgh + gh (for post-processing
    //  and for the 6DoF force computation via function objects)
    // ---------------------------------------------------------------
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + gh
    );

    turbulence->validate();

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                // Do any mesh changes
                mesh.controlledUpdate();

                if (mesh.changing())
                {
                    MRF.update();

                    // Update hydrostatic head after mesh motion
                    // (cell/face centres have moved).
                    // Adjust p_rgh so that total pressure p = p_rgh + gh
                    // remains continuous across the mesh motion step.
                    // Without this, the change in gh (from cell centres
                    // moving with the 6DoF rotation) causes a pressure
                    // discontinuity that generates spurious forces/torques
                    // on the body, growing with rotation angle until
                    // divergence.  (Same pattern as buoyantPimpleFoam.)
                    {
                        const volScalarField gh_old("gh_old", gh);
                        gh = (fallingFrameForce & mesh.C());
                        ghf = (fallingFrameForce & mesh.Cf());
                        p_rgh += gh_old - gh;
                    }

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            #include "UEqn.H"   // p_rgh formulation (no body force source)

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"   // solves for p_rgh, reconstructs p
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
