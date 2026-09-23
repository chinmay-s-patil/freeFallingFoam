/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    tFFPimpleFoam (true Free-Falling PimpleFoam)

Description
    Transient solver for incompressible, turbulent flow of Newtonian fluids
    in a body-fixed reference frame with dynamic, force-coupled free fall.

    Frame acceleration and velocity emerge naturally from the instantaneous
    balance between gravity and CFD aerodynamic forces (p_rgh pressure +
    viscous stresses on specified body patches):

        m * a_frame = m * g + F_aero
        V_frame^(n+1) = V_frame^n + a_frame^n * dt

    Includes:
      - Uniform registered state fields (frameVelocity, frameAcceleration,
        altitude, gravity, rhoAtmosphere, pressureForce, viscousForce,
        aeroForce, gravityForce, netForce, relativeVelocity).
      - Logging to postProcessing/fallingMotion/motion.dat.
      - Terminal velocity diagnostic detection.
      - Ground interrupt check (groundInterrupt).
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
        " in a body-fixed frame with dynamic fluid-body force coupling."
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

    // Read fallingFrameDict, create uniform registered state fields & hydrostatic fields
    #include "createFallingFrameFields.H"

    // Total pressure field p = p_rgh + gh for visualization & post-processing
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

        // Dynamic update of aerodynamic force, net acceleration, frame velocity, altitude, and BCs
        #include "updateFallingFrameMotion.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                mesh.controlledUpdate();

                if (mesh.changing())
                {
                    MRF.update();

                    {
                        const volScalarField gh_old("gh_old", gh);
                        gh = (fallingFrameForce & mesh.C());
                        ghf = (fallingFrameForce & mesh.Cf());
                        p_rgh += gh_old - gh;
                    }

                    if (correctPhi)
                    {
                        phi = mesh.Sf() & Uf();
                        #include "correctPhi.H"
                        fvc::makeRelative(phi, U);
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            #include "UEqn.H"

            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        // Reconstruct total pressure p = p_rgh + gh
        p = p_rgh + gh;

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
