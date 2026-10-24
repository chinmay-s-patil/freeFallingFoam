/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    tFF6DoFPimpleFoam (true Free-Falling 6DoF PimpleFoam)

Description
    Transient solver for incompressible, turbulent flow in a dynamically
    accelerating falling reference frame with full 6DoF fluid-body coupling
    (3 translational DOFs + 3 rotational DOFs).

    Translational dynamics:
        m * a_frame = m * g + F_aero
        V_frame^(n+1) = V_frame^n + a_frame * dt

    Rotational dynamics (Euler's equations solved in body-fixed principal axes):
        I_xx * d(omega_x)/dt + (I_zz - I_yy) * omega_y * omega_z = M_x
        I_yy * d(omega_y)/dt + (I_xx - I_zz) * omega_z * omega_x = M_y
        I_zz * d(omega_z)/dt + (I_yy - I_xx) * omega_x * omega_y = M_z

    Orientation tracking via unit quaternions:
        d(q)/dt = 0.5 * q * omega_body

    Mesh motion:
        The falling frame remains non-rotating (no fictitious Coriolis or
        centrifugal forces in the fluid momentum equations).
        The inner spherical domain rotates physically across a spherical AMI
        interface via tFFRotationalMotion, driven by attitudeW and attitudeVec.

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
#include "meanFlowModel.H"
#include "externalForces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, turbulent flow"
        " in a falling reference frame with full 6DoF translational"
        " and rotational fluid-body coupling."
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

    // Read fallingFrameDict, create uniform registered translation fields
    #include "createFallingFrameFields.H"

    // Read rotor settings, create uniform registered attitude fields
    #include "createRotorState.H"

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

        // Dynamic update: predictor step for frame velocity, altitude, and BCs
        #include "predictFallingFrameMotion.H"

        // Dynamic update: predictor step for rotor attitude and registered fields
        #include "predictRotorDynamics.H"

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

        // Dynamic update: corrector step for frame velocity, altitude, and forces
        #include "correctFallingFrameMotion.H"

        // Dynamic update: corrector step for rotor moments, angular velocity, and attitude
        #include "correctRotorDynamics.H"

        // Reconstruct total pressure p = p_rgh + gh
        p = p_rgh + gh;

        if (triggerTerminalVelocityWrite)
        {
            Info<< "[DIAGNOSTIC] Saving fields at terminal velocity state (t = "
                << runTime.timeName() << " s)" << nl << endl;
            runTime.writeNow();
            triggerTerminalVelocityWrite = false;

            if (terminalVelocityStop)
            {
                Info<< "[DIAGNOSTIC] Terminal velocity stop condition reached. Ending simulation." << nl << endl;
                break;
            }
        }
        else if (triggerGroundInterrupt)
        {
            Info<< "[GROUND INTERRUPT] Body reached ground altitude. Saving state and ending simulation." << nl << endl;
            runTime.writeNow();
            break;
        }
        else if (triggerMachStopper)
        {
            Info<< "[MACH REGIME STOPPER] Mach cutoff reached. Saving state and ending simulation." << nl << endl;
            runTime.writeNow();
            break;
        }
        else
        {
            runTime.write();
        }

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
