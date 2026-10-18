// ============================================================
//  correctFallingFrameMotion.H — Predictor-Corrector Corrector Step
//
//  Trapezoidal (Crank-Nicolson) corrector using converged aero
//  forces from the PIMPLE / density solve at this timestep.
//  See predictFallingFrameMotion.H for full scheme description.
// ============================================================

#include "OFstream.H"

// ------------------------------------------------------------
//  1. Aerodynamic force integration from converged flow fields
// ------------------------------------------------------------
vector fPress = vector::zero;
vector fVisc = vector::zero;

#if defined(FF_INCOMPRESSIBLE)
tmp<volSymmTensorField> tdevTau = turbulence->devReff();
const volSymmTensorField& devTau = tdevTau();

forAll(patchNames, i)
{
    const label patchi = mesh.boundaryMesh().findPatchID(patchNames[i]);
    if (patchi != -1)
    {
        const fvPatchScalarField& p_rgh_p = p_rgh.boundaryField()[patchi];
        const fvPatchSymmTensorField& devTau_p = devTau.boundaryField()[patchi];
        const vectorField& Sf = mesh.boundary()[patchi].Sf();

        forAll(p_rgh_p, facei)
        {
            // Pressure force: p_rgh is kinematic [m2/s2] in incompressible solver,
            // multiplied by rhoAtmosphere [kg/m3] to obtain force [N].
            fPress += p_rgh_p[facei] * Sf[facei] * rhoAtmosphere.value();

            // Viscous shear force exerted on body wall by fluid: (devTau & Sf) * rho
            fVisc += (devTau_p[facei] & Sf[facei]) * rhoAtmosphere.value();
        }
    }
}
#elif defined(FF_COMPRESSIBLE)
tmp<volSymmTensorField> tdevTau = turbulence->devRhoReff();
const volSymmTensorField& devTau = tdevTau();

forAll(patchNames, i)
{
    const label patchi = mesh.boundaryMesh().findPatchID(patchNames[i]);
    if (patchi != -1)
    {
        const fvPatchScalarField& p_p = p.boundaryField()[patchi];
        const fvPatchSymmTensorField& devTau_p = devTau.boundaryField()[patchi];
        const vectorField& Sf = mesh.boundary()[patchi].Sf();

        forAll(p_p, facei)
        {
            // Pressure force: p is static pressure [Pa] in compressible solver
            fPress += p_p[facei] * Sf[facei];

            // Viscous shear force exerted on body wall by fluid: (devTau & Sf)
            fVisc += (devTau_p[facei] & Sf[facei]);
        }
    }
}
#endif

reduce(fPress, sumOp<vector>());
reduce(fVisc, sumOp<vector>());

pressureForce.value() = fPress;
viscousForce.value()  = fVisc;
aeroForce.value()     = fPress + fVisc;

// ------------------------------------------------------------
//  2. Force balance and new acceleration at timestep n+1
// ------------------------------------------------------------
gravityForce.value() = bodyMass * gravity.value();
netForce.value()     = gravityForce.value() + aeroForce.value();

vector a_new = netForce.value() / bodyMass;
a_new.x() *= freedomMask.x();
a_new.y() *= freedomMask.y();
a_new.z() *= freedomMask.z();

// ------------------------------------------------------------
//  3. Heun Corrector: V^{n+1} = V^n + 0.5 * (a^n + a^{n+1}) * dt
// ------------------------------------------------------------
const vector v_corrected = v_prev + 0.5 * (a_prev + a_new) * dt;

const scalar V_up_new = -(v_corrected & gDir);
const scalar alt_corrected = alt_prev + 0.5 * (V_up_prev + V_up_new) * dt;

frameAcceleration.value() = a_new;
frameVelocity.value()     = v_corrected;
altitude.value()          = alt_corrected;

// Update density and gravity at corrected altitude
if (densityModel == "altitudeBased")
{
    const scalar effectiveAlt = max(alt_corrected - groundAltitude, 0.0);
    rhoAtmosphere.value() = rhoRef * std::exp(-effectiveAlt / scaleHeight);
}
else
{
    rhoAtmosphere.value() = rhoRef;
}

if (gravityModel == "altitudeBased")
{
    const scalar magG0 = max(mag(initialAcc), SMALL);
    const scalar effectiveAlt = max(alt_corrected - groundAltitude, 0.0);
    const scalar gMag = magG0 * sqr(earthRadius / (earthRadius + effectiveAlt));
    gravity.value() = (initialAcc / magG0) * gMag;
}
else
{
    gravity.value() = initialAcc;
}

// Update mean atmospheric flow and relative freestream velocity
if (meanFlowPtr.valid())
{
    meanFlowVelocity.value() = meanFlowPtr->velocity(altitude.value());
}
relativeVelocity.value() = meanFlowVelocity.value() - frameVelocity.value();

// Update pseudo-force
fallingFrameForce.value() = -frameAcceleration.value();

#if defined(FF_INCOMPRESSIBLE)
gh = (fallingFrameForce & mesh.C());
ghf = (fallingFrameForce & mesh.Cf());
#endif

// Correct boundary conditions to final converged frame velocity
U.correctBoundaryConditions();

// Normalized velocity magnitude for post-processing: |U| / max(|V_frame|, 1.0)
{
    const scalar Uref = max(mag(frameVelocity.value()), 1.0);
    magUStar = mag(U) / dimensionedScalar("Uref", dimVelocity, Uref);
}

// ------------------------------------------------------------
//  4. Speed of sound and Mach number calculations
// ------------------------------------------------------------
#if defined(FF_INCOMPRESSIBLE)
{
    const scalar rhoVal = max(rhoAtmosphere.value(), SMALL);
    const scalar pAtm = pAtmSeaLevel * (rhoVal / max(rhoRef, SMALL));
    const scalar a_sound = Foam::sqrt(gammaCoeff * (pAtm / rhoVal));
    speedOfSound.value() = a_sound;

    const scalar vRelMag = mag(relativeVelocity.value());
    machNumber.value() = vRelMag / max(a_sound, SMALL);

    if (machNumber.value() > 0.3)
    {
        Info<< "  [DIAGNOSTIC] Compressibility warning: Mach number ("
            << machNumber.value() << ") > 0.3 in incompressible solver. "
            << "Consider using tFFRhoPimpleFoam or tFFRhoCentralFoam." << nl;
    }
}
#elif defined(FF_COMPRESSIBLE)
{
    // Compute speed of sound directly from thermodynamic state
    const volScalarField c
    (
        "c",
        Foam::sqrt(thermo.Cp() / thermo.Cv() * (1.0 / psi))
    );
    speedOfSound.value() = gAverage(c);

    const scalar vRelMag = mag(relativeVelocity.value());
    machNumber.value() = vRelMag / max(speedOfSound.value(), SMALL);

    #if defined(FF_PIMPLE_COMPRESSIBLE)
    volScalarField localMach
    (
        "localMach",
        mag(U) / max(c, dimensionedScalar("cMin", dimVelocity, SMALL))
    );
    maxLocalMachNumber.value() = gMax(localMach);

    // Dynamic Mach regime switch or stopper
    if (machRegime == "auto" || machRegime == "autoSwitch")
    {
        if (maxLocalMachNumber.value() >= machRegimeSwitchCutOff)
        {
            if (machRegimeAction == "switch" && !machRegimeSwitched)
            {
                machRegimeSwitched = true;
                if (mesh.foundObject<IOdictionary>("fvSolution"))
                {
                    IOdictionary& solDict = mesh.lookupObjectRef<IOdictionary>("fvSolution");
                    if (solDict.found("PIMPLE"))
                    {
                        dictionary& pimpleDict = solDict.subDict("PIMPLE");
                        pimpleDict.set("transonic", true);
                    }
                }
                Info<< "\n[MACH REGIME SWITCH] Max local Mach number (" << maxLocalMachNumber.value()
                    << ") reached cutoff (" << machRegimeSwitchCutOff
                    << "). Dynamically enabling PIMPLE transonic mode in memory." << nl << endl;
            }
            else if (machRegimeAction == "stop")
            {
                Info<< "\n[MACH REGIME STOPPER] Max local Mach number (" << maxLocalMachNumber.value()
                    << ") reached cutoff (" << machRegimeSwitchCutOff
                    << "). Stopper triggered. Writing output and ending simulation." << nl << endl;
                runTime.writeNow();
                runTime.setTime(runTime.endTime().value(), runTime.timeIndex());
            }
        }
    }
    #endif
}
#endif

// ------------------------------------------------------------
//  5. Terminal velocity detection (O(1) circular ring buffer)
// ------------------------------------------------------------
{
    velHistory[velHistoryIndex] = mag(frameVelocity.value());
    velHistoryIndex = (velHistoryIndex + 1) % terminalVelocitySamples;
    if (velHistoryCount < terminalVelocitySamples)
    {
        velHistoryCount++;
    }

    if (velHistoryCount >= terminalVelocitySamples)
    {
        scalar minV = velHistory[0];
        scalar maxV = velHistory[0];
        for (label vi = 1; vi < terminalVelocitySamples; ++vi)
        {
            minV = min(minV, velHistory[vi]);
            maxV = max(maxV, velHistory[vi]);
        }
        const scalar vVar = maxV - minV;

        if (mag(frameAcceleration.value()) < terminalAccelerationTolerance && vVar < terminalVelocityTolerance)
        {
            if (!terminalVelocityTriggered)
            {
                terminalVelocityTriggered = true;

                Info<< "\n[DIAGNOSTIC] Natural terminal velocity detected:" << nl
                    << "    Vt       = " << mag(frameVelocity.value()) << " m/s" << nl
                    << "    altitude = " << altitude.value() << " m" << nl
                    << "    rho      = " << rhoAtmosphere.value() << " kg/m3" << nl
                    << "    a_sound  = " << speedOfSound.value() << " m/s" << nl
                    << "    Mach     = " << machNumber.value() << nl
                    << "    Faero    = " << aeroForce.value() << " N" << nl
                    << "    Fnet     = " << netForce.value() << " N" << nl << endl;

                if (terminalVelocityWrite)
                {
                    triggerTerminalVelocityWrite = true;
                }
            }
        }
        else
        {
            terminalVelocityTriggered = false;
        }
    }
}

// ------------------------------------------------------------
//  6. Optional console stdout printing based on print sub-dictionary
// ------------------------------------------------------------
if (fallingFrameDict.found("print"))
{
    const dictionary& printDict = fallingFrameDict.subDict("print");

    if (printDict.getOrDefault<bool>("frameVel", true))
    {
        Info<< "V = " << mag(frameVelocity.value()) << " m/s" << endl;
    }
    if (printDict.getOrDefault<bool>("frameAcc", true))
    {
        Info<< "a = " << mag(frameAcceleration.value()) << " m/s2" << endl;
    }
    if (printDict.getOrDefault<bool>("altitude", true))
    {
        Info<< "Altitude = " << altitude.value() << " m" << endl;
    }
    if (printDict.getOrDefault<bool>("meanFlow", false))
    {
        Info<< "Umean = " << meanFlowVelocity.value() << " m/s" << endl;
    }
    if (printDict.getOrDefault<bool>("relativeVel", false))
    {
        Info<< "Urel = " << relativeVelocity.value() << " m/s" << endl;
    }
    if (printDict.getOrDefault<bool>("speedOfSound", true))
    {
        Info<< "Speed of Sound = " << speedOfSound.value() << " m/s" << endl;
    }
    if (printDict.getOrDefault<bool>("machNumber", true))
    {
        Info<< "Mach Number = " << machNumber.value() << endl;
    }
    #if defined(FF_PIMPLE_COMPRESSIBLE)
    if (printDict.getOrDefault<bool>("maxLocalMachNumber", true))
    {
        Info<< "Max Local Mach = " << maxLocalMachNumber.value() << endl;
    }
    #endif
    if (printDict.getOrDefault<bool>("dynamicPressure", false))
    {
        const scalar qDyn = 0.5 * rhoAtmosphere.value() * magSqr(relativeVelocity.value());
        Info<< "Dynamic Pressure = " << qDyn << " Pa" << endl;
    }
    if (printDict.getOrDefault<bool>("pressureForce", false))
    {
        Info<< "Fpressure = " << pressureForce.value() << " N" << endl;
    }
    if (printDict.getOrDefault<bool>("viscousForce", false))
    {
        Info<< "Fviscous = " << viscousForce.value() << " N" << endl;
    }
    if (printDict.getOrDefault<bool>("gravityForce", true))
    {
        Info<< "Fgravity = " << gravityForce.value().y() << " N" << endl;
    }
    if (printDict.getOrDefault<bool>("aeroForce", true))
    {
        Info<< "Faero = " << (aeroForce.value().y() >= 0 ? "+" : "") << aeroForce.value().y() << " N" << endl;
    }
    if (printDict.getOrDefault<bool>("netForce", true))
    {
        Info<< "Fnet = " << netForce.value() << " N" << endl;
    }
}

// ------------------------------------------------------------
//  7. Check ground interrupt
// ------------------------------------------------------------
if (groundInterrupt && (altitude.value() <= groundAltitude) && (runTime.timeIndex() > 0))
{
    Info<< "\n[GROUND INTERRUPT] Body reached ground altitude ("
        << altitude.value() << " <= " << groundAltitude << " m). Terminating simulation." << nl << endl;
    runTime.writeNow();
    runTime.setTime(runTime.endTime().value(), runTime.timeIndex());
}

// ------------------------------------------------------------
//  8. Write postProcessing motion data (unified 33 columns)
//
//     Controlled by "writeMotionLog" (default: true) in fallingFrameDict.
//     Set to false when using the fallingMotion function object to
//     avoid duplicate writes to postProcessing/fallingMotion/motion.dat.
// ------------------------------------------------------------
if (Pstream::master() && fallingFrameDict.getOrDefault<bool>("writeMotionLog", true))
{
    fileName motionDir = runTime.path()/"postProcessing"/"fallingMotion";
    mkDir(motionDir);

    static autoPtr<OFstream> motionFilePtr;
    if (!motionFilePtr.valid())
    {
        const fileName motionPath = motionDir/"motion.dat";
        const bool append = (runTime.timeIndex() > runTime.startTimeIndex()) || (isFile(motionPath) && runTime.value() > 0);

        motionFilePtr.reset
        (
            new OFstream
            (
                motionPath,
                IOstreamOption(),
                (append ? IOstreamOption::APPEND : IOstreamOption::NO_APPEND)
            )
        );

        if (!append)
        {
            motionFilePtr()
                << "# Time Altitude Vx Vy Vz Ax Ay Az g_x g_y g_z rho "
                << "Umean_x Umean_y Umean_z Urel_x Urel_y Urel_z "
                << "Fpress_x Fpress_y Fpress_z Fvisc_x Fvisc_y Fvisc_z "
                << "Faero_x Faero_y Faero_z Fgrav_x Fgrav_y Fgrav_z "
                << "Fnet_x Fnet_y Fnet_z" << endl;
        }
    }

    OFstream& motionFile = motionFilePtr();
    motionFile << runTime.timeName() << " "
               << altitude.value() << " "
               << frameVelocity.value().x() << " " << frameVelocity.value().y() << " " << frameVelocity.value().z() << " "
               << frameAcceleration.value().x() << " " << frameAcceleration.value().y() << " " << frameAcceleration.value().z() << " "
               << gravity.value().x() << " " << gravity.value().y() << " " << gravity.value().z() << " "
               << rhoAtmosphere.value() << " "
               << meanFlowVelocity.value().x() << " " << meanFlowVelocity.value().y() << " " << meanFlowVelocity.value().z() << " "
               << relativeVelocity.value().x() << " " << relativeVelocity.value().y() << " " << relativeVelocity.value().z() << " "
               << pressureForce.value().x() << " " << pressureForce.value().y() << " " << pressureForce.value().z() << " "
               << viscousForce.value().x() << " " << viscousForce.value().y() << " " << viscousForce.value().z() << " "
               << aeroForce.value().x() << " " << aeroForce.value().y() << " " << aeroForce.value().z() << " "
               << gravityForce.value().x() << " " << gravityForce.value().y() << " " << gravityForce.value().z() << " "
               << netForce.value().x() << " " << netForce.value().y() << " " << netForce.value().z() << endl;
}
