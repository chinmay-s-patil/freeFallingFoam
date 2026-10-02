// ============================================================
//  updateFallingFrameMotion.H — tFFRhoCentralFoam
// ============================================================

#include "OFstream.H"

// ------------------------------------------------------------
//  1. Aerodynamic force integration from static pressure p and devRhoReff
// ------------------------------------------------------------
vector fPress = vector::zero;
vector fVisc = vector::zero;

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
            // Pressure force: p is static pressure [Pa] in compressible solver.
            // Sf points outwards from fluid domain into body wall.
            fPress += p_p[facei] * Sf[facei];

            // Viscous shear force exerted on body wall by fluid = (devTau & Sf)
            fVisc += (devTau_p[facei] & Sf[facei]);
        }
    }
}

reduce(fPress, sumOp<vector>());
reduce(fVisc, sumOp<vector>());

pressureForce.value() = fPress;
viscousForce.value()  = fVisc;
aeroForce.value()     = fPress + fVisc;

// ------------------------------------------------------------
//  2. Force balance and explicit acceleration calculation
// ------------------------------------------------------------
gravityForce.value() = bodyMass * gravity.value();
netForce.value()     = gravityForce.value() + aeroForce.value();

// Acceleration at timestep n
vector a_n = netForce.value() / bodyMass;

// Apply motion freedom constraint (e.g. freedom = (0 1 0) for 1D drop)
a_n.x() *= freedomMask.x();
a_n.y() *= freedomMask.y();
a_n.z() *= freedomMask.z();

frameAcceleration.value() = a_n;

// ------------------------------------------------------------
//  3. Integrate velocity and altitude for timestep n+1
// ------------------------------------------------------------
const scalar dt = runTime.deltaTValue();
frameVelocity.value() += a_n * dt;

// Determine vertical direction from gravity vector (-g / |g|)
vector gDir = gravity.value();
const scalar magG = mag(gDir);
if (magG > SMALL)
{
    gDir /= magG;
}
else
{
    gDir = vector(0, -1, 0);
}

// Vertical velocity component relative to ground (upward = -gDir)
const scalar V_up = -(frameVelocity.value() & gDir);
altitude.value() += V_up * dt;

// Evaluate mean atmospheric flow U_mean(z)
vector Umean = vector::zero;
if (fallingFrameDict.found("meanFlow"))
{
    const dictionary& meanFlowDict = fallingFrameDict.subDict("meanFlow");
    const word model = meanFlowDict.getOrDefault<word>("model", "constant");

    if (model == "constant")
    {
        Umean = meanFlowDict.getOrDefault<vector>("velocity", vector::zero);
    }
    else if (model == "powerLaw")
    {
        const scalar zRef = meanFlowDict.getOrDefault<scalar>("referenceHeight", 10.0);
        const vector Uref = meanFlowDict.getOrDefault<vector>("referenceVelocity", vector::zero);
        const scalar alpha = meanFlowDict.getOrDefault<scalar>("exponent", 0.14);
        const scalar minH  = meanFlowDict.getOrDefault<scalar>("minHeight", 0.1);

        const scalar currentH = max(altitude.value(), minH);
        const scalar ratio = max(currentH / max(zRef, SMALL), SMALL);
        Umean = Uref * Foam::pow(ratio, alpha);
    }
    else if (model == "tabulated")
    {
        List<Tuple2<scalar, vector>> table;
        if (meanFlowDict.readIfPresent("table", table) && table.size() > 0)
        {
            const scalar zCurr = altitude.value();
            if (zCurr <= table.first().first())
            {
                Umean = table.first().second();
            }
            else if (zCurr >= table.last().first())
            {
                Umean = table.last().second();
            }
            else
            {
                for (label i = 0; i < table.size() - 1; ++i)
                {
                    if (zCurr >= table[i].first() && zCurr <= table[i+1].first())
                    {
                        const scalar z0_val = table[i].first();
                        const scalar z1_val = table[i+1].first();
                        const vector& u0_val = table[i].second();
                        const vector& u1_val = table[i+1].second();
                        const scalar frac = (zCurr - z0_val) / max(z1_val - z0_val, SMALL);
                        Umean = u0_val + frac * (u1_val - u0_val);
                        break;
                    }
                }
            }
        }
    }
}

meanFlowVelocity.value() = Umean;

// Update relative freestream velocity field: U_rel = U_mean - V_frame
relativeVelocity.value() = meanFlowVelocity.value() - frameVelocity.value();

// Update speed of sound and Mach number using cell-averaged pressure and density
{
    const scalar gammaCoeff = 1.4;
    const scalar pMean = max(gAverage(p), 1.0);
    const scalar rhoMean = max(gAverage(rho), SMALL);
    const scalar a_sound = Foam::sqrt(gammaCoeff * (pMean / rhoMean));
    speedOfSound.value() = a_sound;

    const scalar vRelMag = mag(relativeVelocity.value());
    machNumber.value() = vRelMag / max(a_sound, SMALL);
}

// ------------------------------------------------------------
//  4. Update pseudo-force
// ------------------------------------------------------------
fallingFrameForce.value() = -frameAcceleration.value();

// Update boundary conditions with new frame velocity
U.correctBoundaryConditions();

// Update magUStar for post-processing scale
{
    const scalar Uref = max(mag(frameVelocity.value()), 1.0);
    magUStar = mag(U) / dimensionedScalar("Uref", dimVelocity, Uref);
}

// ------------------------------------------------------------
//  5. Terminal velocity diagnostic check
// ------------------------------------------------------------
{
    static DynamicList<scalar> velHistory;
    static bool terminalVelocityTriggered = false;

    velHistory.append(mag(frameVelocity.value()));
    if (velHistory.size() > terminalVelocitySamples)
    {
        velHistory.remove(0);
    }

    if (velHistory.size() >= terminalVelocitySamples)
    {
        scalar minV = velHistory[0];
        scalar maxV = velHistory[0];
        forAll(velHistory, vi)
        {
            minV = min(minV, velHistory[vi]);
            maxV = max(maxV, velHistory[vi]);
        }
        scalar vVar = maxV - minV;

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
//  5b. Optional console stdout printing based on print sub-dictionary
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
    if (printDict.getOrDefault<bool>("speedOfSound", true))
    {
        Info<< "Speed of Sound = " << speedOfSound.value() << " m/s" << endl;
    }
    if (printDict.getOrDefault<bool>("machNumber", true))
    {
        Info<< "Mach Number = " << machNumber.value() << endl;
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
//  6. Ground interrupt check
// ------------------------------------------------------------
if (groundInterrupt && altitude.value() <= groundAltitude)
{
    Info<< "\n[GROUND INTERRUPT] Body reached ground altitude ("
        << altitude.value() << " <= " << groundAltitude << " m). Terminating simulation." << nl << endl;
    runTime.writeNow();
    runTime.setTime(runTime.endTime().value(), runTime.timeIndex());
}

// ------------------------------------------------------------
//  7. Write postProcessing motion data
// ------------------------------------------------------------
if (Pstream::master())
{
    fileName motionDir = runTime.path()/"postProcessing"/"fallingMotion";
    mkDir(motionDir);

    static autoPtr<OFstream> motionFilePtr;
    if (!motionFilePtr.valid())
    {
        motionFilePtr.reset(new OFstream(motionDir/"motion.dat"));
        motionFilePtr() << "# Time(s) Altitude(m) V_frame_x V_frame_y V_frame_z "
                        << "a_frame_x a_frame_y a_frame_z "
                        << "F_press_x F_press_y F_press_z "
                        << "F_visc_x F_visc_y F_visc_z "
                        << "F_aero_x F_aero_y F_aero_z "
                        << "F_net_x F_net_y F_net_z" << endl;
    }

    OFstream& motionFile = motionFilePtr();
    motionFile << runTime.timeName() << " "
               << altitude.value() << " "
               << frameVelocity.value().x() << " "
               << frameVelocity.value().y() << " "
               << frameVelocity.value().z() << " "
               << frameAcceleration.value().x() << " "
               << frameAcceleration.value().y() << " "
               << frameAcceleration.value().z() << " "
               << pressureForce.value().x() << " "
               << pressureForce.value().y() << " "
               << pressureForce.value().z() << " "
               << viscousForce.value().x() << " "
               << viscousForce.value().y() << " "
               << viscousForce.value().z() << " "
               << aeroForce.value().x() << " "
               << aeroForce.value().y() << " "
               << aeroForce.value().z() << " "
               << netForce.value().x() << " "
               << netForce.value().y() << " "
               << netForce.value().z() << endl;
}
