/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    writeFFDict

Description
    Utility to generate a constant/fallingFrameDict dictionary file for
    freeFallingFoam cases with full high-accuracy defaults and comments.
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "OFstream.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Generate a constant/fallingFrameDict file for freeFallingFoam simulations."
    );

    argList args(argc, argv);

    fileName casePath = args.path();
    fileName dictPath = casePath / "constant" / "fallingFrameDict";

    Info<< "Writing fallingFrameDict to " << dictPath << endl;

    mkDir(casePath / "constant");

    OFstream os(dictPath);

    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot open file " << dictPath << " for writing"
            << exit(FatalError);
    }

    os  << "/*--------------------------------*- C++ -*----------------------------------*\\\n"
        << "| =========                 |                                                 |\n"
        << "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n"
        << "|  \\\\    /   O peration     | Version:  v2412                                 |\n"
        << "|   \\\\  /    A nd           | Website:  www.openfoam.com                      |\n"
        << "|    \\\\/     M anipulation  |                                                 |\n"
        << "\\*---------------------------------------------------------------------------*/\n"
        << "FoamFile\n"
        << "{\n"
        << "    version     2.0;\n"
        << "    format      ascii;\n"
        << "    class       dictionary;\n"
        << "    location    \"constant\";\n"
        << "    object      fallingFrameDict;\n"
        << "}\n"
        << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n"
        << "// -----------------------------------------------------------------------------\n"
        << "// 1. Initial State & Mass Properties\n"
        << "// -----------------------------------------------------------------------------\n"
        << "u0                  ( 0 0 0 );          // Initial frame velocity [m/s]\n"
        << "acceleration        ( 0 -9.81 0 );      // Gravity acceleration vector [m/s^2]\n"
        << "freedom             ( 0 1 0 );          // Motion mask: (0 1 0) for 1D drop, (1 1 1) for 3D\n"
        << "mass                98018;              // Body mass [kg]\n"
        << "patches             ( sphere );         // Wall patch(es) for aerodynamic force integration\n\n\n"
        << "// -----------------------------------------------------------------------------\n"
        << "// 2. High-Fidelity Atmospheric & Gravity Models\n"
        << "// -----------------------------------------------------------------------------\n"
        << "rhoInf              1.204;              // Sea-level reference density [kg/m^3]\n"
        << "densityModel        altitudeBased;      // Barometric density variation with altitude (options: rhoInf, altitudeBased)\n"
        << "gravityModel        altitudeBased;      // Gravity variation g(z) with altitude (options: constant, altitudeBased)\n\n\n"
        << "// -----------------------------------------------------------------------------\n"
        << "// 3. Trajectory & Ground Collision Controls\n"
        << "// -----------------------------------------------------------------------------\n"
        << "initialAltitude     10000;              // Drop start altitude [m] (e.g. 10 km)\n"
        << "groundAltitude      0;                  // Target ground altitude [m]\n"
        << "groundInterrupt     on;                 // Automatically write & stop at ground impact (options: on, off)\n\n\n"
        << "// -----------------------------------------------------------------------------\n"
        << "// 4. Natural Terminal Velocity & Autosave Toggles\n"
        << "// -----------------------------------------------------------------------------\n"
        << "terminalAccelerationTolerance  1e-3;   // Acceleration threshold for steady terminal state [m/s^2]\n"
        << "terminalVelocityTolerance      1e-3;   // Velocity variation tolerance over sampling window [m/s]\n"
        << "terminalVelocitySamples        20;     // Number of timesteps in moving window\n"
        << "terminalVelocityWrite          on;     // Force-save 3D fields when terminal velocity is reached (options: on, off)\n"
        << "terminalVelocityStop           off;    // Set 'on' to stop run on terminal velocity, or 'off' to continue\n\n\n"
        << "// -----------------------------------------------------------------------------\n"
        << "// 5. Mach Regime & Dynamic Transonic Controls (tFFRhoPimpleFoam)\n"
        << "// -----------------------------------------------------------------------------\n"
        << "machRegime              auto;          // Mode: subsonic (transonic off), transonic (transonic on), auto (auto switch/stop)\n"
        << "machRegimeSwitchCutOff  0.7;           // Peak local Mach threshold for regime trigger (default: 0.7)\n"
        << "machRegimeAction        switch;        // Action in auto mode: 'switch' (dynamically enable transonic) or 'stop' (halt simulation)\n\n\n"
        << "// -----------------------------------------------------------------------------\n"
        << "// 6. Tabulated Atmospheric Crosswind Profile\n"
        << "// -----------------------------------------------------------------------------\n"
        << "meanFlow\n"
        << "{\n"
        << "    model               tabulated;      // Mean flow model (options: constant, powerLaw, tabulated)\n"
        << "    table\n"
        << "    (\n"
        << "        (     0.0       (  5.0  0.0  0.0 ) )   // 5 m/s wind at ground level (z = 0 m)\n"
        << "        (  1000.0       ( 12.0  0.0  0.0 ) )   // 12 m/s wind at 1 km altitude\n"
        << "        (  5000.0       ( 25.0  0.0  0.0 ) )   // 25 m/s wind at 5 km altitude\n"
        << "        ( 10000.0       ( 40.0  0.0  0.0 ) )   // 40 m/s jet stream at 10 km altitude\n"
        << "    );\n"
        << "}\n\n\n"
        << "// -----------------------------------------------------------------------------\n"
        << "// 7. Console Output Logging\n"
        << "// -----------------------------------------------------------------------------\n"
        << "print\n"
        << "{\n"
        << "    frameVel            on;\n"
        << "    frameAcc            on;\n"
        << "    altitude            on;\n"
        << "    pressureForce       off;\n"
        << "    viscousForce        off;\n"
        << "    aeroForce           on;\n"
        << "    gravityForce        on;\n"
        << "    netForce            on;\n"
        << "    meanFlow            on;\n"
        << "    relativeVel         on;\n"
        << "    speedOfSound        on;\n"
        << "    machNumber          on;\n"
        << "    maxLocalMachNumber  on;\n"
        << "    dynamicPressure     off;\n"
        << "}\n\n"
        << "// ************************************************************************* //\n";

    Info<< "Successfully generated " << dictPath << endl;

    return 0;
}

// ************************************************************************* //
