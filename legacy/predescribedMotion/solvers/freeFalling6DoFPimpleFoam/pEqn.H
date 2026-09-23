// ============================================================
//  pEqn.H  (p_rgh formulation for falling reference frame)
//
//  Solves for p_rgh (hydrostatic-reduced kinematic pressure).
//  Total pressure is reconstructed: p = p_rgh + gh
// ============================================================

volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(HbyA)
);

if (pimple.ddtCorr())
{
    phiHbyA += MRF.zeroFilter(fvc::interpolate(rAU)*fvc::ddtCorr(U, phi, Uf));
}

MRF.makeRelative(phiHbyA);

if (p_rgh.needReference())
{
    fvc::makeRelative(phiHbyA, U);
    adjustPhi(phiHbyA, U, p_rgh);
    fvc::makeAbsolute(phiHbyA, U);
}

tmp<volScalarField> rAtU(rAU);

if (pimple.consistent())
{
    rAtU = 1.0/max(1.0/rAU - UEqn.H1(), 0.1/rAU);
    phiHbyA +=
        fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p_rgh)*mesh.magSf();
    HbyA -= (rAU - rAtU())*fvc::grad(p_rgh);
}

// Add the body force flux contribution.
// In the p_rgh formulation, the momentum equation is:
//   UEqn == -grad(p_rgh) - fallingFrameForce
// Uses rAtU (== rAU when consistent mode is off) and snGrad(gh) for
// consistency with the pressure solve and standard buoyant solvers.
phiHbyA += fvc::interpolate(rAtU())*fvc::snGrad(gh)*mesh.magSf();

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_rgh, U, phiHbyA, rAtU(), MRF);

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix p_rghEqn
    (
        fvm::laplacian(rAtU(), p_rgh) == fvc::div(phiHbyA)
    );

    p_rghEqn.setReference(pRefCell, pRefValue);

    p_rghEqn.solve(p_rgh.select(pimple.finalInnerIter()));

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - p_rghEqn.flux();
    }
}

#include "continuityErrs.H"

// Reconstruct total pressure for post-processing
p = p_rgh + gh;
p.correctBoundaryConditions();

// Explicitly relax p_rgh for momentum corrector
p_rgh.relax();

// Velocity correction
U = HbyA - rAtU()*fvc::grad(p_rgh) - rAtU()*fallingFrameForce;
U.correctBoundaryConditions();
fvOptions.correct(U);

// Correct Uf if the mesh is moving
fvc::correctUf(Uf, U, phi);

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi, U);
