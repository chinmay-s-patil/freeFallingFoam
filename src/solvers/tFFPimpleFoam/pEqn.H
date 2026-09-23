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

// Add hydrostatic head flux contribution
phiHbyA += fvc::interpolate(rAtU())*fvc::snGrad(gh)*mesh.magSf();

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

constrainPressure(p_rgh, U, phiHbyA, rAtU(), MRF);

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

p_rgh.relax();

U = HbyA - rAtU()*fvc::grad(p_rgh) + rAtU()*fallingFrameForce;
U.correctBoundaryConditions();
fvOptions.correct(U);

fvc::correctUf(Uf, U, phi);

fvc::makeRelative(phi, U);
