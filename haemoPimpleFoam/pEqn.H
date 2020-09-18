#ifdef OPENFOAMESIORFOUNDATION
volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
surfaceScalarField phiHbyA
(
"phiHbyA",
fvc::flux(HbyA)
+ MRF.zeroFilter(fvc::interpolate(rAU)*fvc::ddtCorr(U, phi, Uf))
);

MRF.makeRelative(phiHbyA);

if (p.needReference())
{
fvc::makeRelative(phiHbyA, U);
adjustPhi(phiHbyA, U, p);
fvc::makeAbsolute(phiHbyA, U);
}

tmp<volScalarField> rAtU(rAU);

if (pimple.consistent())
{
rAtU = 1.0/max(1.0/rAU - UEqn.H1(), 0.1/rAU);
phiHbyA +=
    fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
HbyA -= (rAU - rAtU())*fvc::grad(p);
}

if (pimple.nCorrPiso() <= 1)
{
tUEqn.clear();
}

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAtU(), MRF);

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
fvScalarMatrix pEqn
(
    fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
);

pEqn.setReference(pRefCell, pRefValue);

pEqn.solve();

if (pimple.finalNonOrthogonalIter())
{
    phi = phiHbyA - pEqn.flux();
}
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
    p.relax();
    
U = HbyA - rAtU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);

// Correct Uf if the mesh is moving
fvc::correctUf(Uf, U, phi);

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi, U);
#else
    {
    p.boundaryField().updateCoeffs();

    // Prepare clean Ap without time derivative contribution and
    // without contribution from under-relaxation
    // HJ, 26/Oct/2015
    aU = HUEqn.A();

    // Store velocity under-relaxation point before using U for the flux
    // precursor
    U.storePrevIter();

    U = HUEqn.H()/aU;
    phi = (fvc::interpolate(U) & mesh.Sf());

    // Global flux balance
    adjustPhi(phi, U, p);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
            (
             fvm::laplacian(1/aU, p) == fvc::div(phi)
            );

        pEqn.setReference(pRefCell, pRefValue);
        pEqn.solve
            (
             mesh.solutionDict().solver(p.select(pimple.finalInnerIter()))
            );

        if (pimple.finalNonOrthogonalIter())
        {
            phi -= pEqn.flux();
        }
    }

    // Explicitly relax pressure for momentum corrector except for last corrector
    if (!pimple.finalIter())
    {
        p.relax();
    }

#   include "continuityErrs.H"

    U = UUrf*
        (
         1.0/(aU + ddtUEqn.A())*
         (
          U*aU - fvc::grad(p) + ddtUEqn.H()
         )
        )
        + (1 - UUrf)*U.prevIter();
    U.correctBoundaryConditions();
    }
#endif
