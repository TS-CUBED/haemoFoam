#ifdef OPENFOAMESIORFOUNDATION
{
    volScalarField rAU(1.0/UEqn.A());
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
    MRF.makeRelative(phiHbyA);
    adjustPhi(phiHbyA, U, p);

    tmp<volScalarField> rAtU(rAU);

    if (simple.consistent())
    {
        rAtU = 1.0/(1.0/rAU - UEqn.H1());
        phiHbyA +=
            fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
        HbyA -= (rAU - rAtU())*fvc::grad(p);
    }

    tUEqn.clear();

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, U, phiHbyA, rAtU(), MRF);

    // Non-orthogonal pressure corrector loop
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux();
        }
    }

    #include "continuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    // Momentum corrector
    U = HbyA - rAtU()*fvc::grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U);
}

#else
{
    p.boundaryField().updateCoeffs();

    // Prepare clean 1/Ap without contribution from under-relaxation
    // HJ, 26/Oct/2015
    volScalarField rAU
        (
         "(1|A(U))",
         1/HUEqn().A()
        );

    // Store velocity under-relaxation point before using U for
    // the flux precursor
    U.storePrevIter();

    U = rAU*HUEqn().H();
    HUEqn.clear();

    phi = fvc::interpolate(U) & mesh.Sf();
    adjustPhi(phi, U, p);

    // Non-orthogonal pressure corrector loop
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
            (
             fvm::laplacian(rAU, p) == fvc::div(phi)
            );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi -= pEqn.flux();
        }
    }

#   include "continuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    // Momentum corrector
    // Note: since under-relaxation does not change aU, H/a in U can be
    // re-used.  HJ, 22/Jan/2016
    U = UUrf*(U - rAU*fvc::grad(p)) + (1 - UUrf)*U.prevIter();
    U.correctBoundaryConditions();
}
#endif
