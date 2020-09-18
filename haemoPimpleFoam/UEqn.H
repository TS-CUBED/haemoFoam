#ifdef OPENFOAMESIORFOUNDATION
    // Solve the Momentum equation

    MRF.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(U) + fvm::div(phi, U)
      + MRF.DDt(U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));

        fvOptions.correct(U);
    }
#else
    // Convection-diffusion matrix
    fvVectorMatrix HUEqn
    (
     fvm::div(phi, U)
     + turbulence->divDevReff()
    );

    // Time derivative matrix
    fvVectorMatrix ddtUEqn(fvm::ddt(U));

    // Get under-relaxation factor
    scalar UUrf =
    mesh.solutionDict().equationRelaxationFactor(U.select(pimple.finalIter()));

    if (pimple.momentumPredictor())
    {
        // Solve momentum predictor
        solve
            (
             ddtUEqn
             + relax(HUEqn, UUrf)
             ==
             - fvc::grad(p),
             mesh.solutionDict().solver((U.select(pimple.finalIter())))
            );
    }
    else
    {
        U = (ddtUEqn.H() + HUEqn.H() - fvc::grad(p))/(HUEqn.A() + ddtUEqn.A());
        U.correctBoundaryConditions();
    }
#endif
