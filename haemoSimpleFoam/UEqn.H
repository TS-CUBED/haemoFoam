    // Solve the momentum equation
#ifdef OPENFOAMESIORFOUNDATION
    MRF.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (simple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));

        fvOptions.correct(U);
    }

#else
    tmp<fvVectorMatrix> HUEqn
        (
         fvm::div(phi, U)
         + turbulence->divDevReff()
        );

    // Get under-relaxation factor
    const scalar UUrf = mesh.solutionDict().equationRelaxationFactor(U.name());

    // Momentum solution
    solve
        (
         relax(HUEqn(), UUrf)
         ==
         -fvc::grad(p)
        );
#endif
