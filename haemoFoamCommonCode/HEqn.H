// Haematocrit transport equation


gammaDot = pow(2,0.5)*mag(symm(fvc::grad(U)));
// that is the usual definition of gammaDot

        fvScalarMatrix HEqn
        (
            fvm::ddt(H)
            + fvm::div(phi, H)
            - fvm::laplacian(linKc*Kc*sqr(a)*H*gammaDot, H)
            - fvm::laplacian(nonlinKc*nonlinKcRatio*Kmu*H*sqr(a)*H*gammaDot, H)
            ==
            fvc::laplacian(linKc*Kc*sqr(a)*sqr(H), gammaDot)
            + fvc::laplacian(nonlinKc*nonlinKcRatio*Kmu*H*sqr(a)*sqr(H), gammaDot)
            + fvc::laplacian(Kmu*sqr(a)
            * gammaDot*sqr(H)/laminarTransport.nu(),
            laminarTransport.nu())
        );

        if (haemoSwitch.value() == 0 || runTime < haemoSwitchTime)
        {
            Info<< "Not Solving for H, Migration Model is inactive 1" << nl << endl;
        } else
        {
            HEqn.relax();                   // use these two lines
            HEqn.solve().initialResidual(); // for underrelaxed solver
            // HEqn.solve(); // or this one for no underrelaxation (unstable)
        }
