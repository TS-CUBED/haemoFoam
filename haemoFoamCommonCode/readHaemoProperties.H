    // Variables for the haemo-particle-transport model

    Info << "Reading haemoTransportProperties\n" << endl;

    IOdictionary haemoTransportProperties
        (
         IOobject
         (
          "haemoTransportProperties", // name of the dictionary
          runTime.constant(), // location in the case - this one is in constant
          mesh, // needs the mesh object reference to do some voodoo - unimportant now
          IOobject::MUST_READ, // the file will be re-read if it gets modified during time stepping
          IOobject::NO_WRITE // read-only
         )
        );

#ifdef OPENFOAMESIORFOUNDATION
    dimensionedScalar a   
    ("a",     dimLength, haemoTransportProperties.lookup("a"));
    dimensionedScalar Kc  
    ("Kc",    dimless, haemoTransportProperties.lookup("Kc"));
    dimensionedScalar Kmu 
    ("Kmu",   dimless, haemoTransportProperties.lookup("Kmu"));
    
    dimensionedScalar haemoSwitch 
    ("haemoSwitch",   dimless, haemoTransportProperties.lookup("haemoSwitch"));
    dimensionedScalar haemoSwitchTime
    ("haemoSwitchTime",   dimTime, haemoTransportProperties.lookup("haemoSwitchTime"));
    
    dimensionedScalar haemoSw = haemoSwitch;
    
    dimensionedScalar linKc 
    ("linKc",   dimless, haemoTransportProperties.lookup("linKc"));
    dimensionedScalar nonlinKc
    ("nonlinKc",   dimless, haemoTransportProperties.lookup("nonlinKc"));
    dimensionedScalar nonlinKcRatio 
    ("nonlinKcRatio",   dimless, haemoTransportProperties.lookup("nonlinKcRatio"));
#else
    dimensionedScalar a
        (haemoTransportProperties.lookup("a"));
    dimensionedScalar Kc
        (haemoTransportProperties.lookup("Kc"));
    dimensionedScalar Kmu
        (haemoTransportProperties.lookup("Kmu"));

    dimensionedScalar haemoSwitch
        (haemoTransportProperties.lookup("haemoSwitch"));
    dimensionedScalar haemoSwitchTime
        (haemoTransportProperties.lookup("haemoSwitchTime"));

    dimensionedScalar haemoSw = haemoSwitch;

    dimensionedScalar linKc
        (haemoTransportProperties.lookup("linKc"));
    dimensionedScalar nonlinKc
        (haemoTransportProperties.lookup("nonlinKc"));
    dimensionedScalar nonlinKcRatio
        (haemoTransportProperties.lookup("nonlinKcRatio"));
#endif

    if (linKc.value() == 1.0 && nonlinKc.value() == 0.0) {
        Info << "linear Kc selected" << nl << Kmu
            << nl << Kc << nl << endl;
    }
    else if (linKc.value() == 0.0 && nonlinKc.value() == 1.0) {
        Info << "non-linear Kc selected" << nl << Kmu << nl
            << nonlinKcRatio << nl << endl;
    }
    else {
        Info << "Inconsistent Kc settings! If this is what you want, that's fine." << nl <<
            "But it may not give useful results:" << nl
            << Kmu << nl
            << Kc << nl
            << linKc << nl
            << nonlinKc << nl
            << nonlinKcRatio << nl << endl;
    }


    // Back to original code
