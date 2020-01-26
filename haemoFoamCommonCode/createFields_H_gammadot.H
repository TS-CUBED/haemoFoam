    // This is the additional scalar for haematocrit concentration

    Info<< "Reading field H\n" <<endl;
    volScalarField H
        (
         IOobject
         (
          "H",
          runTime.timeName(),
          mesh,
          IOobject::MUST_READ,
          IOobject::AUTO_WRITE
         ),
         mesh
        );

    Info<< "Creating field gammaDot\n" <<endl;
    volScalarField gammaDot
        (
         IOobject
         (
          "gammaDot",
          runTime.timeName(),
          mesh,
          IOobject::READ_IF_PRESENT,
          IOobject::AUTO_WRITE
         ),
         mesh,
         dimensionedScalar("gammaDot", dimless/dimTime, scalar(0.0))
        );

    //  and back to the original code
