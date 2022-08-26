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
          // IOobject::READ_IF_PRESENT,
          IOobject::MUST_READ,
          IOobject::AUTO_WRITE
         ),
         mesh,
         dimensionedScalar("gammaDot", dimless/dimTime, scalar(0.0))
        );

    
    Info<< "Creating field diffH\n" <<endl;
    volScalarField diffH
        (
         IOobject
         (
          "diffH",
          runTime.timeName(),
          mesh,
          IOobject::READ_IF_PRESENT,
          IOobject::AUTO_WRITE
         ),
         mesh,
         dimensionedScalar("diffH", dimless/dimTime, scalar(0.0))
        );
    
    Info<< "Creating field diffHnonlin\n" <<endl;
    volScalarField diffHnonlin
        (
         IOobject
         (
          "diffHnonlin",
          runTime.timeName(),
          mesh,
          IOobject::READ_IF_PRESENT,
          IOobject::AUTO_WRITE
         ),
         mesh,
         dimensionedScalar("diffHnonlin", dimless/dimTime, scalar(0.0))
        );

    Info<< "Creating field sourceC\n" <<endl;
    volScalarField sourceC
        (
         IOobject
         (
          "sourceC",
          runTime.timeName(),
          mesh,
          IOobject::READ_IF_PRESENT,
          IOobject::AUTO_WRITE
         ),
         mesh,
         dimensionedScalar("sourceC", dimless/dimTime, scalar(0.0))
        );
    
    Info<< "Creating field sourceV\n" <<endl;
    volScalarField sourceV
        (
         IOobject
         (
          "sourceV",
          runTime.timeName(),
          mesh,
          IOobject::READ_IF_PRESENT,
          IOobject::AUTO_WRITE
         ),
         mesh,
         dimensionedScalar("sourceV", dimless/dimTime, scalar(0.0))
        );
    
    Info<< "Creating field sourceCnonlin\n" <<endl;
    volScalarField sourceCnonlin
        (
         IOobject
         (
          "sourceCnonlin",
          runTime.timeName(),
          mesh,
          IOobject::READ_IF_PRESENT,
          IOobject::AUTO_WRITE
         ),
         mesh,
         dimensionedScalar("sourceCnonlin", dimless/dimTime, scalar(0.0))
        );
    //  and back to the original code
