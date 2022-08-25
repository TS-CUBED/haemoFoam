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

    Info<< "Creating field sourceCorig\n" <<endl;
    volScalarField sourceCorig
        (
         IOobject
         (
          "sourceCorig",
          runTime.timeName(),
          mesh,
          IOobject::READ_IF_PRESENT,
          IOobject::AUTO_WRITE
         ),
         mesh,
         dimensionedScalar("sourceCorig", dimless/dimTime, scalar(0.0))
        );
    Info<< "Creating field sourceCHgamma\n" <<endl;
    volScalarField sourceCHgamma
        (
         IOobject
         (
          "sourceCHgamma",
          runTime.timeName(),
          mesh,
          IOobject::READ_IF_PRESENT,
          IOobject::AUTO_WRITE
         ),
         mesh,
         dimensionedScalar("sourceCHgamma", dimless/dimTime, scalar(0.0))
        );
    Info<< "Creating field sourceCHgammaGrad\n" <<endl;
    volScalarField sourceCHgammaGrad
        (
         IOobject
         (
          "sourceCHgammaGrad",
          runTime.timeName(),
          mesh,
          IOobject::READ_IF_PRESENT,
          IOobject::AUTO_WRITE
         ),
         mesh,
         dimensionedScalar("sourceCHgammaGrad", dimless/dimTime, scalar(0.0))
        );
    Info<< "Creating field sourceCHgammaZ\n" <<endl;
    volScalarField sourceCHgammaZ
        (
         IOobject
         (
          "sourceCHgammaZ",
          runTime.timeName(),
          mesh,
          IOobject::READ_IF_PRESENT,
          IOobject::AUTO_WRITE
         ),
         mesh,
         dimensionedScalar("sourceCHgammaZ", dimless/dimTime, scalar(0.0))
        );
    Info<< "Creating field sourceCfvcgrad\n" <<endl;
    volScalarField sourceCgammaZ
        (
         IOobject
         (
          "sourceCgammaZ",
          runTime.timeName(),
          mesh,
          IOobject::READ_IF_PRESENT,
          IOobject::AUTO_WRITE
         ),
         mesh,
         dimensionedScalar("sourceCgammaZ", dimless/dimTime, scalar(0.0))
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

    Info<< "Creating field sourceVorig\n" <<endl;
    volScalarField sourceVorig
        (
         IOobject
         (
          "sourceVorig",
          runTime.timeName(),
          mesh,
          IOobject::READ_IF_PRESENT,
          IOobject::AUTO_WRITE
         ),
         mesh,
         dimensionedScalar("sourceVorig", dimless/dimTime, scalar(0.0))
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

    volTensorField Ei
    (
        IOobject
        (
            "Ei",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedTensor("Ei", dimless, tensor::zero)
    );

    volVectorField Ev
    (
        IOobject
        (
            "Ev",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("Ev", dimless/dimTime, vector::zero)
    );


    volVectorField Evmag
    (
        IOobject
        (
            "Evmag",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("Evmag", dimless/dimTime, vector::zero)
    );


    volVectorField d1
    (
        IOobject
        (
            "d1",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("d1", dimless, vector::zero)
    );

    volVectorField d2
    (
        IOobject
        (
            "d2",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("d2", dimless, vector::zero)
    );

    volVectorField d3
    (
        IOobject
        (
            "d3",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("d3", dimless, vector::zero)
    );


    volVectorField gradGammaDot
    (
        IOobject
        (
            "gradGammaDot",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("gradGammaDot", dimless/dimTime/dimLength, vector::zero)
    );

    volVectorField divGammaDotZ
    (
        IOobject
        (
            "divGammaDotZ",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("divGammaDotZ", dimless/dimTime/dimLength, vector::zero)
    );

    volScalarField minEv
        (
         IOobject
         (
          "minEv",
          runTime.timeName(),
          mesh,
          IOobject::READ_IF_PRESENT,
          IOobject::AUTO_WRITE
         ),
         mesh,
         dimensionedScalar("minEv", dimless/dimTime, scalar(0.0))
        );


    // Info<< "Creating field E\n" <<endl;
    // volTensorField EigenVectors
    //     (
    //      IOobject
    //      (
    //       "Eigenvectors",
    //       runTime.timeName(),
    //       mesh,
    //       IOobject::READ_IF_PRESENT,
    //       IOobject::AUTO_WRITE
    //      ),
    //      mesh,
    //      dimensionedTensor("EigenVectors", dimless/dimTime, tensor::zero)
    //     );
    //  and back to the original code
