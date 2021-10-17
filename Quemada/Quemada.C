/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1812                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*\
Description
    Non-Newtonian model for blood based on:

    QUEMADA, D., 1978. Rheology of concentrated disperse systems III.
        General features of the proposed non-newtonian model.
        Comparison with experimental data.
        Rheologica Acta, 17(6), pp. 643-653.

Author and Copyright
    Dr Torsten Schenkel
    Department Engineering and Mathematics
    Material and Engineering Research Institute MERI
    Sheffield Hallam University
    February 2019
    All Rights Reserved
\*---------------------------------------------------------------------------*/

#include "Quemada.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace viscosityModels
    {
        defineTypeNameAndDebug(Quemada, 0);
        addToRunTimeSelectionTable
            (
             viscosityModel,
             Quemada,
             dictionary
            );
    }
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::Quemada::calcNu() const
{
    const volScalarField& H= U_.mesh().lookupObject<volScalarField>("H");

    return muPlasma_/rho_ *
        pow
        (
         (
          1.0-0.5*
          // K: Quemada Parameter
          (
           (
            (a0_+2/(a1_+H))                                   // k0
            +
            exp(b0_+b1_*H+b2_*pow(H,2.0)+b3_*pow(H,3.0))        // kInf
            *sqrt
            (
             strainRate()/
             (
              exp(c0_+c1_*H+c2_*pow(H,2.0)+c3_*pow(H,3.0))
              *gammaC0_                                        // gammaC
             )
            )
           )
           /
           (
            1+sqrt
            (
             strainRate()/
             (
              exp(c0_+c1_*H+c2_*pow(H,2.0)+c3_*pow(H,3.0))
              *gammaC0_                                        // gammaC
             )
            )
           )
           )
           // end K
           *H
           ),-2
           )
           ;
}





// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //



    Foam::viscosityModels::Quemada::Quemada
(
 const word& name,
 const dictionary& viscosityProperties,
 const volVectorField& U,
 const surfaceScalarField& phi
 )
    :
        viscosityModel(name, viscosityProperties, U, phi),
        QuemadaCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),

#ifdef OPENFOAMESIORFOUNDATION
        a0_("a0", dimless, QuemadaCoeffs_),
        a1_("a1", dimless, QuemadaCoeffs_),
        b0_("b0", dimless, QuemadaCoeffs_),
        b1_("b1", dimless, QuemadaCoeffs_),
        b2_("b2", dimless, QuemadaCoeffs_),
        b3_("b3", dimless, QuemadaCoeffs_),
        c0_("c0", dimless, QuemadaCoeffs_),
        c1_("c1", dimless, QuemadaCoeffs_),
        c2_("c2", dimless, QuemadaCoeffs_),
        c3_("c3", dimless, QuemadaCoeffs_),

        gammaC0_("gammaC0", dimless/dimTime, QuemadaCoeffs_),
        muPlasma_("muPlasma", dimViscosity*dimDensity, QuemadaCoeffs_),
        
        rho_("rho", dimDensity, viscosityProperties),
#else
        a0_(QuemadaCoeffs_.lookup("a0")),
        a1_(QuemadaCoeffs_.lookup("a1")),
        b0_(QuemadaCoeffs_.lookup("b0")),
        b1_(QuemadaCoeffs_.lookup("b1")),
        b2_(QuemadaCoeffs_.lookup("b2")),
        b3_(QuemadaCoeffs_.lookup("b3")),
        c0_(QuemadaCoeffs_.lookup("c0")),
        c1_(QuemadaCoeffs_.lookup("c1")),
        c2_(QuemadaCoeffs_.lookup("c2")),
        c3_(QuemadaCoeffs_.lookup("c3")),

        gammaC0_(QuemadaCoeffs_.lookup("gammaC0")),
        muPlasma_(QuemadaCoeffs_.lookup("muPlasma")),

        rho_(viscosityProperties.lookup("rho")),
#endif


        nu_
        (
         IOobject
         (
          "nu",
          U_.time().timeName(),
          U_.db(),
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
         ),
         calcNu()
        )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

    bool Foam::viscosityModels::Quemada::read
(
 const dictionary& viscosityProperties
 )
{
    viscosityModel::read(viscosityProperties);

    QuemadaCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    QuemadaCoeffs_.lookup("a0") >> a0_;
    QuemadaCoeffs_.lookup("a1") >> a1_;
    QuemadaCoeffs_.lookup("b0") >> b0_;
    QuemadaCoeffs_.lookup("b1") >> b1_;
    QuemadaCoeffs_.lookup("b2") >> b2_;
    QuemadaCoeffs_.lookup("b3") >> b3_;
    QuemadaCoeffs_.lookup("c0") >> c0_;
    QuemadaCoeffs_.lookup("c1") >> c1_;
    QuemadaCoeffs_.lookup("c2") >> c2_;
    QuemadaCoeffs_.lookup("c3") >> c3_;

    QuemadaCoeffs_.lookup("gammaC0") >> gammaC0_;
    QuemadaCoeffs_.lookup("muPlasma") >> muPlasma_;

    // QuemadaCoeffs_.lookup("rho") >> rho_;

    viscosityProperties.lookup("rho") >> rho_;

    return true;
}


// ************************************************************************* //
