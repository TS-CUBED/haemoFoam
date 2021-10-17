/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1812                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*\
Description
    An incompressible non-Newtonian viscosity
    model following

    K. K. Yeleswarapu, M. V. Kameneva, K. R. Rajagopal, and J. F. Antaki,
    ‘The flow of blood in tubes:theory and experiment’,
    Mechanics Research Communications, vol. 25, no. 3, pp. 257–262, May 1998.

    W.-T. Wu, F. Yang, J. F. Antaki, N. Aubry, and M. Massoudi,
    ‘Study of blood flow in several benchmark micro-channels using a two-fluid approach’,
    Int J Eng Sci, vol. 95, pp. 49–59, Oct. 2015.

    Author and Copyright
    Dr Torsten Schenkel
    Department Engineering and Mathematics
    Material and Engineering Research Institute MERI
    Sheffield Hallam University
    April 2021
    All Rights Reserved
\*---------------------------------------------------------------------------*/

#include "Yeleswarapu.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace viscosityModels
    {
        defineTypeNameAndDebug(Yeleswarapu, 0);
        addToRunTimeSelectionTable
            (
             viscosityModel,
             Yeleswarapu,
             dictionary
            );
    }
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::Yeleswarapu::calcNu() const
{
    const volScalarField& H= U_.mesh().lookupObject<volScalarField>("H");

//    const volScalarField& nu_0 = a3 * pow(H,3) + a2 * pow(H,2) + a1 * H;
//    const volScalarField& nu_inf = b3 * pow(H,3) + b2 * pow(H,2) + b1 * H;
    
//    const volScalarField& nuS = (nu_inf + (nu_0 - nu_inf)* (1+log(1+k*strainRate()))/(1+k*strainRate()));
    return muPlasma_/rho_ * 
            (1-H) + H * 
            muPlasma_/muPlasma_.value() / 
            (rho_/rho_.value()) *
            ((b3_ * pow(H,3) + b2_ * pow(H,2) + b1_ * H) + 
            ((a3_ * pow(H,3) + a2_ * pow(H,2) + a1_ * H) - 
            (b3_ * pow(H,3) + b2_ * pow(H,2) + b1_ * H)) * 
            (1+log(1+k_*strainRate()))/(1+k_*strainRate()));  
            // + nuS * H * nuPlasma/nuPlasma.value();                                 
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    Foam::viscosityModels::Yeleswarapu::Yeleswarapu
(
 const word& name,
 const dictionary& viscosityProperties,
 const volVectorField& U,
 const surfaceScalarField& phi
 )
    :
        viscosityModel(name, viscosityProperties, U, phi),
        YeleswarapuCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),

#ifdef OPENFOAMESIORFOUNDATION
        a1_("a1", dimless, YeleswarapuCoeffs_),
        a2_("a2", dimless, YeleswarapuCoeffs_),
        a3_("a3", dimless, YeleswarapuCoeffs_),
        b1_("b1", dimless, YeleswarapuCoeffs_),
        b2_("b2", dimless, YeleswarapuCoeffs_),
        b3_("b3", dimless, YeleswarapuCoeffs_),

        k_("k", dimless, YeleswarapuCoeffs_),

        muPlasma_("muPlasma", dimViscosity*dimDensity, YeleswarapuCoeffs_),
        
        rho_("rho", dimDensity, viscosityProperties),
#else
        a1_(YeleswarapuCoeffs_.lookup("a1")),
        a2_(YeleswarapuCoeffs_.lookup("a2")),
        a3_(YeleswarapuCoeffs_.lookup("a3")),
        b1_(YeleswarapuCoeffs_.lookup("b1")),
        b2_(YeleswarapuCoeffs_.lookup("b2")),
        b3_(YeleswarapuCoeffs_.lookup("b3")),

        k_(YeleswarapuCoeffs_.lookup("k")),

        muPlasma_(YeleswarapuCoeffs_.lookup("muPlasma")),

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

    bool Foam::viscosityModels::Yeleswarapu::read
(
 const dictionary& viscosityProperties
 )
{
    viscosityModel::read(viscosityProperties);

    YeleswarapuCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    YeleswarapuCoeffs_.lookup("a1") >> a1_;
    YeleswarapuCoeffs_.lookup("a2") >> a2_;
    YeleswarapuCoeffs_.lookup("a3") >> a3_;
    YeleswarapuCoeffs_.lookup("b1") >> b1_;
    YeleswarapuCoeffs_.lookup("b2") >> b2_;
    YeleswarapuCoeffs_.lookup("b3") >> b3_;
    
    YeleswarapuCoeffs_.lookup("k") >> k_;

    YeleswarapuCoeffs_.lookup("muPlasma") >> muPlasma_;

    // YeleswarapuCoeffs_.lookup("rho") >> rho_;

    viscosityProperties.lookup("rho") >> rho_;

    return true;
}


// ************************************************************************* //
