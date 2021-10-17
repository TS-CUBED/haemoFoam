/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1812                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*\
Description
    An incompressible Krieger-Dougherty non-Newtonian viscosity model.

Author and Copyright
    Dr Torsten Schenkel
    Department Engineering and Mathematics
    Material and Engineering Research Institute MERI
    Sheffield Hallam University
    February 2019
    All Rights Reserved
\*---------------------------------------------------------------------------*/

#include "KriegerDougherty.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace viscosityModels
    {
        defineTypeNameAndDebug(KriegerDougherty, 0);
        addToRunTimeSelectionTable
            (
             viscosityModel,
             KriegerDougherty,
             dictionary
            );
    }
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::KriegerDougherty::calcNu() const
{
    const volScalarField& H= U_.mesh().lookupObject<volScalarField>("H");

    return muPlasma_/rho_*pow((1-H/Hcrit_),(-1*n_));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    Foam::viscosityModels::KriegerDougherty::KriegerDougherty
(
 const word& name,
 const dictionary& viscosityProperties,
 const volVectorField& U,
 const surfaceScalarField& phi
 )
    :
        viscosityModel(name, viscosityProperties, U, phi),
        KriegerDoughertyCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),

#ifdef OPENFOAMESIORFOUNDATION
        n_("n", dimless, KriegerDoughertyCoeffs_),
        muPlasma_("muPlasma", dimViscosity*dimDensity, KriegerDoughertyCoeffs_),

        Hcrit_("Hcrit", dimless, KriegerDoughertyCoeffs_),
        rho_("rho", dimDensity, viscosityProperties),
#else
        n_(KriegerDoughertyCoeffs_.lookup("n")),
        muPlasma_(KriegerDoughertyCoeffs_.lookup("muPlasma")),

        Hcrit_(KriegerDoughertyCoeffs_.lookup("Hcrit")),
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

    bool Foam::viscosityModels::KriegerDougherty::read
(
 const dictionary& viscosityProperties
 )
{
    viscosityModel::read(viscosityProperties);

    KriegerDoughertyCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    KriegerDoughertyCoeffs_.lookup("n") >> n_;
    KriegerDoughertyCoeffs_.lookup("muPlasma") >> muPlasma_;

    KriegerDoughertyCoeffs_.lookup("Hcrit") >> Hcrit_;

    viscosityProperties.lookup("rho") >> rho_;

    return true;
}


// ************************************************************************* //
