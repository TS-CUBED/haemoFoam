/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1812                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/


#include "Carreau.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace viscosityModels
    {
        defineTypeNameAndDebug(Carreau, 0);
        addToRunTimeSelectionTable
            (
             viscosityModel,
             Carreau,
             dictionary
            );
    }
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::Carreau::calcNu() const
{
    //    const volScalarField& H= U_.mesh().lookupObject<volScalarField>("H");

    return muInf_/rho_ +
        (mu0_ - muInf_)/rho_ * pow((1.0 + pow(strainRate()*lambda_,2)),((n_ - 1)/2));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    Foam::viscosityModels::Carreau::Carreau
(
 const word& name,
 const dictionary& viscosityProperties,
 const volVectorField& U,
 const surfaceScalarField& phi
 )
    :
        viscosityModel(name, viscosityProperties, U, phi),
        CarreauCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),

#ifdef OPENFOAMESIORFOUNDATION
        n_("n", dimless, CarreauCoeffs_),

        lambda_("lambda", dimTime, CarreauCoeffs_),

        mu0_("mu0", dimViscosity*dimDensity, CarreauCoeffs_),
        muInf_("muInf", dimViscosity*dimDensity, CarreauCoeffs_),
        
        rho_("rho", dimDensity, viscosityProperties),
#else
        n_(CarreauCoeffs_.lookup("n")),

        lambda_(CarreauCoeffs_.lookup("lambda")),

        mu0_(CarreauCoeffs_.lookup("mu0")),
        muInf_(CarreauCoeffs_.lookup("muInf")),

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

    bool Foam::viscosityModels::Carreau::read
(
 const dictionary& viscosityProperties
 )
{
    viscosityModel::read(viscosityProperties);

    CarreauCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    CarreauCoeffs_.lookup("n") >> n_;

    CarreauCoeffs_.lookup("lambda") >> lambda_;

    CarreauCoeffs_.lookup("mu0") >> mu0_;
    CarreauCoeffs_.lookup("muInf") >> muInf_;

    CarreauCoeffs_.lookup("rho") >> rho_;

    return true;
}


// ************************************************************************* //
