/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1812                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*\
Description
    An incompressible 5-Parameter modified Krieger non-Newtonian viscosity
    model following

    Hund, S.; Kamenenva, M.; Antaki, J.
    A Quasi-Mechanistic Mathematical Representation for Blood Viscosity
    Fluids 2017, 2, 10; doi:10.3390/fluids2010010


Created by
    Dr Torsten Schenkel
    Department Engineering and Mathematics
    Material and Engineering Research Institute MERI
    Sheffield Hallam University
    December 2018
\*---------------------------------------------------------------------------*/

#include "Krieger5.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace viscosityModels
    {
        defineTypeNameAndDebug(Krieger5, 0);
        addToRunTimeSelectionTable
            (
             viscosityModel,
             Krieger5,
             dictionary
            );
    }
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::Krieger5::calcNu() const
{
    const volScalarField& H= U_.mesh().lookupObject<volScalarField>("H");

    return muPlasma_/rho_ *
        pow((1-H/Hcrit_),-1*((a_+b_*exp(-1*c_*H)+beta_*pow((1+pow(lambda_*strainRate(),2)),(-1*nuK_)))));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    Foam::viscosityModels::Krieger5::Krieger5
(
 const word& name,
 const dictionary& viscosityProperties,
 const volVectorField& U,
 const surfaceScalarField& phi
 )
    :
        viscosityModel(name, viscosityProperties, U, phi),
        Krieger5Coeffs_(viscosityProperties.subDict(typeName + "Coeffs")),

        
#ifdef OPENFOAMESIORFOUNDATION
        a_("a", dimless, Krieger5Coeffs_),
        b_("b", dimless, Krieger5Coeffs_),
        c_("c", dimless, Krieger5Coeffs_),

        beta_("beta", dimless, Krieger5Coeffs_),
        lambda_("lambda", dimTime, Krieger5Coeffs_),

        nuK_("nuK", dimless, Krieger5Coeffs_),

        muPlasma_("muPlasma", dimViscosity*dimDensity, Krieger5Coeffs_),
        
        Hcrit_("Hcrit", dimless, Krieger5Coeffs_),

        rho_("rho", dimDensity, viscosityProperties),
#else
        a_(Krieger5Coeffs_.lookup("a")),
        b_(Krieger5Coeffs_.lookup("b")),
        c_(Krieger5Coeffs_.lookup("c")),

        beta_(Krieger5Coeffs_.lookup("beta")),
        lambda_(Krieger5Coeffs_.lookup("lambda")),

        nuK_(Krieger5Coeffs_.lookup("nuK")),

        muPlasma_(Krieger5Coeffs_.lookup("muPlasma")),

        Hcrit_(Krieger5Coeffs_.lookup("Hcrit")),

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

    bool Foam::viscosityModels::Krieger5::read
(
 const dictionary& viscosityProperties
 )
{
    viscosityModel::read(viscosityProperties);

    Krieger5Coeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    Krieger5Coeffs_.lookup("a") >> a_;
    Krieger5Coeffs_.lookup("b") >> b_;
    Krieger5Coeffs_.lookup("c") >> c_;

    Krieger5Coeffs_.lookup("beta") >> beta_;
    Krieger5Coeffs_.lookup("lambda") >> lambda_;

    Krieger5Coeffs_.lookup("nuK") >> nuK_;

    Krieger5Coeffs_.lookup("muPlasma") >> muPlasma_;


    Krieger5Coeffs_.lookup("Hcrit") >> Hcrit_;

    viscosityProperties.lookup("rho") >> rho_;

    return true;
}


// ************************************************************************* //
