/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

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

        n_(CarreauCoeffs_.lookup("n")),

        lambda_(CarreauCoeffs_.lookup("lambda")),

        mu0_(CarreauCoeffs_.lookup("mu0")),
        muInf_(CarreauCoeffs_.lookup("muInf")),

        rho_(CarreauCoeffs_.lookup("rho")),

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
