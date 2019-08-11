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
    
    n_(KriegerDoughertyCoeffs_.lookup("n")),
    muPlasma_(KriegerDoughertyCoeffs_.lookup("muPlasma")),
    rho_(KriegerDoughertyCoeffs_.lookup("rho")),
    
    Hcrit_(KriegerDoughertyCoeffs_.lookup("Hcrit")),
    
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
    KriegerDoughertyCoeffs_.lookup("rho") >> rho_;
    
    KriegerDoughertyCoeffs_.lookup("Hcrit") >> Hcrit_;
    
    return true;
}


// ************************************************************************* //
