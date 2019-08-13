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
    
    a_(Krieger5Coeffs_.lookup("a")),
    b_(Krieger5Coeffs_.lookup("b")),
    c_(Krieger5Coeffs_.lookup("c")),
    
    beta_(Krieger5Coeffs_.lookup("beta")),
    lambda_(Krieger5Coeffs_.lookup("lambda")),
    
    nuK_(Krieger5Coeffs_.lookup("nuK")),
    
    muPlasma_(Krieger5Coeffs_.lookup("muPlasma")),
    
    Hcrit_(Krieger5Coeffs_.lookup("Hcrit")),
    
    rho_(viscosityProperties.lookup("rho")),
    
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
