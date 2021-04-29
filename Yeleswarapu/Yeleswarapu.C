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
