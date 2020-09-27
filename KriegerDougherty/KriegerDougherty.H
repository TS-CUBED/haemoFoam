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

Class
    Foam::viscosityModels::KriegerDougherty

Description
    An incompressible Krieger-Dougherty non-Newtonian viscosity model.

Author and Copyright
    Dr Torsten Schenkel
    Department Engineering and Mathematics
    Material and Engineering Research Institute MERI
    Sheffield Hallam University
    February 2019
    All Rights Reserved

SourceFiles
    KriegerDougherty.C

\*---------------------------------------------------------------------------*/

#ifndef KriegerDougherty_H
#define KriegerDougherty_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace viscosityModels
    {

        /*---------------------------------------------------------------------------*\
          Class KriegerDougherty Declaration
          \*---------------------------------------------------------------------------*/

        class KriegerDougherty
            :
                public viscosityModel
        {
            // Private data

            dictionary KriegerDoughertyCoeffs_;

            dimensionedScalar n_;
            dimensionedScalar muPlasma_;
            dimensionedScalar Hcrit_;
            dimensionedScalar rho_;
            volScalarField nu_;


            // Private Member Functions

            //- Calculate and return the laminar viscosity
            tmp<volScalarField> calcNu() const;


            public:

            //- Runtime type information
            TypeName("KriegerDougherty");


            // Constructors

            //- construct from components
            KriegerDougherty
                (
                 const word& name,
                 const dictionary& viscosityProperties,
                 const volVectorField& U,
                 const surfaceScalarField& phi
                );


            // Destructor

            virtual ~KriegerDougherty()
            {}


            // Member Functions
#ifdef OPENFOAMESIORFOUNDATION
            //- Return the laminar viscosity
            virtual tmp<volScalarField> nu() const
            {
                return nu_;
            }
            
            //- Return the laminar viscosity for patch
            virtual tmp<scalarField> nu(const label patchi) const
            {
                return nu_.boundaryField()[patchi];
            }
#else
            //- Return the laminar viscosity
            virtual const volScalarField& nu() const
            {
                return nu_;
            }
#endif

            //- Correct the laminar viscosity
            virtual void correct()
            {
                nu_ = calcNu();
            }

            //- Read transportProperties dictionary
            virtual bool read(const dictionary& viscosityProperties);
        };


        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
