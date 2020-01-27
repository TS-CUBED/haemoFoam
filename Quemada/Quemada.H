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
    Foam::viscosityModels::Quemada

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

SourceFiles
    Quemada.C

\*---------------------------------------------------------------------------*/

#ifndef Quemada_H
#define Quemada_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace viscosityModels
    {

        /*---------------------------------------------------------------------------*\
          Class Quemada Declaration
          \*---------------------------------------------------------------------------*/

        class Quemada
            :
                public viscosityModel
        {
            // Private data

            dictionary QuemadaCoeffs_;

            dimensionedScalar a0_;
            dimensionedScalar a1_;
            dimensionedScalar b0_;
            dimensionedScalar b1_;
            dimensionedScalar b2_;
            dimensionedScalar b3_;
            dimensionedScalar c0_;
            dimensionedScalar c1_;
            dimensionedScalar c2_;
            dimensionedScalar c3_;

            dimensionedScalar gammaC0_;
            dimensionedScalar muPlasma_;
            dimensionedScalar rho_;

            volScalarField nu_;


            // Private Member Functions

            //- Calculate and return the laminar viscosity
            tmp<volScalarField> calcNu() const;


            public:

            //- Runtime type information
            TypeName("Quemada");


            // Constructors

            //- construct from components
            Quemada
                (
                 const word& name,
                 const dictionary& viscosityProperties,
                 const volVectorField& U,
                 const surfaceScalarField& phi
                );


            // Destructor

            virtual ~Quemada()
            {}


            // Member Functions

            //- Return the laminar viscosity
            virtual const volScalarField& nu() const
            {
                return nu_;
            }

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
