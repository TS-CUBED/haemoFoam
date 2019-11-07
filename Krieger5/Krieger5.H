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
    Foam::viscosityModels::Krieger5

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

SourceFiles
    Krieger5.C
\*---------------------------------------------------------------------------*/

#ifndef Krieger5_H
#define Krieger5_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class Krieger5 Declaration
\*---------------------------------------------------------------------------*/

class Krieger5
:
    public viscosityModel
{
    // Private data

        dictionary Krieger5Coeffs_;
        
        dimensionedScalar a_, b_, c_;
        dimensionedScalar beta_, lambda_, nuK_;
        dimensionedScalar muPlasma_;
		dimensionedScalar Hcrit_;
		dimensionedScalar rho_;
        volScalarField nu_;


    // Private Member Functions

        //- Calculate and return the laminar viscosity
        tmp<volScalarField> calcNu() const;


public:

    //- Runtime type information
    TypeName("Krieger5");


    // Constructors

        //- construct from components
        Krieger5
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Destructor

        virtual ~Krieger5()
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
