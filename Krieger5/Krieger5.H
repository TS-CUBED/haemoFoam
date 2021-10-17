/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1812                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*\
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
