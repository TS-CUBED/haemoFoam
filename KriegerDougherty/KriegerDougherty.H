/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1812                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*\
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
