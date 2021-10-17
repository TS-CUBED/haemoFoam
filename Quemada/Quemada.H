/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1812                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*\
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
