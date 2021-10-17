/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1812                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*\
Class
    Foam::viscosityModels::Yeleswarapu

Description
    An incompressible non-Newtonian viscosity
    model following

    K. K. Yeleswarapu, M. V. Kameneva, K. R. Rajagopal, and J. F. Antaki,
    ‘The flow of blood in tubes:theory and experiment’,
    Mechanics Research Communications, vol. 25, no. 3, pp. 257–262, May 1998.

    W.-T. Wu, F. Yang, J. F. Antaki, N. Aubry, and M. Massoudi,
    ‘Study of blood flow in several benchmark micro-channels using a two-fluid approach’,
    Int J Eng Sci, vol. 95, pp. 49–59, Oct. 2015.

Author and Copyright
    Dr Torsten Schenkel
    Department Engineering and Mathematics
    Material and Engineering Research Institute MERI
    Sheffield Hallam University
    April 2021
    All Rights Reserved

SourceFiles
    Yeleswarapu.C
\*---------------------------------------------------------------------------*/

#ifndef Yeleswarapu_H
#define Yeleswarapu_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace viscosityModels
    {

        /*---------------------------------------------------------------------------*\
          Class Yeleswarapu Declaration
          \*---------------------------------------------------------------------------*/

        class Yeleswarapu
            :
                public viscosityModel
        {
            // Private data

            dictionary YeleswarapuCoeffs_;

            dimensionedScalar a1_;
            dimensionedScalar a2_;
            dimensionedScalar a3_;
            dimensionedScalar b1_;
            dimensionedScalar b2_;
            dimensionedScalar b3_;
            
            dimensionedScalar k_;

            dimensionedScalar muPlasma_;
            dimensionedScalar rho_;

            volScalarField nu_;


            // Private Member Functions

            //- Calculate and return the laminar viscosity
            tmp<volScalarField> calcNu() const;


            public:

            //- Runtime type information
            TypeName("Yeleswarapu");


            // Constructors

            //- construct from components
            Yeleswarapu
                (
                 const word& name,
                 const dictionary& viscosityProperties,
                 const volVectorField& U,
                 const surfaceScalarField& phi
                );


            // Destructor

            virtual ~Yeleswarapu()
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
