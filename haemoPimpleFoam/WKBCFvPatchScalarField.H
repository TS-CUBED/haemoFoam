/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    WKBCFvPatchVectorField

Description
    Windkessel pressure outflow boundary condition, 3 element. 

SourceFiles
    WKBCFvPatchVectorField.C

Author
    Andris Piebalgs, Imperial College London, 2017  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef WKBCFvPatchScalarField_H
#define WKBCFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class pulsatileVelocityFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class WKBCFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

    // Name of the outlet
    scalar index_;



public:

    //- Runtime type information
    TypeName("WKBC");


    // Constructors

        //- Construct from patch and internal field
        WKBCFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        WKBCFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given WKBCFvPatchScalarField
        //  onto a new patch
        WKBCFvPatchScalarField
        (
            const WKBCFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new WKBCFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        WKBCFvPatchScalarField
        (
            const WKBCFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new WKBCFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

     // Return outlet index

  scalar& index()
  {
    return index_;
  }


        //- Update coefficients
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
