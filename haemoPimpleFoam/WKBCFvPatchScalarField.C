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

\*---------------------------------------------------------------------------*/

#include "WKBCFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "Time.H"
#include "mathematicalConstants.H"
#include "scalarIOList.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

WKBCFvPatchScalarField::WKBCFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    index_(0)

{}


WKBCFvPatchScalarField::WKBCFvPatchScalarField
(
    const WKBCFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    index_(ptf.index_)

{}


WKBCFvPatchScalarField::WKBCFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    index_(readScalar(dict.lookup("index")))
{

    evaluate();
}


WKBCFvPatchScalarField::WKBCFvPatchScalarField
(
    const WKBCFvPatchScalarField& fcvpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(fcvpvf, iF),
    index_(fcvpvf.index_)

{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void WKBCFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    
    /*Creating a patch field of same size as the boundary field*/
    const fvPatch& p = this->patch();
    scalarField report(p.size());


    /* Accessing the variables stored in mesh */
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const scalarIOList& store = mesh.lookupObject<scalarIOList>("store");

    const scalar current_pressure = store[index_];

	//Info << "store in WKBC " << store << nl << endl;
	//Info << "current_pressure" << current_pressure << nl << endl;

    /*Applying the pressure to each face on the outlet*/
    forAll(report,it)
    {
		report[it] = current_pressure/1060;
		
    }

	//Info << "BC P " << average(report)*1060 << nl << endl;

    /*Assigning the operator to the new patch field*/
    scalarField::operator= 
    (
     report
    );
    
}


// Write
void WKBCFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("index") 
        << index_ << token::END_STATEMENT << nl;
#ifdef OPENFOAMESI
    writeEntry("value", os);
#else
    writeEntry(os, "value", *this);
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, WKBCFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
