/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "targetReFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "viscoelasticModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::targetReFvPatchVectorField::targetReFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    targetRe_(),
    Re_(0.0),
    patchName_(),
    length_(0.0),
    weightedAverage_(false)
{}


Foam::targetReFvPatchVectorField::targetReFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict, false),
    Re_(0.0),
	patchName_(dict.get<word>("patch")),
    length_(dict.get<scalar>("length")),
    weightedAverage_(dict.getOrDefault<Switch>("weighted", false))
{
	targetRe_ = Function1<scalar>::New("target", dict, &db());

	// Sanity check
	if (length_ < VSMALL)
	{
		FatalErrorInFunction
			<< "Characteristic length scale = " << length_
			<< " cannot be negative or zero."
			<< abort(FatalError);
	}

    // Value field require if mass based
    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        evaluate(Pstream::commsTypes::blocking);
    }
}


Foam::targetReFvPatchVectorField::targetReFvPatchVectorField
(
    const targetReFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    targetRe_(ptf.targetRe_.clone()),
    Re_(ptf.Re_),
    patchName_(ptf.patchName_),
    length_(ptf.length_),
    weightedAverage_(ptf.weightedAverage_)
{}


Foam::targetReFvPatchVectorField::targetReFvPatchVectorField
(
    const targetReFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    targetRe_(ptf.targetRe_.clone()),
    Re_(ptf.Re_),
    patchName_(ptf.patchName_),
    length_(ptf.length_),
    weightedAverage_(ptf.weightedAverage_)
{}


Foam::targetReFvPatchVectorField::targetReFvPatchVectorField
(
    const targetReFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    targetRe_(ptf.targetRe_.clone()),
    Re_(ptf.Re_),
    patchName_(ptf.patchName_),
    length_(ptf.length_),
    weightedAverage_(ptf.weightedAverage_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::targetReFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

	const scalar t = db().time().timeOutputValue();

	const fvMesh& mesh = static_cast<const fvMesh&>(db());
	const label patchId = mesh.boundary().findPatchID(patchName_);

	// update patch values
	volVectorField& U
	(
		const_cast<volVectorField&>
		(
			db().lookupObject<volVectorField>("U")
		)
	);
	U.boundaryFieldRef()[patchId].updateCoeffs();

	if (patchId < 0)
	{
		FatalErrorInFunction
			<< "Unable to find patch " << patchName_
			<< exit(FatalError);
	}

	if (!db().foundObject<viscoelasticModel>("viscoelasticProperties"))
	{
		FatalErrorInFunction
			<< "No valid model for Reynolds number computation"
			<< exit(FatalError);
	}
	const viscoelasticModel& model =
		db().lookupObject<viscoelasticModel>("viscoelasticProperties");

	const scalar rho = model.rho().value();
	const scalarField mu = model.eta()().boundaryField()[patchId];
	const vectorField& Up = U.boundaryField()[patchId];
	const vectorField n = mesh.boundary()[patchId].nf();
	const scalarField& magSf = mesh.boundary()[patchId].magSf();
	const scalarField nUp = n & Up;

	// NOTE: another option would be to compare with the weighted average of Re
	const scalar phi = gSum(nUp*magSf);
	const scalar Af = gSum(magSf);

	const scalar avgUn = phi/Af;
	scalar avgMu;
	if (weightedAverage_)
	{
		avgMu = gSum(nUp*magSf*mu)/phi;
	}
	else
	{
		avgMu = gSum(magSf*mu)/Af;
	}

	Re_ = rho*avgUn*length_/avgMu;

	if (mag(targetRe_->value(t) - Re_) > VSMALL)
	{
		const vectorField np = patch().nf();
		operator==(-np*(targetRe_->value(t)*avgMu/(rho*length_)));
	}

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::targetReFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    targetRe_->writeData(os);
	os.writeEntry("Re", Re_);
	os.writeEntry("patch", patchName_);

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       targetReFvPatchVectorField
   );
}


// ************************************************************************* //
