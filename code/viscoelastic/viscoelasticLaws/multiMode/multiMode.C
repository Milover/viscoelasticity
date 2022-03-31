/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

\*---------------------------------------------------------------------------*/

#include "multiMode.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiMode, 0);
    addToRunTimeSelectionTable(viscoelasticLaw, multiMode, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiMode::multiMode
(
    const word& name,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const dictionary& dict
)
:
    viscoelasticLaw(name, U, phi),
    tau_
    (
        IOobject
        (
            "tau" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedSymmTensor("zero", dimPressure, Zero)
    ),
    models_()
{
    PtrList<entry> modelEntries(dict.lookup("models"));
    models_.setSize(modelEntries.size());

    forAll (models_, modelI)
    {
        models_.set
        (
            modelI,
            viscoelasticLaw::New
            (
                modelEntries[modelI].keyword(),
                U,
                phi,
                modelEntries[modelI].dict()
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvVectorMatrix> Foam::multiMode::divTau
(
	const volVectorField& U
) const
{
	tmp<fvVectorMatrix> tdivMatrix
	(
		new fvVectorMatrix(U, dimPressure/dimDensity)
	);
	fvVectorMatrix& divMatrix = tdivMatrix.ref();

    for (const auto& m : models_)
    {
        divMatrix += m.divTau(U);
    }

    return tdivMatrix;
}


Foam::dimensionedScalar Foam::multiMode::rho() const
{
	scalar rho = Zero;
	label n = Zero;

	for (const auto& m : models_)
    {
		++n;
        rho += (m.rho().value() - rho) / static_cast<scalar>(n);
    }

	return dimensionedScalar("rho", dimDensity, rho);
}


Foam::tmp<Foam::volSymmTensorField> Foam::multiMode::tau() const
{
    tau_ = symmTensor::zero;

	for (const auto& m : models_)
    {
        tau_ += m.tau();
    }

    return tau_;
}


void Foam::multiMode::correct()
{
	for (auto& m : models_)
    {
		DebugInFunction
			<< "Correcting " << m.name() << nl;

        m.correct();
    }

    tau();
}


// ************************************************************************* //
