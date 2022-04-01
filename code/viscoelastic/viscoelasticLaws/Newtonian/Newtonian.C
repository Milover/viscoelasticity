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

#include "Newtonian.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(Newtonian, 0);
	addToRunTimeSelectionTable(viscoelasticLaw, Newtonian, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Newtonian::Newtonian
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
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedSymmTensor("0", dimPressure, Zero),
        extrapolatedCalculatedFvPatchField<symmTensor>::typeName
    ),
    rho_(dict.get<dimensionedScalar>("rho")),
    eta_(dict.get<dimensionedScalar>("eta"))
{
	correct();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
 
Foam::tmp<Foam::fvVectorMatrix> Foam::Newtonian::divTau
(
  const volVectorField& U
) const
{
	return
	(
		fvm::laplacian( eta_/rho_, U, "laplacian(eta,U)")
	);
}

void Foam::Newtonian::correct()
{
	tau_ = eta_*twoSymm(U());
}


// ************************************************************************* //
