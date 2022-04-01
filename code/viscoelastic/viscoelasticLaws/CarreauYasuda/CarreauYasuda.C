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

#include "CarreauYasuda.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(CarreauYasuda, 0);
    addToRunTimeSelectionTable(viscoelasticLaw, CarreauYasuda, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CarreauYasuda::CarreauYasuda
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
    eta_
    (
        IOobject
        (
            "eta" + name,
            U.time().timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        strainRate()*dimensionedScalar("zero", dimVelocity, Zero) //Just to ensure dimensions and BCs
    ),
    rho_(dict.get<dimensionedScalar>("rho")),
    eta0_(dict.get<dimensionedScalar>("eta0")),
    etaInf_(dict.get<dimensionedScalar>("etaInf")),
    k_(dict.get<dimensionedScalar>("k")),
    a_(dict.get<dimensionedScalar>("a")),
    n_(dict.get<dimensionedScalar>("n"))
{
	correct();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvVectorMatrix> Foam::CarreauYasuda::divTau
(
  const volVectorField& U
) const
{
	return
	(
		fvm::laplacian(eta_/rho_, U, "laplacian(eta,U)")
	  + (fvc::grad(U) & fvc::grad(eta_/rho_))
	);
}


void Foam::CarreauYasuda::correct()
{
	eta_ = etaInf_
		 + (eta0_ - etaInf_)
		 * pow
		   (
		   		scalar(1.0) + pow(k_*strainRate(), a_),
				(n_ - scalar(1.0))/a_
		   );

	tau_ = eta_*twoSymm(U());
}


// ************************************************************************* //
