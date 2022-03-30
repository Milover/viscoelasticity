/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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

#include "pBar.H"
#include "volFieldsFwd.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(pBar, 0);
    addToRunTimeSelectionTable(functionObject, pBar, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::pBar::calc()
{
    if (foundObject<volScalarField>(fieldName_))
    {
		dimensionedScalar rho("rho", dimDensity, Zero);

        if (mesh_.foundObject<viscoelasticModel>("viscoelasticModel"))
        {
            const dictionary& model =
                mesh_.lookupObject<viscoelasticModel>("viscoelasticModel");

			rho = model.rho();
        }
        else
        {
            FatalErrorInFunction
                << "Unable to determine the density"
                << exit(FatalError);
        }

        const volScalarField& p =
            mesh_.lookupObject<volScalarField>(fieldName_);
        const volVectorField& U =
            mesh_.lookupObject<volVectorField>("U");

		const dimensionedScalar small("small", p.dimensions(), SMALL);

        return store(resultName_, p / (rho*magSqr(U) + small));
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::pBar::pBar
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "p", "pBar")
{
    read(dict);
}


// ************************************************************************* //
