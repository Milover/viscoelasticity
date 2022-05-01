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

#include "UBar.H"
#include "volFields.H"
#include "volFieldsFwd.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(UBar, 0);
    addToRunTimeSelectionTable(functionObject, UBar, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::UBar::calc()
{
    if (mesh_.foundObject<volVectorField>(fieldName_))
    {
        const volVectorField& U =
			mesh_.lookupObject<volVectorField>(fieldName_);

        return store
		(
			resultName_,
			U / dimensionedScalar("URef_", dimVelocity, URef_)
		);
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::UBar::UBar
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "U", "UBar")
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::UBar::read(const dictionary& dict)
{
	if (fieldExpression::read(dict))
	{
		dict.readEntry("URef", URef_);

		if (URef_ < VSMALL)
		{
			FatalErrorInFunction
				<< "Reference velocity = " << URef_
				<< " cannot be negative or zero."
				<< abort(FatalError);
		}

		return true;
	}

	return false;
}


// ************************************************************************* //
