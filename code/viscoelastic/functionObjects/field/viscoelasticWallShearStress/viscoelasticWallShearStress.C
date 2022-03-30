/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

#include "viscoelasticWallShearStress.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(viscoelasticWallShearStress, 0);
    addToRunTimeSelectionTable(functionObject, viscoelasticWallShearStress, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::viscoelasticWallShearStress::writeFileHeader
(
	Ostream& os
) const
{
    writeHeader(os, "Viscoelastic wall shear stress");
    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "min");
    writeTabbed(os, "max");
    os  << endl;
}


void Foam::functionObjects::viscoelasticWallShearStress::calcShearStress
(
    const volSymmTensorField& tau,
    volVectorField& shearStress
)
{
    shearStress.dimensions().reset(Reff.dimensions());

    for (const label patchi : patchSet_)
    {
        vectorField& ssp = shearStress.boundaryFieldRef()[patchi];
        const vectorField& Sfp = mesh_.Sf().boundaryField()[patchi];
        const scalarField& magSfp = mesh_.magSf().boundaryField()[patchi];
        const symmTensorField& Reffp = Reff.boundaryField()[patchi];

        ssp = (-Sfp/magSfp) & tau;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::viscoelasticWallShearStress::viscoelasticWallShearStress
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    patchSet_()
{
    read(dict);

    writeFileHeader(file());

    volVectorField* viscoelasticWallShearStressPtr
    (
        new volVectorField
        (
            IOobject
            (
                typeName,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector(sqr(dimLength)/sqr(dimTime), Zero)
        )
    );

    mesh_.objectRegistry::store(viscoelasticWallShearStressPtr);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::viscoelasticWallShearStress::read
(
	const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    patchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            dict.getOrDefault<wordRes>("patches", wordRes())
        );

    Info<< type() << " " << name() << ":" << nl;

    if (patchSet_.empty())
    {
        forAll(pbm, patchi)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                patchSet_.insert(patchi);
            }
        }

        Info<< "    processing all wall patches" << nl << endl;
    }
    else
    {
        Info<< "    processing wall patches: " << nl;
        labelHashSet filteredPatchSet;
        for (const label patchi : patchSet_)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                filteredPatchSet.insert(patchi);
                Info<< "        " << pbm[patchi].name() << endl;
            }
            else
            {
                WarningInFunction
                    << "Requested wall shear stress on non-wall boundary "
                    << "type patch: " << pbm[patchi].name() << endl;
            }
        }

        Info<< endl;

        patchSet_ = filteredPatchSet;
    }

    return true;
}


bool Foam::functionObjects::viscoelasticWallShearStress::execute()
{
    volVectorField& viscoelasticWallShearStress =
        mesh_.lookupObjectRef<volVectorField>(type());

	if (mesh_.foundObject<volSymmTensorField>("tau"))
	{
		calcShearStress
		(
			mesh_.lookupObject<volSymmTensorField>("tau"),
			viscoelasticWallShearStress
		);
		return true;
	}

    FatalErrorInFunction
        << "Unable to find turbulence model in the "
        << "database" << exit(FatalError);

    return false;
}


bool Foam::functionObjects::viscoelasticWallShearStress::write()
{
    const volVectorField& viscoelasticWallShearStress =
        obr_.lookupObject<volVectorField>(type());

    Log << type() << " " << name() << " write:" << nl
        << "    writing field " << viscoelasticWallShearStress.name() << endl;

    viscoelasticWallShearStress.write();

    const fvPatchList& patches = mesh_.boundary();

    for (const label patchi : patchSet_)
    {
        const fvPatch& pp = patches[patchi];

        const vectorField& ssp = viscoelasticWallShearStress.boundaryField()[patchi];

        vector minSsp = gMin(ssp);
        vector maxSsp = gMax(ssp);

        if (Pstream::master())
        {
            writeCurrentTime(file());

            file()
                << token::TAB << pp.name()
                << token::TAB << minSsp
                << token::TAB << maxSsp
                << endl;
        }

        Log << "    min/max(" << pp.name() << ") = "
            << minSsp << ", " << maxSsp << endl;
    }

    return true;
}


// ************************************************************************* //
