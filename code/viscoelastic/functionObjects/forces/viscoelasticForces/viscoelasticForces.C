/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "viscoelasticForces.H"
#include "fvcGrad.H"
#include "viscoelasticModel.H"
#include "addToRunTimeSelectionTable.H"
#include "cartesianCS.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(viscoelasticForces, 0);
    addToRunTimeSelectionTable(functionObject, viscoelasticForces, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::word Foam::functionObjects::viscoelasticForces::fieldName
(
	const word& name
) const
{
    return this->name() + ":" + name;
}


void Foam::functionObjects::viscoelasticForces::createFiles()
{
    if (writeToFile() && !forceFilePtr_)
    {
        forceFilePtr_ = createFile("force");
        writeIntegratedHeader("Force", forceFilePtr_());
        momentFilePtr_ = createFile("moment");
        writeIntegratedHeader("Moment", momentFilePtr_());
    }
}


void Foam::functionObjects::viscoelasticForces::writeIntegratedHeader
(
    const word& header,
    Ostream& os
) const
{
    writeHeader(os, header);
    writeHeaderValue(os, "CofR", coordSys_.origin());
    writeHeader(os, "");
    writeCommented(os, "Time");
    writeTabbed(os, "(total_x total_y total_z)");
    writeTabbed(os, "(pressure_x pressure_y pressure_z)");
    writeTabbed(os, "(viscous_x viscous_y viscous_z)");

    os  << endl;
}


void Foam::functionObjects::viscoelasticForces::setCoordinateSystem
(
    const dictionary& dict,
    const word& e3Name,
    const word& e1Name
)
{
    coordSys_.clear();

    if (dict.readIfPresent<point>("CofR", coordSys_.origin()))
    {
        const vector e3 = e3Name == word::null ?
            vector(0, 0, 1) : dict.get<vector>(e3Name);
        const vector e1 = e1Name == word::null ?
            vector(1, 0, 0) : dict.get<vector>(e1Name);

        coordSys_ =
            coordSystem::cartesian(coordSys_.origin(), e3, e1);
    }
    else
    {
        // The 'coordinateSystem' sub-dictionary is optional,
        // but enforce use of a cartesian system if not found.

        if (dict.found(coordinateSystem::typeName_()))
        {
            // New() for access to indirect (global) coordinate system
            coordSys_ =
                coordinateSystem::New
                (
                    obr_,
                    dict,
                    coordinateSystem::typeName_()
                );
        }
        else
        {
            coordSys_ = coordSystem::cartesian(dict);
        }
    }

}


void Foam::functionObjects::viscoelasticForces::initialise()
{
    if (initialised_)
    {
        return;
    }

	if
	(
		!foundObject<volVectorField>(UName_)
	 || !foundObject<volScalarField>(pName_)

	)
	{
		FatalErrorInFunction
			<< "Could not find U: " << UName_ << " or p:" << pName_
			<< " in database"
			<< exit(FatalError);
	}

	forAll(force_, i)
	{
		force_[i].setSize(1, vector::zero);
		moment_[i].setSize(1, vector::zero);
	}

    initialised_ = true;
}


void Foam::functionObjects::viscoelasticForces::resetFields()
{
    force_[0] = Zero;
    force_[1] = Zero;

    moment_[0] = Zero;
    moment_[1] = Zero;

    if (writeFields_)
    {
        volVectorField& force =
            lookupObjectRef<volVectorField>(fieldName("force"));

        force == dimensionedVector(force.dimensions(), Zero);

        volVectorField& moment =
            lookupObjectRef<volVectorField>(fieldName("moment"));

        moment == dimensionedVector(moment.dimensions(), Zero);
    }
}


Foam::tmp<Foam::volSymmTensorField>
Foam::functionObjects::viscoelasticForces::tau() const
{
	if (foundObject<viscoelasticModel>("viscoelasticProperties"))
	{
		const viscoelasticModel& model =
			lookupObject<viscoelasticModel>("viscoelasticProperties");

		return model.tau();
	}

	FatalErrorInFunction
		<< "No valid model for viscous stress calculation"
		<< exit(FatalError);

	return volSymmTensorField::null();
}


Foam::scalar Foam::functionObjects::viscoelasticForces::rho() const
{
	if (foundObject<viscoelasticModel>("viscoelasticProperties"))
	{
		const viscoelasticModel& model =
			lookupObject<viscoelasticModel>("viscoelasticProperties");

		return model.rho().value();
	}

	FatalErrorInFunction
		<< "No valid model for density"
		<< exit(FatalError);

	return Zero;
}


void Foam::functionObjects::viscoelasticForces::addToFields
(
    const label patchi,
    const vectorField& Md,
    const vectorField& fN,
    const vectorField& fT
)
{
    if (!writeFields_)
    {
        return;
    }

    auto& force = lookupObjectRef<volVectorField>(fieldName("force"));
    vectorField& pf = force.boundaryFieldRef()[patchi];
    pf += fN + fT;

    auto& moment = lookupObjectRef<volVectorField>(fieldName("moment"));
    vectorField& pm = moment.boundaryFieldRef()[patchi];
    pm = Md^pf;
}


void Foam::functionObjects::viscoelasticForces::addToFields
(
    const labelList& cellIDs,
    const vectorField& Md,
    const vectorField& fN,
    const vectorField& fT
)
{
    if (!writeFields_)
    {
        return;
    }

    auto& force = lookupObjectRef<volVectorField>(fieldName("force"));
    auto& moment = lookupObjectRef<volVectorField>(fieldName("moment"));

    forAll(cellIDs, i)
    {
        label celli = cellIDs[i];
        force[celli] += fN[i] + fT[i];
        moment[celli] = Md[i]^force[celli];
    }
}


void Foam::functionObjects::viscoelasticForces::writeIntegratedForceMoment
(
    const string& descriptor,
    const vectorField& fm0,
    const vectorField& fm1,
    autoPtr<OFstream>& osPtr
) const
{
    vector pressure = sum(fm0);
    vector viscous = sum(fm1);
    vector total = pressure + viscous;

    Log << "    Sum of " << descriptor.c_str() << nl
        << "        Total    : " << total << nl
        << "        Pressure : " << pressure << nl
        << "        Viscous  : " << viscous << nl;

    if (writeToFile())
    {
        Ostream& os = osPtr();

        writeCurrentTime(os);

        os  << tab << total
            << tab << pressure
            << tab << viscous;

        os  << endl;
    }
}


void Foam::functionObjects::viscoelasticForces::writeForces()
{
    Log << type() << " " << name() << " write:" << nl;

    writeIntegratedForceMoment
    (
        "viscoelasticForces",
        coordSys_.localVector(force_[0]),
        coordSys_.localVector(force_[1]),
        forceFilePtr_
    );

    writeIntegratedForceMoment
    (
        "moments",
        coordSys_.localVector(moment_[0]),
        coordSys_.localVector(moment_[1]),
        momentFilePtr_
    );

    Log << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::viscoelasticForces::viscoelasticForces
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    bool readFields
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name),
    force_(2),
    moment_(2),
    forceFilePtr_(),
    momentFilePtr_(),
    patchSet_(),
    pName_("p"),
    UName_("U"),
    pRef_(0),
    coordSys_(),
    writeFields_(false),
    initialised_(false)
{
    if (readFields)
    {
        read(dict);
        setCoordinateSystem(dict);
        Log << endl;
    }
}


Foam::functionObjects::viscoelasticForces::viscoelasticForces
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    bool readFields
)
:
    fvMeshFunctionObject(name, obr, dict),
    writeFile(mesh_, name),
    force_(2),
    moment_(2),
    forceFilePtr_(),
    momentFilePtr_(),
    patchSet_(),
    pName_("p"),
    UName_("U"),
    pRef_(0),
    coordSys_(),
    writeFields_(false),
    initialised_(false)
{
    if (readFields)
    {
        read(dict);
        setCoordinateSystem(dict);
        Log << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::viscoelasticForces::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    initialised_ = false;

    Info<< type() << " " << name() << ":" << nl;

    patchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            dict.get<wordRes>("patches")
        );

	// Optional field name entries
	if (dict.readIfPresent<word>("p", pName_))
	{
		Info<< "    p: " << pName_ << endl;
	}
	if (dict.readIfPresent<word>("U", UName_))
	{
		Info<< "    U: " << UName_ << endl;
	}

	// Reference pressure, 0 by default
	if (dict.readIfPresent<scalar>("pRef", pRef_))
	{
		Info<< "    Reference pressure (pRef) set to " << pRef_ << endl;
	}

    writeFields_ = dict.getOrDefault("writeFields", false);

    if (writeFields_)
    {
        Info<< "    Fields will be written" << endl;

        volVectorField* forcePtr
        (
            new volVectorField
            (
                IOobject
                (
                    fieldName("force"),
                    time_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedVector(dimForce, Zero)
            )
        );

        mesh_.objectRegistry::store(forcePtr);

        volVectorField* momentPtr
        (
            new volVectorField
            (
                IOobject
                (
                    fieldName("moment"),
                    time_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedVector(dimForce*dimLength, Zero)
            )
        );

        mesh_.objectRegistry::store(momentPtr);
    }

    return true;
}


void Foam::functionObjects::viscoelasticForces::calcForcesMoment()
{
    initialise();

    resetFields();

	const volScalarField& p = lookupObject<volScalarField>(pName_);

	const surfaceVectorField::Boundary& Sfb = mesh_.Sf().boundaryField();

	tmp<volSymmTensorField> ttau = tau();
	const volSymmTensorField::Boundary& taub
		= ttau().boundaryField();

	// Scale pRef by density for incompressible simulations
	scalar pRef = pRef_/rho();

	for (const label patchi : patchSet_)
	{
		vectorField Md
		(
			mesh_.C().boundaryField()[patchi] - coordSys_.origin()
		);

		vectorField fN
		(
			rho()*Sfb[patchi]*(p.boundaryField()[patchi] - pRef)
		);

		// reverse the normal to get the force on the wall
		vectorField fT(-Sfb[patchi] & taub[patchi]);

		addToFields(patchi, Md, fN, fT);

		force_[0][0] += sum(fN);
		force_[1][0] += sum(fT);
		moment_[0][0] += sum(Md^fN);
		moment_[1][0] += sum(Md^fT);
	}

    Pstream::listCombineGather(force_, plusEqOp<vectorField>());
    Pstream::listCombineGather(moment_, plusEqOp<vectorField>());
    Pstream::listCombineScatter(force_);
    Pstream::listCombineScatter(moment_);
}


Foam::vector Foam::functionObjects::viscoelasticForces::forceEff() const
{
    return sum(force_[0]) + sum(force_[1]);
}


Foam::vector Foam::functionObjects::viscoelasticForces::momentEff() const
{
    return sum(moment_[0]) + sum(moment_[1]);
}


bool Foam::functionObjects::viscoelasticForces::execute()
{
    calcForcesMoment();

    if (Pstream::master())
    {
        createFiles();

        writeForces();

        Log << endl;
    }

    // Write state/results information
    setResult("normalForce", sum(force_[0]));
    setResult("tangentialForce", sum(force_[1]));

    setResult("normalMoment", sum(moment_[0]));
    setResult("tangentialMoment", sum(moment_[1]));

    return true;
}


bool Foam::functionObjects::viscoelasticForces::write()
{
    if (writeFields_)
    {
        lookupObject<volVectorField>(fieldName("force")).write();
        lookupObject<volVectorField>(fieldName("moment")).write();
    }

    return true;
}


// ************************************************************************* //
