Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

autoPtr<viscoelasticModel> visco;

autoPtr<singlePhaseTransportModel> laminarTransport;
autoPtr<incompressible::turbulenceModel> turbulence;

bool viscoelastic;
{
	Foam::fileOperations::uncollatedFileOperation op(false);
	viscoelastic = op.exists(runTime.constantPath() / "viscoelasticProperties");
}

// check solver mode
if (viscoelastic)
{
	visco.reset
	(
		new viscoelasticModel(U, phi)
	);
}
else
{
	laminarTransport.reset
	(
		new singlePhaseTransportModel(U, phi)
	);
	turbulence.reset
	(
		incompressible::turbulenceModel::New(U, phi, laminarTransport.ref())
	);
}

#include "createMRF.H"
#include "createFvOptions.H"
