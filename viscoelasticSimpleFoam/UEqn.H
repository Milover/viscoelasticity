    // Momentum predictor

    MRF.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(U)
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

	// XXX: does putting it here affect fvOptions?
	// viscoelastic mode
	if (viscoelastic)
	{
		UEqn -= visco->divTau(U);
	}
	// simpleFoam mode
	else
	{
		UEqn += turbulence->divDevReff(U);
	}

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (simple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));

        fvOptions.correct(U);
    }

    tUEqn.clear();
