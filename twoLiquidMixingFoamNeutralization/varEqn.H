// Variance equation
fvScalarMatrix varEqn(
    fvm::ddt(var) 
  + fvm::div(phi, var) 
  - fvm::laplacian
    (
        volScalarField("Dab", Dab + alphatab * turbulence->nut()),
        var
    )
  - 2 * volScalarField("DabT", alphatab * turbulence->nut()) * (fvc::grad(alpha1) & fvc::grad(alpha1))
  + fvc::SuSp(Cphi_.CphiValue() * ( turbulence->epsilon() / turbulence->k() ) / 2.0, var)
);

varEqn.relax();
varEqn.solve();