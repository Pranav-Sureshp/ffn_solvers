fvVectorMatrix BEqn
(
    fvm::ddt(B_)
  + fvm::div(phi_, B_)
  - fvc::div(phiB_, U_)
  - fvm::laplacian(DB_, B_)
);

BEqn.solve();

phiB_ = fvc::flux(B_);

// ************************************************************************* //
