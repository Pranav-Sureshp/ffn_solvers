fvVectorMatrix UEqn
(
    fvm::ddt(U_)
  + fvm::div(phi_, U_)
  - fvc::div(phiB_, 2.0*DBU_*B_)
  - fvm::laplacian(nu_, U_)
  + fvc::grad(DBU_*magSqr(B_))
);

if (piso_.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p_));
}

// ************************************************************************* //
