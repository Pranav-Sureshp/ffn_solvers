{
    volScalarField rAB(1.0/BEqn.A());
    surfaceScalarField rABf("rABf", fvc::interpolate(rAB));

    while (bpiso_.correctNonOrthogonal())
    {
        fvScalarMatrix pBEqn
        (
            fvm::laplacian(rABf, pB_) == fvc::div(phiB_)
        );

        pBEqn.solve(mesh_.solver(pB_.select(bpiso_.finalInnerIter())));

        if (bpiso_.finalNonOrthogonalIter())
        {
            phiB_ -= pBEqn.flux();
        }
    }
}

// ************************************************************************* //
