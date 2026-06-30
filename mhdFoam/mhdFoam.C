/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2011-2018 OpenFOAM Foundation (original mhdFoam)           |
|                                                                             |
|    Built on OpenFOAM v2512                                                  |
-------------------------------------------------------------------------------
License
    This file is part of GeN-Foam.

    GeN-Foam is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    GeN-Foam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    This offering is not approved or endorsed by the OpenFOAM Foundation nor
    OpenCFD Limited, producer and distributor of the OpenFOAM(R) software via
    www.openfoam.com, and owner of the OPENFOAM(R) and OpenCFD(R) trademarks.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "mhdFoam.H"
#include "addToRunTimeSelectionTable.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(mhdFoam, 0);
    addToRunTimeSelectionTable
    (
        solver,
        mhdFoam,
        dynamicFvMesh
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::mhdFoam::mhdFoam
(
    dynamicFvMesh& mesh
)
:
    solver(mesh),
    mesh_(mesh),
    piso_(mesh_),
    bpiso_(mesh_, "BPISO"),
    transportProperties_
    (
        IOobject
        (
            "transportProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    rho_
    (
        "rho",
        dimDensity,
        transportProperties_
    ),
    nu_
    (
        "nu",
        dimViscosity,
        transportProperties_
    ),
    mu_
    (
        "mu",
        dimensionSet(1, 1, -2, 0, 0, -2, 0),
        transportProperties_
    ),
    sigma_
    (
        "sigma",
        dimensionSet(-1, -3, 3, 0, 0, 2, 0),
        transportProperties_
    ),
    DB_
    (
        "DB",
        1.0/(mu_*sigma_)
    ),
    DBU_
    (
        "DBU",
        1.0/(2.0*mu_*rho_)
    ),
    p_
    (
        IOobject
        (
            "p",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    U_
    (
        IOobject
        (
            "U",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    phi_
    (
        IOobject
        (
            "phi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::flux(U_)
    ),
    pB_
    (
        IOobject
        (
            "pB",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    B_
    (
        IOobject
        (
            "B",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    phiB_
    (
        IOobject
        (
            "phiB",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::flux(B_)
    ),
    pRefCell_(0),
    pRefValue_(0.0),
    cumulativeContErr_(0),
    sumLocalContErr_(0),
    globalContErr_(0)
{
    setRefCell(p_, piso_.dict(), pRefCell_, pRefValue_);

    mesh_.setFluxRequired(p_.name());
    mesh_.setFluxRequired(pB_.name());

    Info<< endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::mhdFoam::correctPhysics()
{
    correctFluidMechanics();
    correctMagneticField();
}


void Foam::solvers::mhdFoam::correctTightlyCoupledPhysics()
{}


void Foam::solvers::mhdFoam::correctFluidMechanics()
{
    #include "UEqn.H"

    // --- PISO loop
    while (piso_.correct())
    {
        #include "pEqn.H"
    }
}


void Foam::solvers::mhdFoam::correctMagneticField()
{
    // --- B-PISO loop
    while (bpiso_.correct())
    {
        #include "BEqn.H"
        #include "pBEqn.H"
    }

    printMagneticFieldErr();
}


void Foam::solvers::mhdFoam::correctCourant()
{
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    if (mesh_.nInternalFaces())
    {
        scalarField sumPhi
        (
            fvc::surfaceSum(mag(phi_))().primitiveField()
        );

        CoNum = 0.5*gMax(sumPhi/mesh_.V().field())*mesh_.time().deltaTValue();

        meanCoNum =
            0.5*(gSum(sumPhi)/gSum(mesh_.V().field()))*mesh_.time().deltaTValue();
    }

    Info<< "Courant Number mean: " << meanCoNum
        << " max: " << CoNum << endl;
}


void Foam::solvers::mhdFoam::correctContErr()
{
    volScalarField contErr(fvc::div(phi_));

    sumLocalContErr_ =
        mesh_.time().deltaTValue()
       *mag(contErr)().weightedAverage(mesh_.V()).value();

    globalContErr_ =
        mesh_.time().deltaTValue()
       *contErr.weightedAverage(mesh_.V()).value();

    calcCumulContErr();
}


void Foam::solvers::mhdFoam::calcCumulContErr()
{
    cumulativeContErr_ += globalContErr_;
}


void Foam::solvers::mhdFoam::printContErr()
{
    Info<< "time step continuity errors : sum local = " << sumLocalContErr_
        << ", global = " << globalContErr_
        << ", cumulative = " << cumulativeContErr_
        << endl;
}


void Foam::solvers::mhdFoam::printMagneticFieldErr()
{
    Info<< "magnetic flux divergence error = "
        << mesh_.time().deltaTValue()
          *mag(fvc::div(phiB_))().weightedAverage(mesh_.V()).value()
        << endl;
}


Foam::scalar Foam::solvers::mhdFoam::maxDeltaT()
{
    scalar newDeltaT =
        mesh_.time().controlDict().getOrDefault<scalar>("maxDeltaT", GREAT);

    scalar maxCo =
        mesh_.time().controlDict().getOrDefault<scalar>("maxCo", 1);

    scalar CoNum = 0.0;

    if (mesh_.nInternalFaces())
    {
        scalarField sumPhi
        (
            fvc::surfaceSum(mag(phi_))().primitiveField()
        );

        CoNum = 0.5*gMax(sumPhi/mesh_.V().field())*mesh_.time().deltaTValue();
    }

    if (CoNum > SMALL)
    {
        scalar maxDeltaTFact = maxCo/(CoNum + SMALL);
        scalar deltaTFact =
            min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

        newDeltaT =
            min(deltaTFact*mesh_.time().deltaTValue(), newDeltaT);
    }

    return newDeltaT;
}


// ************************************************************************* //
