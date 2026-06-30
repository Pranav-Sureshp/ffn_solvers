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

Class
    Foam::solvers::mhdFoam

Description
    foamForNuclear module wrapping the legacy OpenFOAM mhdFoam application:
    incompressible, laminar flow of a conducting fluid under the influence
    of an applied magnetic field B, with a fictitious magnetic flux
    pressure pB used to keep the magnetic face flux phiB divergence-free.

    Ported from applications/solvers/electromagnetics/mhdFoam (OpenFOAM)
    to the foamForNuclear solver/module run-time-selectable architecture,
    following the same pattern used for compressibleInterFoam.

SourceFiles
    mhdFoam.C
    UEqn.H
    pEqn.H
    BEqn.H
    pBEqn.H

\*---------------------------------------------------------------------------*/

#ifndef mhdFoam_H
#define mhdFoam_H

#include "fvCFD.H"
#include "pisoControl.H"
#include "solver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class mhdFoam Declaration
\*---------------------------------------------------------------------------*/

class mhdFoam
:
    public solver
{
private:

    //- Mesh reference
    fvMesh& mesh_;

    //- PISO control for the flow (U, p)
    pisoControl piso_;

    //- PISO control for the magnetic field (B, pB), subDict "BPISO"
    pisoControl bpiso_;

    //- Transport properties dictionary
    IOdictionary transportProperties_;

    //- Physical properties
    dimensionedScalar rho_;
    dimensionedScalar nu_;
    dimensionedScalar mu_;
    dimensionedScalar sigma_;

    //- Magnetic diffusivity, 1/(mu*sigma)
    dimensionedScalar DB_;

    //- Magnetic body-force coefficient, 1/(2*mu*rho)
    dimensionedScalar DBU_;

    //- Flow fields
    volScalarField p_;
    volVectorField U_;
    surfaceScalarField phi_;

    //- Magnetic fields
    volScalarField pB_;
    volVectorField B_;
    surfaceScalarField phiB_;

    //- Pressure reference cell/value for p_
    label pRefCell_;
    scalar pRefValue_;

    //- Continuity error bookkeeping
    scalar cumulativeContErr_;
    scalar sumLocalContErr_;
    scalar globalContErr_;


public:

    //- Runtime type information
    TypeName("mhdFoam");


    // Constructors

        //- Construct from mesh
        mhdFoam(dynamicFvMesh& mesh);

        //- Disallow default bitwise copy construction
        mhdFoam(const mhdFoam&) = delete;


    //- Destructor
    virtual ~mhdFoam()
    {}


    // Member Functions

        //- Solve the MHD physics: flow PISO loop followed by the
        //  magnetic-field B-PISO loop
        virtual void correctPhysics();

        //- No tightly-coupled sub-physics for this module
        virtual void correctTightlyCoupledPhysics();

        //- Maximum stable time step from the flow Courant number
        virtual scalar maxDeltaT();

        //- Solve the momentum predictor and the U-p PISO loop
        void correctFluidMechanics();

        //- Solve the induction equation and the B-pB PISO loop
        void correctMagneticField();

        //- Update and print the flow Courant number
        void correctCourant();

        //- Update the continuity error of the flow
        void correctContErr();

        //- Print continuity error information
        void printContErr();

        //- Accumulate the cumulative continuity error
        void calcCumulContErr();

        //- Print the magnetic flux divergence error
        void printMagneticFieldErr();

        virtual scalar getResidual()
        {
            return scalar(0);
        }

        virtual void correctBaffleLessFields()
        {}

        virtual void deformMesh()
        {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const mhdFoam&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
