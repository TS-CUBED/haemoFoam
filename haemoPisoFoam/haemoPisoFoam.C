/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    haemoPisoFoam

Description
    Transient solver for incompressible, turbulent flow 
    of non-Newtonian high particle load suspension fluids, using
    the PISO algorithm.
    
    Developed for modelling the shear thinning and hematocrit transport
    in blood flow.
    
    Uses particle transport
    equation based on (Phillips, Armstrong et al. 1992):
    
    PHILLIPS, R.J., ARMSTRONG, R.C., BROWN, R.A., GRAHAM, A.L. and 
    * ABBOTT, J.R., 1992.
       A constitutive equation for concentrated suspensions that accounts for 
       shear‐induced particle migration.
       Physics of Fluids A: Fluid Dynamics, 4(1), pp. 30-40.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.
    * although particle transport model is not developed for turbulent
    * flow models.

    Consistent formulation without time-step and relaxation dependence by Jasak

Author and Copyright
    Dr Torsten Schenkel
    Department Engineering and Mathematics
    Material and Engineering Research Institute MERI
    Sheffield Hallam University
    February 2019
    All Rights Reserved
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    pisoControl piso(mesh);

#   include "createFields.H"
#   include "initContinuityErrs.H"

    // Read properties for haematocrit transp

#   include "readHaemoProperties.H"

    turbulence->correct();


    // Back to original code


    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "CourantNo.H"

        // Pressure-velocity PISO corrector
        {
            // Momentum predictor

            fvVectorMatrix UEqn
                (
                 fvm::ddt(U)
                 + fvm::div(phi, U)
                 + turbulence->divDevReff()
                );

            UEqn.relax();

            if (piso.momentumPredictor())
            {
                solve(UEqn == -fvc::grad(p));
            }

            // --- PISO loop

            while (piso.correct())
            {
                volScalarField rUA = 1.0/UEqn.A();

                U = rUA*UEqn.H();
                phi = (fvc::interpolate(U) & mesh.Sf())
                    + fvc::ddtPhiCorr(rUA, U, phi);

                adjustPhi(phi, U, p);

                // Non-orthogonal pressure corrector loop
                while (piso.correctNonOrthogonal())
                {
                    // Pressure corrector

                    fvScalarMatrix pEqn
                        (
                         fvm::laplacian(rUA, p) == fvc::div(phi)
                        );

                    pEqn.setReference(pRefCell, pRefValue);
                    pEqn.solve
                        (
                         mesh.solutionDict().solver
                         (
                          p.select(piso.finalInnerIter())
                         )
                        );

                    if (piso.finalNonOrthogonalIter())
                    {
                        phi -= pEqn.flux();
                    }

                    //                  Calculate H equation                    

#                   include "HEqn.H"

                    //                  Back to original code

                }

#               include "continuityErrs.H"

                U -= rUA*fvc::grad(p);
                U.correctBoundaryConditions();
            }
        }

        turbulence->correct();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
