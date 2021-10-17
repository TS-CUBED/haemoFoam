/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1812                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*\
Application
    haemoSimpleFoam

Description
    Steady state solver for incompressible, turbulent flow
    of non-Newtonian high particle load suspension fluids, using
    the SIMPLE algorithm.

    Developed for modelling the shear thinning and hematocrit transport
    in blood flow.

    Uses particle transport
    equation based on (Phillips, Armstrong et al. 1992):

    PHILLIPS, R.J., ARMSTRONG, R.C., BROWN, R.A., GRAHAM, A.L. and
    * ABBOTT, J.R., 1992.
       A constitutive equation for concentrated suspensions that accounts for
       shear‐induced particle migration.
       Physics of Fluids A: Fluid Dynamics, 4(1), pp. 30-40.

    QUEMADA, D., 1978. Rheology of concentrated disperse systems III.
        General features of the proposed non-newtonian model.
        Comparison with experimental data.
        Rheologica Acta, 17(6), pp. 643-653.

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
#ifdef OPENFOAMESIORFOUNDATION
#include "turbulentTransportModel.H"
#else
#include "turbulenceModel.H"
#endif
#include "simpleControl.H"
#ifdef OPENFOAMESIORFOUNDATION
#include "fvOptions.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Steady-state solver for incompressible, turbulent flows."
    );

#ifdef OPENFOAMESIORFOUNDATION
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    turbulence->validate();
#else
   include "setRootCase.H"
   include "createTime.H"
   include "createMesh.H"

   simpleControl simple(mesh);

   include "createFields.H"
   include "initContinuityErrs.H"
#endif

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Read properties for haematocrit transp

#   include "readHaemoProperties.H"

    turbulence->correct();

    // Back to original code

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Pressure-velocity SIMPLE corrector
        {
#           include "UEqn.H"
#           include "pEqn.H"


        }

        // Calculate H equation

#	    include "HEqn.H"

        if (haemoSwitch.value() == 0 || runTime < haemoSwitchTime)
        {
            Info<< "Not Solving for H, Migration Model is inactive 1" << nl << endl;
        } else
        {
            HEqn.relax();                   // use these two lines
            HEqn.solve().initialResidual(); // for underrelaxed solver
            // HEqn.solve(); // or this one for no underrelaxation (unstable)
        }

        laminarTransport.correct();
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
