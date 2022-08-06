`haeomFoam` is a Finite Volume code for blood flow simulations. It is
written in \`C++\` using the OpenFOAM framework and includes some
additional features:

  - non-Newtonian rheology models, in particular the Quemada model
    ([4](#citeproc_bib_item_7)) and a modified Krieger-Dougherty model
    ([4](#citeproc_bib_item_1)), as well as the original Krieger model
    ([4](#citeproc_bib_item_2)).
  - haematocrit transport model based on ([4](#citeproc_bib_item_3))
  - Haematocrit Transport Model ([4](#citeproc_bib_item_9))
  - Windkessel BCs ([4](#citeproc_bib_item_10)), adapted from code
    written by colleagues at Imperial College, London
    ([4](#citeproc_bib_item_5))

# Installation

`haemoFoam` is developed on [ESI OpenFOAM](https://www.openfoam.com).
The current development version is v2106. It will not work on Foundation
OpenFOAM, or of FOAM Extend. There are some conditional compiler
statements in the code, which suggest that it does, but these versions
are no longer maintained at the moment.

Installation of `haemoFoam` requires a full installation of OpenFOAM,
including the development libraries, and C++ compiler. It has only been
tested on [Ubuntu Linux LTS](https://www.ubuntu.com) (or derivatives
like [LinuxMint](https://www.linuxmint.com), or [Pop\!-OS
LTS](https://pop.system76.com/) with OpenFOAM installed from the
[PPA](https://develop.openfoam.com/Development/openfoam/-/wikis/precompiled/debian).
It should work on any other Linux based OpenFOAM installation, though
(self-compiled on CentOS7 HPC cluster is one of the test-systems, as
well).

  - Download the latest release of `haemoFoam` (if you want to
    contribute to development, you can also clone the github repository,
    but keep in mind that I will update the code regularly, so things
    may break\!).

  - Unpack the archive into
    `~/OpenFOAM/yourUserName-vYYMM/applications/`

  - Compile `haemoFoam` after initialising OpenFOAM:
    
    ``` bash
    openfoam2106
    
    cd ~/OpenFOAM/$USER-v2106/applications/haemoFoam-0.2.8
    
    ./Allwmake
    ```
    
    This will install `haemoFoam` for the current user in
    `~/OpenFOAM/$USER-v2106/platforms/linux.../bin`

## Binaries

`haemoFoam` comes with three executables:

1.  `haemoSimpleFoam` - steady state solver
2.  `haemoPimpleFoam` - transient solver, can be run in `PIMPLE` or
    `PISO` mode
3.  `haemoPostProcess` - post-processing utility, calculates
    wall-shear-stress derived metrics. Needs `wallShearStress` fields to
    be calculated either at run-time (see Function Objects) or using the
    `-postProces` option\[fn::Note that the postProcessing routine does
    not work with the windkessel boundary condition at the moment, so
    using the function object to write the `wallShearStress` field at
    run-time is the recommended way.)

## Models

`haemoFoam` includes the following models:

### Viscosity models

  - **Carreau** (`libCarreau.so`) - implementation of the Carreau model
    as implemented in Fluent (to allow for one-to-one comparisons with
    Fluent). This model does NOT use the haematocrit transport
    model\!\[1\]
  - **KriegerDougherty** (`libKriegerDougherty.so`) - traditional
    Krieger-Dougherty model ([4](#citeproc_bib_item_2)).
  - **Krieger5** (`libKrieger5.so`) - a modified 5-parameter Krieger
    model after ([4](#citeproc_bib_item_1)).
  - **Yeleswarapu** (`libYeleswarapu.so`) - Carreau-type model after
    ([4](#citeproc_bib_item_12)),([4](#citeproc_bib_item_11)).
  - **Quemada** (`libQuemada.so`) - the main development model. Based on
    ([4](#citeproc_bib_item_6))([4](#citeproc_bib_item_8))([4](#citeproc_bib_item_7)).

All but the Carreau model use the local haemotocrit and shear rate to
calculate the viscosity.

Refer to ([4](#citeproc_bib_item_9)) for details on the implementation.

### Boundary conditions

  - **WKBC** - 3-element windkessel boundary condition based on code
    written by Andris Piebalgs, Boram Gu, Emily Manchester, Imperial
    College London, who kindly let me have their [original
    code](https://github.com/KeepFloyding/OpenFOAM-phys-flow).
    Implemented as part of `haemoPimpleFoam`\[2\].
  - **splitFlowRateOutletVelocity**
    (`libSplitFlowRateOutletVelocity.so`) - A split flow rate outlet
    condition. Not recommended\[3\], but included since there are many
    models out there that use this kind of BC. Forces outlet flow rate
    at the outlet to match a chosen percentage of the inlet flow rate.

# Case setup

`haemoFoam` cases are set up in a similar fashion as a regular
`pimpleFoam` case. I will assume here, that you know how to do that and
will focus on the differences:

## Loading the model libraries

The new models need to be loaded in `system/controlDict` in the `libs`
section. E.g., to load the Quemada model, the flow split, and the
Groovy\[4\] BC, the `controlDict` file will start with:

``` example
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1812                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/


FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

libs
(
    "libQuemada.so"
    "libgroovyBC.so"
    "libsplitFlowRateOutletVelocity.so"
);


application     haemoPimpleFoam;
```

## Haematocrit field

### BCs

The H (haematocrit) field needs an initial and boundary condition file
`0/H`. This contains the volume fraction of haematocrit (dimensionless).

  - Inlet: `fixedValue` - give average H value. Alternatively you can
    give an H profile using something like `groovyBC`.

  - Outlets: `inletOutlet` - should there be backflow on the outlet, we
    need to provide the H value for the inflow (otherwise this will be
    zero, which will cause instability).

  - Walls: I use `slip` to allow arbitrary gradients, `fixedGradients`
    should work as well.

<!-- end list -->

``` example
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2106                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/


FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      H;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 0.45;

boundaryField
{
    ICA
    {
        type            inletOutlet;
        inletValue      uniform 0.45;
        value           uniform 0.45;
    }
    ECA
    {
        type            inletOutlet;
        inletValue      uniform 0.45;
        value           uniform 0.45;
    }
    WALL
    {
        type            slip;
    }
    APEX
    {
        type            slip;
    }
    SINUS
    {
        type            slip;
    }
    CCA
    {
        type            fixedValue;
        value           uniform 0.45;
    }
}


// ************************************************************************* //
```

### Solver settings

The haemotocrit field is simulated using a transport equation. This
means it needs discretisation schemes for all terms set in `fvSchemes`
and convergence criteria in `fvSolution`. Refer to the tutorial case for
examples.

## Boundary conditions

Most OpenFOAM boundary conditions are available in `haemoFoam` as well.

However, there is a limitation in the implementation of the windkessel
boundary condition in `haemoFoam`. Even if this BC is not needed, the
configuration file (`windkesselProperties`) needs to be present in the
`constant` directory\[5\].

Here I only describe the special boundary conditions in `haemoFoam`.

### Windkessel Boundary Conditions

Windkessel BCs were implemented using code developed by Andris Piebalgs
at Imperial College (,Piebalgs [4](#citeproc_bib_item_4)).

The BCs are set in the boundary/initial condition directory `0` and in
the `constant` directory in two files:

  - the `p` boundary condition file (in `0`), where each WK outlet gets
    an entry of the form
    
        outlet_name
          {
              type            WKBC;
              index           0;
              value           uniform 0;
          }

    where the `index` is a number that has to match the entry in
    
      - the `windkesselProperties` file (in `constant`), where the
        parameters are set:

<!-- end list -->

```example
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1812                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/


FoamFile {
  version 4.0;
  format ascii;
  class dictionary;
  object windkesselProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

ICA {
  Z 0.1;    // commonly R_1 or R_proximal
  C 0.5;
  R 2.0;   // commonly R_2 or R_distal

  // NOTE: Set this to true to use physiological units [mmHg.ml^-1.s] and [ml.mmHg^-1]
  //       if this is false, then SI units are used. Physiological units are converted
  //       at runtime.
  physiologicalUnits true;

  outIndex 0; // must equal 'index' value in 0/p
  FDM_order 1;
  // backward difference order: up to 3rd order

  // Initialise WK parameters
  // also useful to set initial condidions to speed up pressure development
  Flowrate_threeStepBefore 0;
  Flowrate_twoStepBefore 0;
  Flowrate_oneStepBefore 0;
  // NOTE: these pressures are in [Pa]!
  Pressure_twoStepBefore 9000;
  Pressure_oneStepBefore 9000;
  Pressure_start 9000;
}

ECA {
  Z 0.2;
  C 1.5;
  R 0.6;
  
  physiologicalUnits true;

  outIndex 1;
  FDM_order 1;
  Flowrate_threeStepBefore 0;
  Flowrate_twoStepBefore 0;
  Flowrate_oneStepBefore 0;
  Pressure_twoStepBefore 9000;
  Pressure_oneStepBefore 9000;
  Pressure_start 9000;
}

// ************************************************************************* //
```

**Important Notes:**

  - The order of the finite difference method to solve the windkessel
    equation can be set to up to third order. However for realistic
    cases, first order is sufficient, and higher orders can become
    unstable.
  - The parameters for the windkessel properties can be given in 
    physiological or in SI units.
    In most publications these will be given in physiological units,
    e.g., [ mm_{Hg}.ml^{-1}.s ]. These need to be
    converted to SI-units\!
  - If the windkessel properties are not available, fixed relative
    pressures at the outlets are a reasonable alternative. A simple
    `fixedValue` BC with `value uniform 0` will give realistic flow
    distributions in most cases (based on diameter scaling).\]

### Flow split boundary conditions

These BCs are implemented to allow comparisons to some published cases
that use Fluent split flow outlet condition. These BCs are typically
unphysiological and can deliver unreliable results. They are typically
employed where there are data available at the outlets, e.g., from echo
or MRI flux measurements.

The flow split boundary condition works in cases where there is a single
inlet\[6\] and multiple outlets. One outlet needs to be set as a free
outlet (`zeroGradient`) and the others are set as
`splitFlowRateOutletVelocity` in `0/U`:

``` example
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1812                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/


FoamFile
{
    version     2.0;
    format      ascii;
    class       volVectorField;
    location    "0";
    object      U;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 1 -1 0 0 0 0];

internalField   uniform (0 0 0.1);

boundaryField
{
    ICA
    {
        type            splitFlowRateOutletVelocity;
        inletPatch      "CCA";
        flowSplit       0.66;
        value           uniform ( 0 0 0.1 );
    }
    ECA
    {
        type            zeroGradient;
    }


....
```

The pressure BCs for this case are `zeroGradient` for the split flow
outlets, and `fixedValue` for the free outlet (in `0/p`):

``` example
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1812                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/


FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      p;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -2 0 0 0 0];

internalField   uniform 8.5;

boundaryField
{
    ICA
    {
        type            zeroGradient;
    }
    ECA
    {
        type            fixedValue;
        value           uniform 0;
    }

....
```

## Function Objects

There are a few things that we routinely look at when simulating blood
flows. To save postprocessing time and file space, many of these can be
written to text files at runtime. These can then be analysed and plotted
without the need to load the field data.

The Wall Shear stress is a required one, the others are optional and are
just the ones that I routinely use:

### Wall Shear Stress (required for `haemoPostProcess`):

Wall shear stress is an important metric for haemodynamic simulation.
OpenFOAM does not calculate WSS, but has a function object that will do
that. This will take turbulence and non-Newtonian models into account
and calculate WSS as
\(\nu_{eff} \frac{\partial\vec{u}}{\partial\vec{n}}\).

``` example
wallShearStress
{
    // Mandatory entries (unmodifiable)
    type            wallShearStress;
    libs            (fieldFunctionObjects);

    // Optional entries (runtime modifiable)
    // patches         (WALL APEX SINUS); // (wall1 "(wall2|wall3)");

    // Optional (inherited) entries
    // writePrecision  8;
    // writeToFile     true;
    // useUserTime     true;
    // region          region0;
    // enabled         true;
    // log             true;
    // timeStart       0;
    // timeEnd         1000;
    // executeControl  timeStep;
    // executeInterval 1;
    writeControl    outputTime;

}
```

*Note: `writeControl outputTime;` synchronises the output of the field
data with the output times of the rest of the data. Without this the
wall shear stress data will be written every time step\!*

`./postProcessing/wallShearStress/0/wallShearStress.dat` will contain
min and max values for all patches - in this case also for the
inlet/outlet patches, which can be ignored when plotting WSS in
Paraview. In some cases these help avoid division by zero problems at
the edges of the walls when calculating the WSS metrics.

``` example
# Wall shear stress
# Time          patch           min             max
0.01            CCA     (-2.285838e-03 -2.083722e-03 -1.215122e-03)     (1.951771e-03 1.997152e-03 1.098803e-03)
0.01            ICA     (-4.724931e-03 -5.578844e-03 -2.057902e-03)     (5.158803e-03 5.917156e-03 1.551159e-03)
0.01            ECA     (-7.802540e-03 -5.654608e-03 -1.776018e-03)     (5.327658e-03 8.974748e-03 3.580220e-04)
0.01            APEX    (-7.673434e-03 -2.424255e-03 -5.991218e-03)     (2.045973e-03 1.787954e-03 -4.718725e-04)
0.01            SINUS   (-1.249938e-03 -1.057508e-03 -7.173346e-03)     (2.467540e-03 1.757495e-03 1.624837e-03)
0.01            WALL    (-4.859152e-03 -2.770018e-03 -1.472179e-02)     (3.020219e-03 1.545451e-03 1.980882e-03)
```

### Average values on patches (e.g. inlet and outlet pressures)

Create a section in the `functions` block of `controlDict` for all
patches of interest; all fields can be averaged.

E.g. averaging over the common carotid arterie outflow (patch name
`CCA`):

``` example
averagesCCA
{
    type            surfaceFieldValue;
    libs            ("libfieldFunctionObjects.so");
    enabled         yes;
    writeControl    timeStep;
    log             yes;
    writeFields     no;
    regionType      patch;
    name            CCA;
    operation       areaAverage;

    fields
    (
        p
        U
        H
    );
}
```

This file will be written to
`./postProcessing/averagesCCA/<startTime>/0/surfaceFieldValue.dat` and
look like this:

``` example
# Region type : patch CCA
# Faces       : 731
# Area        : 3.089305e-05
# Scale factor: 1.000000e+00
# Time          areaAverage(p)  areaAverage(U)  areaAverage(H)
0.001           2.747977e+01    (-2.650558e-09 -3.238202e-09 2.523065e-01)      4.500000e-01
0.002           2.370257e+00    (-2.673728e-09 -3.266508e-09 2.545120e-01)      4.500000e-01
0.003           9.029901e+00    (-2.697549e-09 -3.295612e-09 2.567795e-01)      4.500000e-01
0.004           9.070692e+00    (-2.722026e-09 -3.325515e-09 2.591095e-01)      4.500000e-01
0.005           9.100255e+00    (-2.747160e-09 -3.356221e-09 2.615020e-01)      4.500000e-01
```

Note that the velocity is a vector in parenthesis. These can cause
problems with some import filters. A helpful command for removing these
parenthesis is `sed -e "s/[()]//g" surfaceFieldValue.dat >
surfaceFieldValue_cleaned.dat`, or similar.

### Flow Rates (volume flow)

These are similar to the averages. The flow rate can be calculated as
the `sum` over the scalar flow `phi`. Calculating the outflow from the
external carotid artery (ECA):

``` example
flowRateECA
   {
        type            surfaceFieldValue;
        libs ("libfieldFunctionObjects.so");
        enabled         yes;
        writeControl    timeStep;
        writeInterval   1;
        log             yes;
        writeFields     no;
        regionType      patch;
        name            ECA;
        operation       sum;

        fields
        (
            phi
        );
   }
```

*Note: due to the orientation of the normal vector and since the volume
flow is \((\vec{u} \cdot \vec{n})A\) flow into the volume is negative
and flow out of the volume is positive\!*

This file will be written to
`./postProcessing/flowRateCCA/<startTime>/0/surfaceFieldValue.dat` and
look like this:

``` example
# Region type : patch CCA
# Faces       : 731
# Area        : 3.089305e-05
# Scale factor: 1.000000e+00
# Time          sumphi
0.001           -7.794516e-06
0.002           -7.862651e-06
0.003           -7.932703e-06
0.004           -8.004682e-06
0.005           -8.078593e-06
```

### Solver info

Some additional information, like residuals, can be written to a text
file for each time step. Note that this is only of limited use for the
PIMPLE solver, since it will use the residual information in the first
PIMPLE loop, not the last.

``` example
solverInfo
{
    type            solverInfo;
    libs            ("libutilityFunctionObjects.so");
    fields          (U H p);
    timeStart       3.01;
    writeResidualFields no;
}
```

*Note: Adjust the `timeStart` to be greater than the `haemoSwitch` value
in `haemoTransportProperties`, so the `H` solver will run. Otherwise the
`H` residual will not be recorded\!*

# Post-Processing

## Calculation of WSS metrics

`haemoFOAM` has postprocessing routines for the most important
wall-shear-stress (WSS) derived metrics that are used in atherosclerosis
risk prediction.

These require data over the last cycle of the simulation. In the
examples the cycle length is 1s and data was written every 0.01s, data
for the 6th cycle (5.01s to 6s) is evaluated.

The postprocessing is done in two steps (only one for serial runs):

1.  `reconstructPar` to reconstruct the data from the distributed run
    (only reconstructing last cycle, `processor` directories can be
    removed after this step to save hard drive space):
    
    \#+begin<sub>example</sub> sh reconstructPar -time 5.01:6
    
    \#+end<sub>example</sub>

2.  `haemoPostProcess` calculates the other shear stress metrics:
    
    ``` example
    haemoPostProcess -time 5.01:6
    ```

# Bibliography

<span id="citeproc_bib_item_1"></span>Hund, Samuel, Marina Kameneva, and
James Antaki. 2017. “A Quasi-Mechanistic Mathematical Representation for
Blood Viscosity.” *Fluids* 2 (1): 10.
<https://doi.org/10.3390/fluids2010010>.

<span id="citeproc_bib_item_2"></span>Krieger, Irvin M., and Thomas J.
Dougherty. 1959. “A Mechanism for Non-Newtonian Flow in Suspensions of
Rigid Spheres.” *Transactions of the Society of Rheology* 3 (1): 137–52.
<https://doi.org/10.1122/1.548848>.

<span id="citeproc_bib_item_3"></span>Phillips, Ronald J., Robert C.
Armstrong, Robert A. Brown, Alan L. Graham, and James R. Abbott. 1992.
“A Constitutive Equation for Concentrated Suspensions That Accounts
for Shear-Induced Particle Migration.” *Physics of Fluids a: Fluid
Dynamics* 4 (1): 30–40. <https://doi.org/10.1063/1.858498>.

<span id="citeproc_bib_item_4"></span>Piebalgs, Andris. 2017.
“Development of a Multi-Physics and Multi-Scale Model of Thrombolysis
for Patient-Specific Appplications,” May.
<https://doi.org/10.25560/67674>.

<span id="citeproc_bib_item_5"></span>Pirola, S., Z. Cheng, O. A.
Jarral, D. P. O’Regan, J. R. Pepper, T. Athanasiou, and X. Y. Xu. 2017.
“On the Choice of Outlet Boundary Conditions for Patient-Specific
Analysis of Aortic Flow Using Computational Fluid Dynamics.” *Journal of
Biomechanics* 60 (July): 15–21.
<https://doi.org/10.1016/j.jbiomech.2017.06.005>.

<span id="citeproc_bib_item_6"></span>Quemada, D. 1977. “Rheology of
Concentrated Disperse Systems and Minimum Energy Dissipation Principle -
I. Viscosity-Concentration Relationship.” *Rheologica Acta* 16: 82–94.
<https://doi.org/10.1007/BF01516932>.

<span id="citeproc_bib_item_7"></span>———. 1978a. “Rheology of
Concentrated Disperse Systems III. General Features of the Proposed
Non-Newtonian Model. Comparison with Experimental Data.” *Rheologica
Acta* 17 (6): 643–53. <https://doi.org/10.1007/BF01522037>.

<span id="citeproc_bib_item_8"></span>———. 1978b. “Rheology of
Concentrated Disperse Systems II. A Model for Non-Newtonian Shear
Viscosity in Steady Flows.” *Rheologica Acta* 17 (6): 632–42.
<https://doi.org/10.1007/BF01522036>.

<span id="citeproc_bib_item_9"></span>Schenkel, Torsten, and Ian
Halliday. 2021. “Continuum Scale Non Newtonian Particle Transport Model
for Hæ morheology.” *Mathematics* 9 (17): 2100.
<https://doi.org/10.3390/math9172100>.

<span id="citeproc_bib_item_10"></span>Westerhof, Nico, Jan-Willem
Lankhaar, and Berend E. Westerhof. 2009. “The Arterial Windkessel.”
*Medical & Biological Engineering & Computing* 47 (2): 131–41.
<https://doi.org/10.1007/s11517-008-0359-2>.

<span id="citeproc_bib_item_11"></span>Wu, Wei-Tao, Fang Yang, James F.
Antaki, Nadine Aubry, and Mehrdad Massoudi. 2015. “Study of Blood Flow
in Several Benchmark Micro-Channels Using a Two-Fluid Approach.”
*International Journal of Engineering Science* 95 (Journal Article):
49–59. <https://doi.org/10.1016/j.ijengsci.2015.06.004>.

<span id="citeproc_bib_item_12"></span>Yeleswarapu, K. K., M. V.
Kameneva, K. R. Rajagopal, and J. F. Antaki. 1998. “The Flow of Blood in
Tubes: Theory and Experiment.” *Mechanics Research Communications* 25
(3): 257–62. <https://doi.org/10.1016/S0093-6413(98)00036-6>.

1.  \*Bird-Carreau\* is available as an OpenFOAM model. All other OF
    viscosity models can be used as well. They will, however not link to
    the haematocrit model.

2.  This means that the `windkesselProperties` file must exist, even if
    these BC is not used. Will need to reimplement this in a cleaner
    way.

3.  This BC is unphysiological. I would discourage it's use. It seems to
    have fallen out of favour recently, but there are many simulations -
    including some of mine - that use this BC if there is data from,
    e.g., echo doppler or similar. It can work in these cases, but I
    would recommend fixed pressure BCs if no windkessel data is
    available.

4.  groovyBC is part of the "Swiss Army Knife for Foam" (swak4Foam)
    toolkit, which I recommend.

5.  I will need to address this in an upcoming refactoring, but at the
    moment, I have decided to just leave it as is.

6.  could be extended to multiple inlets
