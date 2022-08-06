/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1812                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*\
Application
    Any custom solver with Windkessel boundary conditions

Description
    Header file for the initialisation of Windkessel properties.
    Loads the WindkesselProperties dictionary
    Extracts the number of outlets in the system
    Updates the global variable N_OUTLETS
    Assigns a scalarList for storage of pressure values
    Finds the timestep and updates the global variable dt
    Runs the function intiialiseWindkessel.C; creates the WK struct

Dependencies
    WKFunctions.C
    windkesselProperties dictionary (when running)
\*---------------------------------------------------------------------------*/


/* Load dictionary of Windkessel properties */

IOdictionary windkesselProperties
(
    IOobject
    (
        "windkesselProperties",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        // IOobject::MUST_READ_IF_MODIFIED,
        IOobject::AUTO_WRITE
    )
);


/* Finding the number of outlets */

const wordList outletNames(windkesselProperties.toc());
int sizeOf = outletNames.size();
N_OUTLETS = sizeOf;
if (N_OUTLETS > 0)
{
    Info << "Number of windkessel outlets are :" << N_OUTLETS << endl;
} else {
    Info << "File constant/windkesselProperties not found. If you use WKBC, please define the properties in that file." << endl;
}

/* Read flag */
int WK_FLAG = N_OUTLETS;



/* Define timestep and density for use with Windkessel */

dt = runTime.deltaTValue();

IOdictionary transportProperties_extra
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar rho
(
    "rho",
    transportProperties_extra.lookup("rho")
);


RHO_0 = rho.value();


/* Updating global array patch_names */

// if (WK_FLAG >0)
// {

 forAll(outletNames, item)
   {

     const word& outletName = outletNames[item];
     const dictionary& subDict = windkesselProperties.subDict(outletName);
     int realIndex = readScalar(subDict.lookup("outIndex"));

     patch_names[realIndex] = outletName;
   }

for (int it_check = 0; it_check < N_OUTLETS; it_check++)
  {

    Info << "Patch array " << it_check << " is " << patch_names[it_check] << endl;

  }

/* Creating scalar list for storage of pressure values*/
if (N_OUTLETS > 0)
{
    Info<< "Reading scalar list WKpressures \n" << endl;
}

scalarIOList WKpressures
(
    IOobject
    (
        "WKpressures",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    N_OUTLETS
);

/* Initialise the WK properties*/
initialise(windkesselProperties);


/*Initialising the storage array*/

for (int it = 0; it <N_OUTLETS; it++)
{
	WKpressures[it] = wk[it].P_0;
	//Info << "First Stored " << store[it] << endl;
}
// }
