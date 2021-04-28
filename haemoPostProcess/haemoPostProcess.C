/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License

    This is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    This code is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    haemoPostProcess

Description
    Calculates and writes WSS derived metrics for haemodynamic cases


    WSS: WSS, gradient of U at the wall, multiplied with local
    * kinematic viscosity * density

    avgWSS: average WSS for the time range

    OSI: Oscillatory Shear Stress Index

    transWSS: Transverse WSS

    All WSS calculations are using kinematic values in
    * accordance to FOAMs stress and pressure definitions and are
    * multiplied with the density rho as defined in transportProperties
    * to get the commonly used dynamic values.

    Be aware that pressures in OpenFoam are kinematic pressures!

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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createMesh.H"

#   include "readHaemoProperties.H"

    IOdictionary transportProperties
        (
         IOobject
         (
          "transportProperties", // name of the dictionary
          runTime.constant(), // location in the case - this one is in constant
          mesh, // needs the mesh object reference to do some voodoo - unimportant now
          IOobject::MUST_READ,
          IOobject::NO_WRITE // read-only
         )
        );

    // const dimensionedScalar rho(transportProperties.lookup("rho"));
    dimensionedScalar rho
    (
        // "rho",
        // viscosityProperties.lookup("rho")
        "rho",
        transportProperties.lookup("rho")
    );

    // Info << rho << endl;


    volVectorField normalVector
        (
         IOobject
         (
          "normalVector",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
         ),
         mesh,
         dimensionedVector
         (
          "normalVector",
          dimLength,
          vector::zero
         )
        );

    volVectorField TAWSS
        (
         IOobject
         (
          "TAWSS",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
         ),
         mesh,
         dimensionedVector
         (
          "TAWSS",
          dimMass/(dimLength*sqr(dimTime)),
          //sqr(dimLength)/sqr(dimTime),
          vector::zero
         )
        );

    volScalarField transWSS
        (
         IOobject
         (
          "transWSS",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
         ),
         mesh,
         dimensionedScalar
         (
          "transWSS",
          dimMass/(dimLength*sqr(dimTime)),
          //sqr(dimLength)/sqr(dimTime),
          0
         )
        );

    volScalarField TAWSSMag
        (
         IOobject
         (
          "TAWSSMag",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
         ),
         mesh,
         dimensionedScalar
         (
          "TAWSSMag",
          dimMass/(dimLength*sqr(dimTime)),
          //sqr(dimLength)/sqr(dimTime),
          0
         )
        );

    volScalarField OSI
        (
         IOobject
         (
          "OSI",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
         ),
         mesh,
         dimensionedScalar
         (
          "OSI",
          dimless,
          0
         )
        );

    volScalarField RRT
        (
         IOobject
         (
          "RRT",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
         ),
         mesh,
         dimensionedScalar
         (
          "RRT",
          (dimLength*sqr(dimTime))/dimMass,
          //sqr(dimTime)/sqr(dimLength),
          0
         )
        );

    volScalarField TAHct
        (
         IOobject
         (
          "TAHct",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
         ),
         mesh,
         dimensionedScalar
         (
          "TAHct",
          dimless,
          0
         )
        );



    int nfield = 0;

    // First run, calculate WSS, avgWSS and OSI

    Info<< "First Run - calculating WSS, average WSS, OSI, RRT" << endl;

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << "    " << endl;

        IOobject Uheader
            (
             "U",
             runTime.timeName(),
             mesh,
             IOobject::MUST_READ
            );
        
        IOobject U_0header
            (
             "U_0",
             runTime.timeName(),
             mesh,
             IOobject::MUST_READ
            );
        
        IOobject nuheader
            (
             "nu",
             runTime.timeName(),
             mesh,
             IOobject::MUST_READ
            );
            
        IOobject Hheader
           (
            "H",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
           );


        // // Check U exists
        // if (Uheader.headerOk())
        // {
        //     if (nuheader.headerOk())
        //     {
                mesh.readUpdate();

                Info<< "    Reading U" << endl;
                volVectorField U(Uheader, mesh);

                Info<< "    Reading nu" << endl;
                volScalarField nu(nuheader, mesh);
                
                Info<< "    Reading H" << endl;
                volScalarField H(Hheader, mesh);

                Info<< "    Calculating WSS" << endl;

                volVectorField WSS
                    (
                     IOobject
                     (
                      "WSS",
                      runTime.timeName(),
                      mesh,
                      IOobject::NO_READ,
                      IOobject::AUTO_WRITE
                     ),
                     mesh,
                     dimensionedVector
                     (
                      "WSS",
                      dimMass/(dimLength*sqr(dimTime)),
                      vector::zero
                     )
                    );

                volScalarField WSSMag
                    (
                     IOobject
                     (
                      "WSSMag",
                      runTime.timeName(),
                      mesh,
                      IOobject::NO_READ,
                      IOobject::AUTO_WRITE
                     ),
                     mesh,
                     dimensionedScalar
                     (
                      "WSSMag",
                      dimMass/(dimLength*sqr(dimTime)),
                      0
                     )
                    );

                forAll(WSS.boundaryField(), patchi)
                {

                    WSS.boundaryFieldRef()[patchi] =
                         -U.boundaryField()[patchi].snGrad()
                         * nu.boundaryField()[patchi]
                         * rho.value();

                    WSSMag.boundaryFieldRef()[patchi] =
                        mag(-U.boundaryField()[patchi].snGrad()
                        * nu.boundaryField()[patchi])
                        * rho.value() ;
                        
                    TAWSS.boundaryFieldRef()[patchi] +=
                         -U.boundaryField()[patchi].snGrad()
                        * nu.boundaryField()[patchi]
                        * rho.value() ;

                    TAWSSMag.boundaryFieldRef()[patchi] +=
                        mag(-U.boundaryField()[patchi].snGrad()
                        * nu.boundaryField()[patchi])
                        * rho.value() ;
                        
                    TAHct.boundaryFieldRef()[patchi] +=
                        H.boundaryField()[patchi] ;
                }
                
                WSS.write();
                WSSMag.write();

                nfield++;
        //     }
        //     else
        //     {
        //         Info<< "    No nu" << endl;
        //     }
        // }
        // else
        // {
        //     Info<< "    No U" << endl;
        // }

        // Check U_0 exists for TWSSG
        // if (U_0header.headerOk())
        // {
        //     if (nuheader.headerOk())
        //     {
                mesh.readUpdate();
    }
    
    Info<< "Writing TAWSS, TAWSSMag, OSI, RRT, normalVectors" << endl;
    
    // devide by the number of added fields
    if(nfield>0){
        Info<< "number of fields added: "<< nfield << endl;
    
        forAll(TAWSS.boundaryField(), patchi)
        {
            TAWSS.boundaryFieldRef()[patchi] /= nfield;
    
            TAWSSMag.boundaryFieldRef()[patchi] /= nfield;
            
            TAHct /= nfield ;
    
            OSI.boundaryFieldRef()[patchi] =
                0.5
                * ( 1 - mag(TAWSS.boundaryField()[patchi])
                        /(TAWSSMag.boundaryField()[patchi]+1e-64)
                  );
    
            RRT.boundaryFieldRef()[patchi] =
                1
                /((
                            (1 - 2.0*OSI.boundaryField()[patchi])
                            *
                            TAWSSMag.boundaryField()[patchi]
                  )+1e-64);
    
            normalVector.boundaryFieldRef()[patchi] =
                mesh.Sf().boundaryField()[patchi]
                /mesh.magSf().boundaryField()[patchi]
                ;
        }
    }
    
    
    
    // Second run, calculate transWSS
    
    Info<< "Second Run - calculating transWSS" << endl;
    
    nfield = 0;
    
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << "    " << endl;
    
        IOobject WSSheader
            (
             "WSS",
             runTime.timeName(),
             mesh,
             IOobject::MUST_READ
            );
    
        // Check WSS exist
        // if (WSSheader.headerOk())
        // {
            mesh.readUpdate();
    
            Info<< "    Reading WSS" << endl;
            volVectorField WSS(WSSheader, mesh);
    
    
            forAll(transWSS.boundaryField(), patchi)
            {
                transWSS.boundaryFieldRef()[patchi] +=
                    mag(
                            WSS.boundaryField()[patchi]
                            &
                            (
                             normalVector.boundaryField()[patchi]
                             ^
                             (
                              TAWSS.boundaryField()[patchi]
                              / (mag(TAWSS.boundaryField()[patchi])+1e-64)
                             )
                            )
                       );
            }

            nfield++;
        }
        // else
        // {
        //     Info<< "    No WSS" << endl;
        // }
    
    
    Info<< "Writing transWSS" << endl;
    
    // devide by the number of added fields
    if(nfield>0){
        Info<< "number of fields added: "<< nfield << endl;
    
        forAll(transWSS.boundaryField(), patchi)
        {
            transWSS.boundaryFieldRef()[patchi] /= nfield;
        }
    }
    
    forAll(timeDirs, timeI)
    {
    transWSS.write();
    TAWSS.write();
    TAWSSMag.write();
    TAHct.write();
    OSI.write();
    RRT.write();
    normalVector.write();
    }

    Info<< "End" << endl;
    return 0;
}


// ************************************************************************* //
