double RHO_0;
double dt;
int N_OUTLETS;
DynamicList<string> patch_names(10); // 10 has been set as the maximum limit of outlets that are expected

/* Windkessel Structure Definition */
typedef struct
{
    double Q_4;
    double Q_3;
    double Q_2;
    double Q_1;
    double Q_0;
    double P_3;
    double P_2;
    double P_1;
    double P_0;
    double P1;
    int id;             /* Windkessel element id */
    double R;       /* Resistance */
    double C;       /* Compliance */
    double Z;       /* Impedance */
    bool physUnits; /* Physiological or SI units for R, C, Z */
    int order;
    double time;    // store time in windkessel struct to reset values
                    // should the time step be run multiple times
                    // e.g. for FSI simulations

} WindKessel;


WindKessel *wk;

void initialise(const dictionary& windkesselProperties)
{

    /* Initialising Windkessel object */

    wk = (WindKessel*)malloc(N_OUTLETS*sizeof(WindKessel));

    /* Retrieving values from windkessel dictionary*/

    const wordList outletNames(windkesselProperties.toc());

    forAll(outletNames, item)
    {
        const word& outletName = outletNames[item];
        Info << "Evaluating properties for " << outletName << endl;

        const dictionary& subDict = windkesselProperties.subDict(outletName);

        scalar Z = readScalar(subDict.lookup("Z"));
        scalar C = readScalar(subDict.lookup("C"));
        scalar R = readScalar(subDict.lookup("R"));
        bool physUnits = readBool(subDict.lookup("physiologicalUnits"));
        scalar real_index = readScalar(subDict.lookup("outIndex"));
        scalar Qpre1 = readScalar(subDict.lookup("Flowrate_oneStepBefore"));
        scalar Qpre2 = readScalar(subDict.lookup("Flowrate_twoStepBefore"));
        scalar Qpre3 = readScalar(subDict.lookup("Flowrate_threeStepBefore"));
        scalar Ppre2 = readScalar(subDict.lookup("Pressure_twoStepBefore"));
        scalar Ppre1 = readScalar(subDict.lookup("Pressure_oneStepBefore"));
        scalar Pnow = readScalar(subDict.lookup("Pressure_start"));
        scalar FDMorder = readScalar(subDict.lookup("FDM_order"));

        int out_index = real_index;
        int order = FDMorder;


        /* e.g., last saved time: 80.00, running from 80.01 ---> 0: 80.01, -1: 80.00, -2: 79.99, -3: 79.98*/
        if (physUnits == true)
            {
            Z = Z * 133.322365e6;
            C = C * 1e-6 / 133.322365;
            R = R * 133.322365e6;
            Info << "Converting Z, C, R from physioligical units to SI units!" <<endl;;
            }

        Info << "Z, C, R and index are " << Z << ", " << C << ", " << R << ", " << out_index << "." <<endl;

        wk[out_index].Z = Z;
        wk[out_index].C = C;
        wk[out_index].R = R;
        wk[out_index].physUnits = physUnits;
        wk[out_index].id = out_index;
        wk[out_index].Q_3 = Qpre3;
        wk[out_index].Q_2 = Qpre2;
        wk[out_index].Q_1 = Qpre1;
        wk[out_index].Q_0 = 0;      // updated from Windkessel
        wk[out_index].P_1 = Ppre1;
        wk[out_index].P_2 = Ppre2;
        wk[out_index].P_0 = Pnow;
        wk[out_index].P1 = 0;       // updated from Windkessel
        wk[out_index].order = order;
        wk[out_index].time = 0.0; // can I store the time step?
    }

}

double calculate_pressure_veryFirst(int i, fvMesh & mesh, volScalarField& p)
{
    double outPressure = 0.0;
    // Pressure calculation for each outlet

    label outletPatch = mesh.boundaryMesh().findPatchID(patch_names[i]);

    double psize = p.boundaryField()[outletPatch].size();
    Info << "psize " << psize << endl;

    if (outletPatch >= 0)
    {
        outPressure = 1060.0 * average(p.boundaryField()[outletPatch]);
    }
    reduce(outPressure, sumOp<scalar>());

    Info << "Pressure for " << patch_names[i] << " :  " << outPressure << endl;

    return outPressure;
}


double calculate_flowrate(int i, fvMesh & mesh, surfaceScalarField& phi)
{
    scalar outflow = 0.0;

    label outletPatch = mesh.boundaryMesh().findPatchID(patch_names[i]);
    if (outletPatch >= 0)
    {
        outflow = sum(phi.boundaryField()[outletPatch]);
    }

    reduce(outflow, sumOp<scalar>());
    //Info << "Flowrate for " << patch_names[i] << " :  " << outflow << endl;

    return outflow;
}

void execute_at_end(fvMesh& mesh, surfaceScalarField& phi, scalarIOList& store, volScalarField& p, dictionary& windkesselProperties, double& time)
{

    bool repeatTS;

    // check if time step is repeated:
    if (wk[0].time == time)
    {
        Info << "time step not advanced, reverting WK properties" << endl;

        repeatTS = true;
    }
    else
    {
        Info << "time step advanced" << endl;
        repeatTS = false;
    }
    
    int i;

    for (i=0;i<N_OUTLETS;i++)
    {


    // Revert back to the previous time step if the time step is repeated (e.g. in FSI calcs using preCICE)
    if (repeatTS)
    { 
        Info << "Time: " << time << ", Stored time: " << wk[0].time << endl;
        wk[i].P_0 = wk[i].P_1;
        wk[i].P_2 = wk[i].P_1;
        wk[i].Q_0 = wk[i].Q_1;
        wk[i].Q_1 = wk[i].Q_2;
        wk[i].Q_2 = wk[i].Q_3;
        wk[i].Q_3 = wk[i].Q_4;

        Info << patch_names[i] << ":" << endl;
        Info << "P(n-1) = " << wk[i].P_0 << endl; 
    }

    wk[i].Q_0 = calculate_flowrate(i, mesh, phi);                   //current Q

    // Windkessel third order
    double Q_source = 0.0;
    double Pgrad_part = 0.0;
    double Pdenom = 0.0;

    /* Identifying the order of finite difference method */
    int flag = wk[i].order;
    /* Discretisation depending on the selected order */
    /* last saved: 80.00, running from 80.01 ---> 0: 80.01, -1: 80.00, -2: 79.99, -3: 79.98*/
    switch (flag)
    {
        case 1:     //1st order
            Q_source = (wk[i].Q_0 / wk[i].C)*(1 + wk[i].Z / wk[i].R) + (wk[i].Z / dt)*(wk[i].Q_0 - wk[i].Q_1);
            Pgrad_part = - wk[i].P_0 / dt;
            Pdenom = 1 / dt + 1 / (wk[i].R*wk[i].C);
            break;
        case 2:     //2nd order
            Q_source = (wk[i].Q_0 / wk[i].C)*(1 + wk[i].Z / wk[i].R) + (wk[i].Z / dt)*( 3/2 * wk[i].Q_0 - 2 * wk[i].Q_1 + 1/2 * wk[i].Q_2 );
            Pgrad_part = -2 * wk[i].P_0 / dt + wk[i].P_1 / (2 * dt);
            Pdenom = 3 / (2*dt) + 1 / (wk[i].R*wk[i].C);
            break;
        case 3:     //3rd order
            Q_source = (wk[i].Q_0 / wk[i].C)*(1 + wk[i].Z / wk[i].R) + (wk[i].Z / dt)*(11 / 6 * wk[i].Q_0 - 3 * wk[i].Q_1 + 3 / 2 * wk[i].Q_2 - 1 / 3 * wk[i].Q_3);
            Pgrad_part = -3 * wk[i].P_0/dt + (3*wk[i].P_1)/(2*dt) - wk[i].P_2/(3*dt);
            Pdenom = 11 / (6 * dt) + 1 / (wk[i].R*wk[i].C);
            break;
        default: //1st order
            Q_source = (wk[i].Q_0 / wk[i].C)*(1 + wk[i].Z / wk[i].R) + (wk[i].Z / dt)*(wk[i].Q_0 - wk[i].Q_1);
            Pgrad_part = -wk[i].P_0 / dt;
            Pdenom = 1 / dt + 1 / (wk[i].R*wk[i].C);
    }
    // calculate pressure
    wk[i].P1 = (Q_source - Pgrad_part) / Pdenom;

    // Saving the pressure in a scalar array
    store[i] = wk[i].P1;

    //WRITE OUT TO constant/windkesselProperties
    dictionary& WK_dict = windkesselProperties.subDict(patch_names[i]);
    Info << patch_names[i] << ": " << endl;
    if (flag >= 3)
    {
        WK_dict.set("Flowrate_threeStepBefore", wk[i].Q_2);
    }
    if (flag >= 2)
    {
        WK_dict.set("Flowrate_twoStepBefore", wk[i].Q_1);
    }
      Info << "Flowrate    " << wk[i].Q_0 << ";" << endl;
        WK_dict.set("Flowrate_oneStepBefore", wk[i].Q_0);

    if (flag >= 3)
    {
        WK_dict.set("Pressure_twoStepBefore", wk[i].P_1);
    }
    if (flag >= 2)
    {
        WK_dict.set("Pressure_oneStepBefore", wk[i].P_0);
    }
        WK_dict.set("Pressure_start", wk[i].P1);

    /* update */
        wk[i].P_2 = wk[i].P_1;
        wk[i].P_1 = wk[i].P_0;
        wk[i].P_0 = wk[i].P1;
        wk[i].Q_3 = wk[i].Q_2;
        wk[i].Q_2 = wk[i].Q_1;
        wk[i].Q_1 = wk[i].Q_0;
        wk[i].time = time;

    Info << "P(n)   = " << wk[i].P1 << endl;

    
    }

}
