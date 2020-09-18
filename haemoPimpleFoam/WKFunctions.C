
double RHO_0;
double dt;
int N_OUTLETS;
//char* patch_names[] = {"OUTLET_ACA","OUTLET_MCA"};
DynamicList<string> patch_names(10); // 10 has been set as the maximum limit of outlets that are expected

/* Windkessel Structure Definition */
typedef struct {
	double Q_current;    	/* Current  time stepproximal flow rate */
	double Q_previous; 	  /* Previous time step proximal flow rate */
	double Q_previous2;  	/* 2 Previous time step proximal flow rate */
	double P_current;		  /* Current time stepproximal pressure */
	double P_previous;		  /* Previous time step proximal pressure */
	double P_previous2;	  /* 2 Previous time step proximal pressure */
	double Pout_current;		/* Current back-pressure */
	double Pout_previous;	/* Previous back-pressure */
	double Pout_previous2;  /* 2 Previous back-pressure */
	double Pc_current;		/* Current back-pressure */
	double Pc_previous;	/* Previous back-pressure */
	double Pc_previous2;  /* 2 Previous back-pressure */
	int id;			/* Windkessel element id */
	double R;		/* Resistance */
	double C;		/* Compliance */
	double Z;		/* Impedance */
} WindKessel;


WindKessel *wk;

void initialise(const dictionary& windkesselProperties)
{


 /* Initialising Windkessel object */

	int i;
	wk = (WindKessel*)malloc(N_OUTLETS*sizeof(WindKessel));
	for (i=0;i<N_OUTLETS;i++)
	{
		/* retrieve from Scheme instantiation : First instantiation variables */
		wk[i].Q_current 	= 0;
		wk[i].Q_previous 	= 0;
		wk[i].Q_previous2 	= 0;
		wk[i].P_current 	= 0; /* Proximal pressure */
		wk[i].P_previous 	= 0;
		wk[i].P_previous2 	= 0;
		wk[i].Pout_current 	= 0;  /* Venous pressure */
		wk[i].Pout_previous = 0;
		wk[i].Pout_previous2 = 0;
		wk[i].Pc_current 	= 0;	/* Intramural pressure */
		wk[i].Pc_previous 	= 0;
		wk[i].Pc_previous2 	= 0;
	}	



	/* Retrieving values from windkessel dictionary*/

 const wordList outletNames(windkesselProperties.toc());

 forAll(outletNames, item)
   {
     const word& outletName = outletNames[item];
     Info << "Evaluating properties for " << outletName << endl;

     const dictionary& subDict = windkesselProperties.subDict(outletName);

     scalar C = readScalar(subDict.lookup("C"));
     scalar R = readScalar(subDict.lookup("R"));
     scalar Z = readScalar(subDict.lookup("Z"));
     scalar real_index = readScalar(subDict.lookup("outIndex"));


     int out_index = real_index;
   
     Info << "C, R, Z and index are " << C << ", " << R << ", " << Z << ", " << out_index << "." <<endl; 

     wk[out_index].R 			= R;
     wk[out_index].C 			= C;
     wk[out_index].Z 			= Z;
     wk[out_index].id 			= out_index;
				
   }


 /* End of file*/

}


double first_order_derivative(double x, double xp)
{
	double derivative;

	derivative = (x - xp)/dt;

	return derivative;
}

double front_first_order_derivative()
{
	double derivative;

	derivative = 1.0/dt;

	return derivative;
}

double back_first_order_derivative(double xp)
{
	double derivative;

	derivative = -xp/dt;

	return derivative;
}

double second_order_derivative(double x, double xp,double xp2)
{
	double derivative;

	derivative = (1.5*x - 2.0*xp + 0.5*xp2)/dt;

	return derivative;
}

double front_second_order_derivative()
{
	double derivative;

	derivative = 1.5/dt;

	return derivative;
}

double back_second_order_derivative(double xp, double xp2)
{
	double derivative;

	derivative = (-2.0*xp + 0.5*xp2)/dt;

	return derivative;
}

double derivative(double x, double xp, double xp2)
{
	double d;

	// d = first_order_derivative(x,xp);

	d = second_order_derivative(x,xp,xp2);

	return d;
}

double front_derivative()
{
	double d;
	//d = front_first_order_derivative();

	d = front_second_order_derivative();

	return d;
}

double back_derivative(double xp, double xp2)
{
		double d;
		//d = back_first_order_derivative(xp);
		d = back_second_order_derivative(xp,xp2);

		return d;
}

double calculate_flow_rate_additional(int i, fvMesh & mesh, volVectorField & U)
{
	scalar flux = 0.0;

	word outlet_name = patch_names[i];

	//access boundary elements 
	const fvPatchList& patches = mesh.boundary(); 

	//loop over boundaries 
	forAll(patches, patchI) 
	{ 
	    const fvPatch& cPatch = patches[patchI]; 

	    //check boundary name 
	    if(cPatch.name() == outlet_name) 
	      { 

		//Access boundary face area vectors 
		const vectorField& inletAreaVectors = cPatch.Sf(); 

		//loop over all faces of selected 
		//boundary 
		forAll(cPatch, faceI) 
		{ 

		    //calculate boundary face mass flux 
                    scalar cFaceFlux = U.boundaryField()[patchI][faceI] & inletAreaVectors[faceI]; 

                     //sum face mass fluxes 
                     flux += cFaceFlux; 
		} 

	      } 

	  } 

	Info<< "FlowRate: " << flux << " m^3/s" << endl; 


	return flux;

}

double calculate_flow_rate(int i, fvMesh & mesh, surfaceScalarField & phi)
{
	scalar outflow = 0.0;

	label outletPatch = mesh.boundaryMesh().findPatchID(patch_names[i]);

	if (outletPatch >=0)
	{
	  outflow = sum(phi.boundaryField()[outletPatch]);
	}

	reduce(outflow, sumOp<scalar>());

	Info << "Flowrate for " << patch_names[i] << " :  " << outflow << endl; 

	return outflow;

}

void Wk_pressure_update(int i, double rho, fvMesh & mesh, surfaceScalarField & phi, scalarIOList & store)
{
  
	scalar p,dpc,dpq;

	//wk[i].Q_current = calculate_flow_rate(i,mesh,U);
	wk[i].Q_current = calculate_flow_rate(i,mesh,phi);
	dpc = derivative(wk[i].Pc_current,wk[i].Pc_previous,wk[i].Pc_previous2);
	dpq = derivative(wk[i].Q_current,wk[i].Q_previous,wk[i].Q_previous2);

	p = wk[i].Q_current
	- wk[i].C*back_derivative(wk[i].P_previous,wk[i].P_previous2)
	+ wk[i].Z*(wk[i].C*dpq+wk[i].Q_current/wk[i].R)
	+ wk[i].Pout_current/wk[i].R+wk[i].C*dpc;

	scalar P_current = p/(1.0/wk[i].R+wk[i].C*front_derivative());

	wk[i].P_current = P_current;
    

	/*Saving the pressure in a scalar array*/
	store[i]=P_current;

	Info<< "Pressure: " << P_current << " Pa" << endl; 

}

void execute_at_end(fvMesh & mesh, surfaceScalarField & phi, scalarIOList & store)
{

  int i;

  for (i=0;i<N_OUTLETS;i++)
    {

      /* Save previous states */
      wk[i].P_previous2 = wk[i].P_previous;
      wk[i].Q_previous2 = wk[i].Q_previous;
      wk[i].Pout_previous2 = wk[i].Pout_previous;
      wk[i].Pc_previous2 = wk[i].Pc_previous;

      wk[i].P_previous = wk[i].P_current;
      wk[i].Q_previous = wk[i].Q_current;
      wk[i].Pout_previous = wk[i].Pout_current;
      wk[i].Pc_previous = wk[i].Pc_current;

      /*Update WindKessel values*/

      Wk_pressure_update(i, RHO_0, mesh, phi,store);

    }

}
