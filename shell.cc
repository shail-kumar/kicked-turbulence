#include "shell.h"
using namespace blitz;

// ------------------- Global constructs --------------------------------------
typedef complex<double> cdouble;

//++++++++++++++++++++++++ Control Parameters +++++++++++++++++++++++++++++++++
const double k_0 = 1.0/16;								
const double q = 0.5*(1+sqrt(5));	// Shell spacing
int N;	// # of Shells
double nu;						// Viscosity		
double dt;						// time step
int T;							// # of steps at an interval of dt
cdouble f_amp;					// Forcing amplitude
int n_f;						// Forced shell
int print_steps;				// Write the data in file at an interval of these many steps
double Dt;						// Time step of printing field
int IC_scheme;					// Initial condition scheme: 0 for random; 1 from  a file 	
int nic;						// Total number of initial conditions
double Omega;					// Rotation Strength

// =====================I/O Files =============================================
ifstream para_file;fstream initial_field_file;ofstream final_field_file;
ofstream out_file1;ofstream out_file2;ofstream out_file3;

// =====================main program ==========================================
int main()
{

	initial_field_file.open("initial_field.d");
	int attempt;
	cout<<"Which attempt?\t";
	cin>>attempt;

	stringstream final_u;
	final_u<<"../outdata/final_field_"<<attempt<<".d"; 
	final_field_file.open(final_u.str().c_str());
	// out_file2.open("../out_files/velocity_vs_t.d");

  /*================================Read input================================*/

	//--------Read parameters----------
	para_file.open("parameters.d");
	// parameter_reading(para_file, N, nu, dt, T, f_amp, n_f, print_steps, IC_scheme);
	parameter_reading(para_file);
	para_file.close();
	cout<<"# Number of Shells: "<<N<<", Viscosity: "<<nu<<", time step-size: "<<dt<<\
	", Total # of steps: "<<T<<", force amplitude: "<<f_amp<<", forced shellnumber: "<<n_f\
	<<", printing steps: "<<print_steps<<",	IC Scheme: "<<IC_scheme<<", # of ICs"<<nic<<\
	", Rot. Strength: "<<Omega<<endl;

	final_field_file<<"# Number of Shells: "<<N<<", Viscosity: "<<nu<<", time step-size: "<<dt<<\
	", Total # of steps: "<<T<<", force amplitude: "<<f_amp<<", forced shellnumber: "<<n_f\
	<<", printing steps: "<<print_steps<<",	IC Scheme: "<<IC_scheme<<", # of ICs: "<<nic<<\
	", Rot. Strength: "<<Omega<<endl;
	
		//-------Construction of shells----
	Array<cdouble,1> *u;
	u= new Array<cdouble,1>(N+4); 
	Array<double,1> *kn;	//,FortranArray<1>());
	kn=new Array<double,1>(N+2);
	shell_wavenumbers(kn);
	// out_file1<<setw(10)<<0<<"\t\t"; writeRealData(kn,out_file1);
	// cout<<"shell_wavenumbers"<<*kn<<endl;


	//------Force--------------------- 
	Array<cdouble,1> *force;
	force=new Array<cdouble,1>(N+4);
	// *force=0.0;
	// cout<<"force"<<*force<<endl;

	// ------------Create directories and file------------
	stringstream dir_field;
	dir_field<<"../outdata/velocity_field_"<<attempt;
	int dir = mkdir(dir_field.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if (-1==dir){cout<<"../outdata/velocity_field_"<<attempt<<" exists beforehand."<<endl;}
	else {cout<<"../outdata/velocity_field_"<<attempt<<" created."<<endl;}

	stringstream energy_evolution;
	energy_evolution<<"../outdata/energy_vs_t_"<<attempt<<".d"; 
	out_file3.open(energy_evolution.str().c_str());


	// ----------For output-------------
	Array<double,1> *mod_u;
	mod_u= new Array<double,1>(N+2);
	Dt = dt*print_steps; 
	int icn;							// icn --> IC number 


	clock_t begin = clock();
//***********************************************************************
//*							Loop over ICs 								*
//***********************************************************************
	for(icn=1;icn<=nic; icn++)
	{
		// ----------Create and open files-----------------------------
		stringstream u_field;
		u_field<<"../outdata/velocity_field_"<<attempt<<"/mod_u_"<<icn<<".d"; 
		out_file1.open(u_field.str().c_str());
		// Write wavenumbers for a given IC
		out_file1<<setw(10)<<0<<"\t\t"; writeRealData(kn,out_file1);
		//-------Read initial conditions---
		switch(IC_scheme)
		{
			case 0 :
				u_initial(u, kn); //By generating randonm numbers
				break;
			case 1 :
				initial_field(u, initial_field_file);			// Reading from a file
				break;
			default :
				initial_field(u, initial_field_file);			// Reading from a file
		}
		
		cout<<"# intial field\t"<<*u<<endl;

	  /*=================Processing and writing post-processed output=============*/
		*mod_u = 0.0;
		// double decay_rate = 0;
		double total_energy=0; 	
		for(int Dt_index = 1; Dt_index <= T; Dt_index++)
		{
			for (int dt_index = 1; dt_index <= print_steps; ++dt_index)
			{
				Time_evolution(u,force, kn);
			}
			// out_file2<<Dt_index*Dt<<"\t"<<*u<<endl;
			// 
			for (int j = 2; j < N+2; ++j)
			{
				(*mod_u)(j)=abs((*u)(j));
				total_energy += pow((*mod_u)(j),2)/2.0;
			}
			out_file1<<setw(10)<<Dt_index*Dt<<"\t\t"; writeRealData(mod_u,out_file1);
			out_file3<<icn<<"\t\t"<<((icn-1)*T+Dt_index)*Dt<<"\t\t"<<total_energy/Dt_index<<endl;
			cout<<icn<<"\t\t"<<Dt_index<<endl;
			// cout<<icn<<"\t\t"<<Dt_index*Dt<<"\t\t"<<total_energy/Dt_index<<"\n";
			// double forcing = decay_rate/((*u)(n_f+2).real() + (*u)(n_f+2).imag());

			if (isnan(abs((*u)(2))))
			{
				cout<<"Termination step: "<<Dt_index<<"\n NUMERICAL OVERFLOW! Program exits."; exit(0);
			}
			// cout<<i*dt<<"\t"<<*u<<endl;	
		}
		out_file1.close();
		final_field_file<<icn<<"\t\t"<<*u<<endl;
	}
  /*==============Processing and post processing ends here====================*/

	clock_t end = clock(); 
	cout<<"# Computation time = "<<double(end - begin)/CLOCKS_PER_SEC/60<<" minutes"<<endl;
	// initial_field_file<<"final_field u\t"<<*u<<endl;
	final_field_file.close();
	initial_field_file.close();
	// out_file2.close();
	out_file3.close();
	cout<<"# DONE!";	

	return 0;
}
/*===============================Main Program ends here====================*/