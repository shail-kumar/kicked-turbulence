#include "shell.h"
using namespace blitz;

typedef complex<double> cdouble;
const long double PI = 3.141592653589793238L;


/*void parameter_reading(ifstream& para_file, int &no_of_shells, double &nu, double &dt, int &T,\
 cdouble &f_amp, int &n_f, int &print_steps, int &IC_scheme)*/
void parameter_reading(ifstream& para_file)
{
	string useless_line;
	if (para_file.is_open())
	{
		getline(para_file, useless_line);
		para_file>>N>>nu>>dt>>T>>f_amp>>n_f>>print_steps>>IC_scheme>>nic>>Omega;
		// para_file>>no_of_shells>>nu>>dt>>T>>f_amp>>n_f>>print_steps>>IC_scheme;
		// cout<<no_of_shells<<"		"<<nu<<"		"<<dt<<"		"<<T<<"		"<<endl;
	}
	else
	{
		cout<<"Parameter file is not open"<<endl;
	}

}

void shell_wavenumbers(Array<double,1>* k_n)
{
	for (int i = 0; i < N+2; ++i)
	{
		(*k_n)(i)=k_0*pow(q,i-1);			//i=2 corressponds to 1st shell
	}
	// cout<<1<<(*k_n)<<endl;
	// return (*k_n);
}

void u_initial(Array<cdouble,1>* u0, Array<double,1>* k_n)
{
	srand(time(0));
	for(int i=2;i<N+2;i++){
		// u(i)=1.0*rand()/(RAND_MAX +1.0);
		// (*u0)(i) = pow((*k_n)(i),2);
		double theta = 2*PI*rand()/(RAND_MAX + 1.0);
		if (i<5)
		{
			(*u0)(i) = cdouble (cos (theta), sin (theta))*sqrt((*k_n)(i));
		}
		else
		{
			(*u0)(i) = cdouble (cos (theta), sin (theta))*sqrt((*k_n)(i))*exp(-pow((*k_n)(i),2));
		}

	}

	// return (*u0);	
}

void initial_field(Array<cdouble,1>* u0, fstream& initial_field_file)
{
	initial_field_file >> (*u0);
	cout  << "Reading of field configurations ended successfully" << endl; 
}