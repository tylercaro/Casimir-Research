#include <iostream>
#include <cmath>
#include <ctime>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "odesolver.h"

using namespace std;

/*
 * Constructor initializing current position as initial position
 * starting time set to 0, y[8]
 * and creates a generator chosen by the 
 * environment variable GSL_RNG_TYPE
 */
Langevin_Problem::Langevin_Problem(langevin_params params, double *y)
{
	p = params;
	for( int i=0; i < 5; ++i)
	{
		yi[i]=y[i];
	}
	
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	gsl_rng_set( r, static_cast<unsigned long int>(time(NULL)) );
}

/* Deconstructor freeing up allocated space to generate the
 * random numbers
 */
Langevin_Problem::~Langevin_Problem()
{
	gsl_rng_free( r);
}
/*
 * argv[run,depth,strength,power,repulsion]
 */
int main(int argc, char* argv[])
{

  langevin_params params;
  params.m1 = 1.3791591749259193e-15;
  params.m2 = 3.723729772299982e-14;
  params.gamma1 = 2.310186846116964e7;
  params.gamma2 = 241806.9967755277;
  params.kbt = 4.14195e-9;
  params.pdepth = strtod(argv[1],NULL);
  params.pmin = 3;
  params.delta_t = .0001;
  params.strength= strtod(argv[2],NULL);
  params.power= strtod(argv[3],NULL);
  params.rep_type= strtod(argv[4],NULL);

  double y[5];
  y[0] =3;
  y[1] =0;
  y[2] =0;
  y[3] =0;
  y[4] =0;
   
  Langevin_Problem langevin_problem( params, y);
  langevin_problem.solve();
  return 0;
 }
