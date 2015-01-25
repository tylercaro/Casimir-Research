#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#ifndef ODESOLVER
#define ODESOLVER

struct langevin_params
{
    double m1;
    double m2;
    double gamma1;
    double gamma2;
    double kbt;
    double pdepth;
    double pmin;
    double delta_t;
    double Rand[4];
    double strength;
    double power;
    double rep_type;
};

class Langevin_Problem
{
	private:
	 langevin_params p;
	 double yi[5];//current step
	 double ym[4];//half step
	 double k1[4];//2nd order r-k
	 double k2[4];//2nd order r-k
	 double yn[5];//full step
	 double force[4];
	 double sep;
	 
	 /* 
	  * gsl rng objects
	  */
	 const gsl_rng_type * T;
	 gsl_rng * r;
	   
	  /*
	   * force function
	   */
	 void force_calc(double *data);
        public:
	  Langevin_Problem( langevin_params params, double *y);
	  ~Langevin_Problem();
	  
	  void solve();
	  void step_solution();
  };
#endif 
