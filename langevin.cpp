#include "odesolver.h"
#include <cmath>
#include <iostream>

using namespace std;

/*
 * yi[0->4] = x1,y1,x2,y2,t
 * void fucntion of distance, xy positions of bead 1 and 2
 * and conitains information about potential where
 * depth in kbt (U0)
 * micro meter position of min (r0)
 * passes forces into force member function
 */
void Langevin_Problem::force_calc(double *y)
{
	double dist = sqrt( (y[0]-y[2])*(y[0]-y[2]) + 
                  (y[1]-y[3])*(y[1]-y[3]) );
	sep = dist;
	double pdepth= p.pdepth;
        double pmin = p.pmin;
        double rodivdist = pmin/dist;
	double uprime;
	double repulsion;

	/*
	 * choose repulsion type
	 */

	if( p.rep_type== 0){
	  repulsion= 0;
	}
	if( p.rep_type== 1){
	  repulsion= p.strength;
	}
	if( p.rep_type== 2){
	  repulsion= 2*p.strength*dist;
	}
	if( p.rep_type== 3 ){
	  repulsion= (-3*p.strength)/(dist*dist);
	}
	if( p.rep_type== 4){
	  repulsion= (-2*p.strength)/(dist*dist*dist);
	}

	//uprime = 5*dist*dist*dist*dist;
	// uprime = -(pdepth/pmin)*8*(pow(rodivdist,9)
	//         - pow(rodivdist,5));
        // uprime = 0;
	//if( dist < pmin ){
	//  uprime += -(pdepth/pmin)*8*(pow(rodivdist,9)
	//			     - pow(rodivdist,5));
	//}
	// double uprime = 0.0;
        // uprime += 0.25*dist;
	// double uprime = 0.25*dist;
	uprime= -(pdepth/pmin)*(2*p.power)*(pow(rodivdist,((2*p.power)+1))
					    -pow(rodivdist,((p.power)+1)))
	  + repulsion;
	force[0] = -(y[0]-y[2])/dist*uprime;
	force[1] = -(y[1]-y[3])/dist*uprime;
	force[2] = (y[0]-y[2])/dist*uprime;
	force[3] = (y[1]-y[3])/dist*uprime;
	return;
};

/*
 * yi[0->4] = x1,y1,x2,y2,t
 * second order Runge-Kutta to find new position
 * step_choice and following if statements choose appropriate
 * variable timesteps for each run under a uniform distribution.
 * Once conditions are met, prints to standard out positions and 
 * time while passing solved positions into yi[i](i:0->4) member
 */
void Langevin_Problem::step_solution ()
{  

    double m1 = p.m1;
    double m2 = p.m2;
    double gamma1 = p.gamma1;
    double gamma2 = p.gamma2;
    double kbt = p.kbt;
    double delta_t = p.delta_t;
    double t_frame;

    double step_choice = gsl_rng_uniform(r);
    if(step_choice >= 0 && step_choice <=.01){
      t_frame = .2;}
    if(step_choice > .01 && step_choice <=.02){
      t_frame =.18;}
    if(step_choice > .02 && step_choice <= .05){
      t_frame =.19;}
    if(step_choice > .05 && step_choice <= .5){
      t_frame = .09;}
    if(step_choice > .5 && step_choice<= 1){
      t_frame = .1;}

    int threshold = static_cast<int>(t_frame/delta_t);

	for(int i=0; i<threshold; ++i){
		

		for( int j=0; j<4; ++j){
			p.Rand[j] = gsl_ran_gaussian( r, 1./sqrt( delta_t));
		}
		int broken;
		force_calc(yi);

		k1[0]=(kbt/(gamma1*m1))*force[0]
		          +sqrt(2*kbt/(gamma1*m1))*p.Rand[0];
		k1[1]=(kbt/(gamma1*m1))*force[1]
			  +sqrt(2*kbt/(gamma1*m1))*p.Rand[1];
		k1[2]=(kbt/(gamma2*m2))*force[2]
			  +sqrt(2*kbt/(gamma2*m2))*p.Rand[2];
		k1[3]=(kbt/(gamma2*m2))*force[3]
			  +sqrt(2*kbt/(gamma2*m2))*p.Rand[3];
 
		ym[0] = yi[0] + k1[0]*(0.5*delta_t); 
		ym[1] = yi[1] + k1[1]*(0.5*delta_t);
		ym[2] = yi[2] + k1[2]*(0.5*delta_t);
		ym[3] = yi[3] + k1[3]*(0.5*delta_t);

		if(sep <= 2){
		  broken= 1;
		}
		else{
		  broken=0;
		}

		force_calc(ym);

		k2[0]=(kbt/(gamma1*m1))*force[0]
			  +sqrt(2*kbt/(gamma1*m1))*p.Rand[0];
		k2[1]=(kbt/(gamma1*m1))*force[1]
			  +sqrt(2*kbt/(gamma1*m1))*p.Rand[1];
		k2[2]=(kbt/(gamma2*m2))*force[2]
			  +sqrt(2*kbt/(gamma2*m2))*p.Rand[2];
		k2[3]=(kbt/(gamma2*m2))*force[3]
			  +sqrt(2*kbt/(gamma2*m2))*p.Rand[3];
		
	
		yn[0]=yi[0]+delta_t*0.5*(k1[0]+k2[0]);
		yn[1]=yi[1]+delta_t*0.5*(k1[1]+k2[1]);
		yn[2]=yi[2]+delta_t*0.5*(k1[2]+k2[2]);
		yn[3]=yi[3]+delta_t*0.5*(k1[3]+k2[3]);
	       
		if(sep <= 2){
		  //double shift2 = yn[2]+sqrt(4-((yn[1]-yn[3])*(yn[1]-yn[3])));
		  //yn[0]=shift2;
		  broken=1;
		}
		else{
		  broken=0;
		}
		//if(broken==0){
		yn[4]=yn[4]+delta_t;	
	        for(int k=0; k<5; ++k){
		    yi[k]=yn[k];
	        }
		  //}
		  //	else{
		  // i=i-1;
		  // }
		  
	}
	cout << yi[4] << " "
		 << yi[0] << " "
		 << yi[1] << " "
		 << yi[2] << " "
		 << yi[3] << " " << endl;
	
	return;
}

/*
 * Calls step_solution to drive problem
 */
void Langevin_Problem::solve()
{
	do{
		step_solution ();
	}while(yi[4] <= 180);
	return;
}
