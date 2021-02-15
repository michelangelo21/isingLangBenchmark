//Jacek
#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>


double calculate_v(double x)
{
	return ( 0.5*pow(x,2.0) + 0.1*pow(x,3.0) );
}

const double calculate_step_size(const double start_step, const double stop_step, const uint64_t n)
{
	return (stop_step - start_step)/(n-1); 
}


double rng(void)
{
	return (double)rand()/RAND_MAX;
}


uint64_t rng_range(uint64_t max)
{
	return (uint64_t)rand()%(max+1);
}

double calculate_dt(double* phi_arr, double dphi_trial, double dphi_trial_sqr, double step_size, uint64_t idx, uint64_t n)
{
	double dt;
	if( (idx != 0) && (idx != n-1))
		dt = (dphi_trial_sqr - dphi_trial*(phi_arr[idx+1] + phi_arr[idx-1]))/pow(step_size, 2.0);
	if(idx == 0)
		dt = (dphi_trial_sqr - dphi_trial*(phi_arr[idx+1]))/pow(step_size, 2.0);
	if(idx == n-1)
		dt = (dphi_trial_sqr - dphi_trial*(phi_arr[idx-1]))/pow(step_size, 2.0);
	return dt;
}

int main(void)
{
	//Simulation parameters
	const double start_step = -3.0;
	const double stop_step = 3.0;
	const uint64_t n = 101;
	const double c = 0.5;
	const double dphi = 0.1;
	const uint64_t mc_step = 100000;
	//-----------------------------
	
	const double step_size = calculate_step_size(start_step, stop_step, n);	
	printf("step_size: %f\n", step_size);
	double* phi_arr = malloc(sizeof(double)*n);
	FILE *fp_energy = fopen("C:\\Users\\jatsekku\\Desktop\\MC2\\MC2_energy.txt", "w");
	FILE *fp_fi = fopen("C:\\Users\\jatsekku\\Desktop\\MC2\\MC2_f.txt", "w");
	double U = 0;
	double T = 0;
	double denominator = 0;
	
	uint64_t it;
	double x;
	for(it = 0, x = start_step; it < n; it++, x+= step_size) {
		phi_arr[it] = c;
		U += pow(phi_arr[it],2.0)*calculate_v(x);
		denominator += pow(phi_arr[it],2.0);
	}
	
	printf("U_0: %f\n", U);
	printf("Denominator_0: %f\n", denominator);
	T += phi_arr[0] * (2*phi_arr[0] - phi_arr[1]);
	T += phi_arr[n-1] * (2*phi_arr[n-1] - phi_arr[n-2]);
	T = 0.5*T/pow(step_size, 2.0);
	printf("T_0: %f\n", T);
	
	double numerator = U + T;
	printf("Numerator_0: %f\n", numerator);
	double E_start = numerator/denominator;	
	printf("E_start: %f \n", E_start);
	
	uint64_t x_idx;
	double phi_trial;
	double dphi_trial;
	double dU,dT;
	double E_curr;
	double dphi_trial_sqr;
	
	
	for(uint64_t its = 0; its < mc_step; its++ )
	{
		for(uint64_t itp = 0; itp < n; itp++)
		{
			x_idx = rng_range(n);
			phi_trial = phi_arr[x_idx] + (rng()-0.5)*dphi;
			
			dphi_trial = phi_trial - phi_arr[x_idx];
			dphi_trial_sqr = pow(phi_trial,2.0) - pow(phi_arr[x_idx], 2.0);
			
		
			dT = calculate_dt(phi_arr, dphi_trial, dphi_trial_sqr, step_size, x_idx, n);
			dU = dphi_trial_sqr*calculate_v(start_step + x_idx*step_size);
					
			E_curr = (numerator + dU + dT)/(denominator +  dphi_trial_sqr);
			
			if(E_curr < E_start)
			{
				phi_arr[x_idx] += dphi_trial;
				E_start = E_curr;
				numerator = numerator + dU + dT;
				denominator = denominator + dphi_trial_sqr;
			}
		}
		fprintf(fp_energy, "%f %d\n", E_start, its);	
	}
	
	double denominator2 = 0;
	
	for(int itn = 0; itn < n; itn++)
		denominator2 += pow(phi_arr[itn], 2.0) + pow(phi_arr[itn+1],2)/2 * step_size; 
	
	printf("denominator2 %f \n", denominator2);
	
	for(int itn = 0; itn < n; itn++)
		phi_arr[itn] = phi_arr[itn]/pow(denominator2,0.5);
	
	double integral_test = 0;
	for(int iti = 0; iti < n; iti++)
		integral_test += pow(phi_arr[iti]*step_size,2.0);
	
	printf("integral_test %f \n", integral_test);

	for(int itn = 0; itn < n; itn++)
		fprintf(fp_fi,"%f %f\n", phi_arr[itn], start_step + itn*step_size);
	
}