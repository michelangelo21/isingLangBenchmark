#include <stdio.h>
#include <malloc.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif


uint32_t x, y, z, w;
 

uint32_t xorshift128(void){
    uint32_t t = x;
    t ^= t << 11;
    t ^= t >> 8;
    x = y; y = z; z = w;
    w ^= w >> 19;
    w ^= t;
    return w;
}

typedef struct spin_array{
	int16_t* array;		//pointer to data
	uint16_t l_size;		//side length
	uint16_t t_size;	//total size
	uint16_t* ni;
	uint16_t* pi;
} spin_array_t;

void SpinArrayConstructor(spin_array_t* sa, uint8_t l){
	sa->l_size = l; 
	sa->t_size = l*l;
	sa->array = (int16_t*) malloc(sizeof(int16_t) * sa->t_size);
	sa->ni = (uint16_t*) malloc(sizeof(uint16_t) * sa->l_size);
	sa->pi = (uint16_t*) malloc(sizeof(uint16_t) * sa->l_size);
	if( (sa->array == NULL) || (sa->ni == NULL) || (sa->pi == NULL) )
		printf("SpinArrayConstructor:\n(ERROR) Memory allocation error \n");
	
	for(uint16_t it = 0; it < sa->l_size; it++){
		sa->ni[it] = it-1;
		sa->pi[it] = it+1; 
	}
	sa->ni[0] = sa->l_size - 1;
	sa->pi[sa->l_size - 1] = 0;
}

void SpinArrayRandomInit(spin_array_t* sa){
	if(sa->array == NULL){
		printf("SpinArrayRandomInit:\n(ERROR) Memory pointer is null (array) \n");
	}else{
		for(uint16_t it = 0; it < sa->t_size; it++)
			sa->array[it] = (xorshift128() % 2) ? 1 : -1;
	}
}

inline int8_t SpinArrayGetSpinAt(spin_array_t* sa, uint16_t x, uint16_t y){
	return sa->array[ (y*sa->l_size) + x ];
}

void SpinArrayDisplay(spin_array_t* sa){
	for(uint16_t it_y = 0; it_y < sa->l_size; it_y++){
		for(uint16_t it_x = 0; it_x < sa->l_size; it_x++){
			printf(" %d ", (SpinArrayGetSpinAt(sa, it_x, it_y) > 0) ? 1 : 0);
		}
		printf("\n");
	}	
}

inline int8_t SpinArrayTrialEChange(spin_array_t* sa, uint16_t x, uint16_t y){
	return 2*SpinArrayGetSpinAt(sa,x,y)
	*(SpinArrayGetSpinAt(sa,sa->pi[x],y) + SpinArrayGetSpinAt(sa,sa->ni[x],y) + 
	SpinArrayGetSpinAt(sa,x,sa->ni[y]) + SpinArrayGetSpinAt(sa,x,sa->pi[y]));
}

inline uint8_t SpinArrayGetLSize(spin_array_t* sa){
	return sa->l_size;
}

inline uint16_t SpinArrayGetTSize(spin_array_t* sa){
	return sa->t_size;
}

inline void SpinArrayFlipSpinAt(spin_array_t* sa, uint16_t x, uint16_t y){
	sa->array[ (y*sa->l_size) + x ] = - sa->array[ (y*sa->l_size) + x ];
}

double* TableBoltzmanFactorConstructor()
{
	double* tab = (double*) malloc(sizeof(double)*9);
	if(tab == NULL)
		printf("TableBoltzmanFactorConstructor:\n(ERROR) Memory allocation error\n");
	return tab;
}

inline double* TableBoltzmanFactorInit(double* tab, double t_reduced){ 
	tab[0] = 1.00;  
	tab[4] = exp(-4.0/t_reduced);
	tab[8] = exp(-8.0/t_reduced);
	return tab; 
}	

void TableBoltzmanFactorDisplay(double* boltzman_table){
	printf(" dE = 0, Bf = %e \n", boltzman_table[0]);
	printf(" dE = 4, Bf = %e \n", boltzman_table[4]); 
	printf(" dE = 8, Bf = %e \n", boltzman_table[8]); 	
}

inline double RngUniformDistribution(void){
	return xorshift128()/(double)(4294967295U);
}

inline void SpinArrayCalculateME(spin_array_t* sa, double* magnetizatation_res, double* energy_res)
{
	double magnetization = 0;
	double energy = 0;
	uint16_t l_size = sa->l_size;
	uint16_t t_size = sa->t_size;
	for(uint16_t it_y = 0; it_y < l_size; it_y++){
		for(uint16_t it_x = 0; it_x < l_size; it_x++){
			magnetization += SpinArrayGetSpinAt(sa,it_x,it_y);
			energy += 0.5* (double)(SpinArrayGetSpinAt(sa,it_x,it_y)*
			(SpinArrayGetSpinAt(sa,sa->pi[it_x],it_y) + SpinArrayGetSpinAt(sa,sa->ni[it_x],it_y) 
			+ SpinArrayGetSpinAt(sa,it_x,sa->ni[it_y]) + SpinArrayGetSpinAt(sa,it_x,sa->pi[it_y])));
		}
	}
		
	*magnetizatation_res = fabs(magnetization/(double)t_size);
	*energy_res = energy;
}

void SpinArrayCalculateSpinOccurance(spin_array_t* sa){
	uint16_t n_spin_cnt = 0;
	uint16_t p_spin_cnt = 0;
	for(uint16_t  it = 0; it < sa->t_size; it++)
		if(sa->array[it] == 1) 
			p_spin_cnt++;
		else
			n_spin_cnt++; 
	printf("SpinArrayCalculateSpinOccurance:\n(INFO) Negative spin: %d Positive spin: %d\n", n_spin_cnt, p_spin_cnt);
}

void SpinArraySaveSpins(spin_array_t* sa, FILE* file)
{
	uint16_t l_size = SpinArrayGetLSize(sa);
	for(uint16_t y = 0; y < l_size; y++){
		for(uint16_t x = 0; x < l_size; x++){
			fprintf(file, "%d ", SpinArrayGetSpinAt(sa,x,y));
		}
		fprintf(file, "\n");
	}
	
}


void simulation(spin_array_t* model,
				double* boltzman_table,
				const uint64_t mcs_steps,
				const uint64_t mcs_observable_step,
				double temp,
				FILE* file){
	
	TableBoltzmanFactorInit(boltzman_table, temp);
	uint16_t it_x,it_y;
	uint16_t l_size = SpinArrayGetLSize(model);
	uint16_t t_size = SpinArrayGetTSize(model);
	double magnetization_sum = 0.0;
	double magnetization_sum_square = 0.0;
	double energy_sum = 0.0;
	double energy_sum_square = 0.0;
	uint64_t config_cnt = 0;
	double m,e,w;
	int8_t d_energy; 
		
	for(uint64_t it_mcs = 0; it_mcs < mcs_steps; it_mcs++)
	{		
		for(it_y = 0; it_y < l_size; it_y++){
			for(it_x = 0; it_x < l_size; it_x++){
				d_energy = SpinArrayTrialEChange(model, it_x, it_y);
				if(d_energy < 0){
					SpinArrayFlipSpinAt(model, it_x, it_y);		
				}else{
					w = min(1.0, boltzman_table[d_energy]);
					if(RngUniformDistribution() <= w){
						SpinArrayFlipSpinAt(model, it_x, it_y);
					}
				}
			}	
		}
		
	}	
		if(!(it_mcs % mcs_observable_step) && (it_mcs != 0))
		{
			config_cnt ++;
			SpinArrayCalculateME(model, &m, &e);
			magnetization_sum += m;
			magnetization_sum_square += (m*m);
			energy_sum += e;
			energy_sum_square += (e*e);
		}	
	}
	
	if(config_cnt != 0)
	{
		double average_magnetization = magnetization_sum/(double)config_cnt;
		
		double magnetic_susceptibility = ((double)t_size/(double)temp) *
		((magnetization_sum_square/((double)config_cnt)) - 
		( ((magnetization_sum/((double)config_cnt))) * ((magnetization_sum/((double)config_cnt)))));
	
		double thermal_capacity = (1.0/((double)t_size*temp*temp))*
		( (energy_sum_square/(double)config_cnt) - 
		((energy_sum/(double)config_cnt) * (energy_sum/(double)config_cnt)));
	
		fprintf(file, "%f %f %f %f \n", temp, average_magnetization, magnetic_susceptibility, thermal_capacity);
	}
	*/
}

int main(int argc, char* argv[]){
	x = rand();
	y = rand();
	z = rand();
	w = rand();
	const uint16_t l_size = atoi(argv[1]);
	
	const uint64_t mcs_steps = 230000;
	const uint64_t mcs_observable_step = 100;
	char file_name[256];
	strcpy(file_name, "...");
	strcat(file_name, argv[1]);
	strcat(file_name, ".txt");
	FILE *file = fopen(file_name, "w");

	spin_array_t model;
	double* boltzman_table = TableBoltzmanFactorConstructor();
	SpinArrayConstructor(&model, l_size);
	SpinArrayRandomInit(&model);
	SpinArrayCalculateSpinOccurance(&model);

	simulation(&model, boltzman_table, mcs_steps, mcs_observable_step, 8.0, file);
	SpinArraySaveSpins(&model, file);	
}