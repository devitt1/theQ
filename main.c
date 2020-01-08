//
//  main.c
//  QuBOX
//
//  Created by Simon Devitt on 2/8/19.
//  Copyright © 2019 Simon Devitt. All rights reserved.
//

#include <stdlib.h>
#include <time.h>
#include "stdio.h"
#include <math.h>
#include "sim.h"
#include "norm.h"
#include <string.h>
//#include <stdint.h>

void bv(int qs, FILE *out){
	int i;
	printf("Qubits = %d\n", qs);	
	printf("Qubits = %d\n", qs);	
	printf("Size of double: %ld bytes\n",sizeof(ds_Complex));
	printf("Size of double: %ld bytes\n",sizeof(double));
	//printf("SIZE_MAX       = %zu\n", SIZE_MAX);
        ds_Register ds_reg;
	double err = 0.01;
	ds_initialize_simulator(0);
	ds_reg = ds_create_register(qs, err, 0);
	ds_set_state(ds_reg, 0, 1, 0);
	
	for(i=0;i<qs;i++){
		ds_Hadamard(ds_reg, i, 0);		
	}
	for(i=0;i<qs;i++){
		ds_Z(ds_reg, i, 1);		
	}
	for(i=0;i<qs;i++){
		ds_Hadamard(ds_reg, i, 2);		
	}
	ds_print(ds_reg);

	int *qubits_to_measure = (int *)calloc(qs, sizeof(int));

	for(int j=0;j<qs;j++){
		qubits_to_measure[j] = qs - j - 1;
	}

        double correct_result = ds_set_measure(ds_reg, qs, qubits_to_measure, (1 << qs) - 1);
	printf("result = %f\n", correct_result);
	ds_destroy_register(ds_reg);

}

int main(){
    
    int i,j,k,N, qubit;

    FILE *out=fopen("out.dat","w");
    clock_t start, end;
    double cpu_time_used;
	for(i=4; i<5;i++){	
	start = clock();
	int qs = i;
	bv(qs, out);
    /*
    for(i=1; i<=30; i++){
        printf("Qubits = %d\n", i);
        for(j=1; j<=100; j++){
        
            start = clock();
            initialize_simulator();
            reg = create_register(i, 0, 0, out);
            set_state(reg, 0, 1, 0);
    
            for(k=0; k<j; k++){
				//void yrot(Register reg, int q, double theta, int time)
                //yrot(reg,rand() % i,Pi*uniform(),0);
				//void Hadamard(Register reg, int q, int time)

                yrot(reg,rand() % i,Pi*uniform(),0);
            }
            
            destroy_register(reg);
    
            end = clock();
            cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
            fprintf(out, "%d, %d, %f\n", i, j, cpu_time_used);
        }
    }
	*/
    
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	fprintf(out, "%d, %f\n", qs, cpu_time_used);
    }	
    fclose(out);
    return 0;
}

;
/*-----------------------end-------------------------------*/

