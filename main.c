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
#include "Simulator/sim.h"
#include "Simulator/norm.h"
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

void ibm_u3_gate_test(int qs, FILE *out){
	printf("Qubits = %d\n", qs);	
	printf("Qubits = %d\n", qs);	
	printf("Size of double: %ld bytes\n",sizeof(ds_Complex));
	printf("Size of double: %ld bytes\n",sizeof(double));
	//printf("SIZE_MAX       = %zu\n", SIZE_MAX);
	double err = 0.0;
        ds_Register ds_reg;

	ds_initialize_simulator(0);
	ds_reg = ds_create_register(qs, err, 0);
	ds_set_state(ds_reg, 0, 1, 0);

	//ds_u2(ds_reg.state+1, ds_reg.state+0, 0, ds_Pi);
	//ds_u2_gate(ds_reg.state+1, ds_reg.state+0, 0, ds_Pi);
	

	//ds_Hadamard(ds_reg, 0, 1);

	//ds_U2(ds_reg, 0, 1, 0, ds_Pi);
	ds_Hadamard(ds_reg, 0, 1);
	ds_U1(ds_reg, 0, 1, ds_Pi/2);

	ds_print(ds_reg);
        ds_destroy_register(ds_reg);
}

int main(){
    
    int qs = 1;

    FILE *out=fopen("out.dat","w");
    clock_t start, end;
    double cpu_time_used;

    start = clock();
    ibm_u3_gate_test(1,out);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    fprintf(out, "%d, %f\n", qs, cpu_time_used);
    fclose(out);
    return 0;
}

;
/*-----------------------end-------------------------------*/

