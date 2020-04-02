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


int main(){
    
    int qs = 1;
    FILE *out=fopen("out.dat","w");
    clock_t start, end;
    double cpu_time_used;
    start = clock();
	
	
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    fprintf(out, "%d, %f\n", qs, cpu_time_used);
    fclose(out);
    return 0;
}

;
/*-----------------------end-------------------------------*/

