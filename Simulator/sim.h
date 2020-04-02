#ifndef SIM_H
#define SIM_H

#include <stdio.h>

typedef struct {
   double x,y ;
} ds_Complex;

typedef struct {
   ds_Complex *state;
   int *steps;
   int nq, nc;
   int *n_errors;
   double err, sigma;
} ds_Register;

extern double ds_Pi, ds_Pio2, ds_Pio4, ds_Pio8, ds_root2_2;
extern ds_Register ds_reg;

ds_Register ds_create_register(int nq_L, double err_L, double sigma_L);
void ds_destroy_register(ds_Register reg);
void ds_equate_registers(ds_Register reg1, ds_Register reg2);
void ds_initialize_simulator(long ds_seed);
void ds_clearreg(ds_Register reg);
void ds_set_state(ds_Register reg, int n, double x, double y);
int ds_query_state(ds_Register reg, int n, double tol);
void ds_print(ds_Register reg);
void ds_update(ds_Register reg, int q1, int q2);
void ds_global_update(ds_Register reg);

ds_Complex ds_eitheta(double theta);
ds_Complex ds_add(ds_Complex z1, ds_Complex z2);
ds_Complex ds_multiply(ds_Complex z1, ds_Complex z2);
ds_Complex ds_zstarz(ds_Complex z1, ds_Complex z2);
double ds_modsq(ds_Complex z);
double ds_inner_product(ds_Register reg1, ds_Register reg2);

void ds_esigx(ds_Complex *zPtr1, ds_Complex *zPtr2, double theta);
void ds_esigy(ds_Complex *zPtr1, ds_Complex *zPtr2, double theta);
void ds_esigz(ds_Complex *zPtr1, ds_Complex *zPtr2, double theta);
void ds_unitary(ds_Complex *zPtr1, ds_Complex *zPtr2,
             double alpha, double beta, double theta);

void ds_one_qubit_indices(int n, int q, int *iPtr, int *jPtr);
void ds_controlled_indices(int n, int q1, int q2, int *iPtr, int *jPtr);
void ds_swap_indices(int n, int q1, int q2, int *iPtr, int *jPtr);
void ds_two_qubit_indices(int n, int q1, int q2, int *i1Ptr, int *i2Ptr, int *i3Ptr, int*i4Ptr);

void ds_xrot(ds_Register reg, int q, double theta, int time);
void ds_yrot(ds_Register reg, int q, double theta, int time);
void ds_zrot(ds_Register reg, int q, double theta, int time);

void ds_X(ds_Register reg, int q, int time);
void ds_Z(ds_Register reg, int q, int time);
void ds_XZ(ds_Register reg, int q, int time);

void ds_cnot(ds_Register reg, int qcont, int qtarg, int time);
void ds_swap(ds_Register reg, int q1, int q2, int time);
void ds_lerr(ds_Register reg, int q, int time);
void ds_cphase(ds_Register reg, int qcont, int qtarg, double theta, int time);
void ds_Hadamard(ds_Register reg, int q, int time);

int ds_measure(ds_Register reg, int nq2m, int *lq2m);

double ds_set_measure(ds_Register reg, int nq2m, int *lq2m, int val);

#endif
