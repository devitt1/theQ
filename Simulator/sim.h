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
void ds_change_err(ds_Register *regPtr, double err_L);
void ds_change_sigma(ds_Register *regPtr, double sigma_L);
void ds_update(ds_Register reg, int q1, int q2);
void ds_global_update(ds_Register reg);

void dosmth(void);

ds_Complex ds_eitheta(double theta);
ds_Complex ds_add(ds_Complex z1, ds_Complex z2);
ds_Complex ds_multiply(ds_Complex z1, ds_Complex z2);
ds_Complex ds_zstarz(ds_Complex z1, ds_Complex z2);
double ds_modsq(ds_Complex z);
double ds_inner_product(ds_Register reg1, ds_Register reg2);

void ds_esigi(ds_Complex *zPtr1, ds_Complex *zPtr2, double theta);
void ds_esigx(ds_Complex *zPtr1, ds_Complex *zPtr2, double theta);
void ds_esigy(ds_Complex *zPtr1, ds_Complex *zPtr2, double theta);
void ds_esigz(ds_Complex *zPtr1, ds_Complex *zPtr2, double theta);
void ds_unitary(ds_Complex *zPtr1, ds_Complex *zPtr2,
             double alpha, double beta, double theta);

void ds_one_qubit_indices(int n, int q, int *iPtr, int *jPtr);
void ds_controlled_indices(int n, int q1, int q2, int *iPtr, int *jPtr);
void ds_swap_indices(int n, int q1, int q2, int *iPtr, int *jPtr);
void ds_two_qubit_indices(int n, int q1, int q2, int *i1Ptr, int *i2Ptr, int *i3Ptr, int*i4Ptr);

void ds_irot(ds_Register reg, int q, double theta, int time);
void ds_xrot(ds_Register reg, int q, double theta, int time);
void ds_yrot(ds_Register reg, int q, double theta, int time);
void ds_zrot(ds_Register reg, int q, double theta, int time);

void ds_X(ds_Register reg, int q, int time);
void ds_Z(ds_Register reg, int q, int time);
void ds_XZ(ds_Register reg, int q, int time);

void ds_X_no_error(ds_Register reg, int q, int time);
void ds_Z_no_error(ds_Register reg, int q, int time);
void ds_XZ_no_error(ds_Register reg, int q, int time);

void ds_cnot(ds_Register reg, int qcont, int qtarg, int time);
void ds_swap(ds_Register reg, int q1, int q2, int time);
void ds_cnots(ds_Register reg, int qcont, int qtarg, int time);
void ds_phase(ds_Register reg, int q, double theta, int time);
void ds_lerr(ds_Register reg, int q, int time);
void ds_global_error(ds_Register reg);
void ds_cphase(ds_Register reg, int qcont, int qtarg, double theta, int time);
void ds_cphases(ds_Register reg, int qcont, int qtarg, double theta, int time);
void ds_Hadamard(ds_Register reg, int q, int time);

/* IBM physical gates */
void ds_u3_gate(ds_Complex *zPtr0, ds_Complex *zPtr1, double theta, double phi, double lambda);
void ds_u2_gate(ds_Complex *zPtr0, ds_Complex *zPtr1, double phi, double lambda);
void ds_u1_gate(ds_Complex *zPtr0, ds_Complex *zPtr1, double lambda);

void ds_U3(ds_Register reg, int q, int time, double theta, double phi, double lambda);
void ds_U2(ds_Register reg, int q, int time, double phi, double lambda);
void ds_U1(ds_Register reg, int q, int time, double lambda);

/* compound operators act left to right : ds_Hcnot == H then ds_cnot */
void ds_Hcnot(ds_Register reg, int qcont, int qtarg, int time);
void ds_cnotH(ds_Register reg, int qcont, int qtarg, int time);
void ds_Hcnots(ds_Register reg, int qcont, int qtarg, int time);
void ds_scnotH(ds_Register reg, int qcont, int qtarg, int time);


/*-------General Shor's Algorithm---*/

void ds_generalqft(ds_Register reg,int N,int tq,int bq,int maxangle);
void ds_generalqftinv(ds_Register reg,int N, int tq, int bq,int maxangle);
void ds_addergen(ds_Register reg,int N,int initial, int inverse,int tq,int bq);
void ds_contaddergen(ds_Register reg,int N, int initial,int inverse,int tq,int bq,int control);
void ds_addmodgen(ds_Register reg, int N, int a, int mod, int inverse,int tq,int bq,int ancilla,int control,int maxangle);
void ds_toffoli(ds_Register reg,int qcont1,int qcont2,int qtarg);
void ds_contmultgen(ds_Register reg,int N, int a, int mod, int inverse, int tq, int bq,int control,int maxangle);
void ds_contswap(ds_Register reg, int tq,int bq,int control);
void ds_cswap(ds_Register reg,int control,int qtarg1,int qtarg2);
void ds_ugategen(ds_Register reg, int N, int a, int mod, int ainv, int tq,int bq,int control,int maxangle);
void ds_shorgen(ds_Register reg,int N, int a, int mod,int tq,int bq,int maxangle);
void ds_shortrickgen(ds_Register reg,int N,int a,int mod,int tq,int bq,int maxangle,int j,FILE *out);


/*-------Linear Shor's Algorithm-----------------*/

void ds_scnot(ds_Register reg, int qcont, int qtarg, int time);
void ds_sqrtnot(ds_Register reg, int qcont, int qtarg, int dagger, int time);
void ds_phases(ds_Register reg,int qcont, int q, double theta, int time);
void ds_linearqft(ds_Register reg,int N,int inverse,int tq,int bq,int maxangle);
void ds_adder(ds_Register reg, int initial, int N, int inverse,int ci,int cf,int carry, int time);
void ds_contadder(ds_Register reg, int initial, int N, int inverse, int ci, int cf,int up);
void ds_addmod1(ds_Register reg, int N,int a,int mod, int tq, int bq, int maxangle);
void ds_conmodinv(ds_Register reg,int N,int a,int mod,int tq,int bq,int maxangle);
void ds_conmod1(ds_Register reg, int N,int a, int mod,int tq,int bq,int maxangle);
void ds_conswaptoffoli(ds_Register reg, int qcont, int qtarg1, int qtarg2);
void ds_linearmesh(ds_Register reg, int tqmesh, int bqmesh,int inverse);
void ds_swapupdown(ds_Register reg, int qstart, int number,int up);
void ds_linearconswap(ds_Register reg, int N, int tq, int bq);
void ds_addmod1inv(ds_Register reg, int N,int a,int mod, int tq, int bq, int maxangle);
void ds_ugate1(ds_Register reg,int N,int a,int mod,int ainv,int tq,int bq,int maxangle);
void ds_shortrick(ds_Register reg,int N,int a,int mod,int tq,int bq,int j, FILE *out,int maxangle);

/*--------------------------------------------------------------------*/

void ds_countstate(ds_Register reg);
void ds_printmod(ds_Register reg, FILE *out);
int ds_main2(int num, int a, char run[5], int linear,int maxangle,int j,double p,double sigma);

void ds_qec5en(ds_Register reg, int q4, int q3, int q2, int q1, int q0);
void ds_qec5dec(ds_Register reg, int q4, int q3, int q2, int q1, int q0);
int ds_qec5cor(ds_Register reg, int q4, int q3, int q2, int q1, int q0);

void ds_qec7en0(int q1, int q2, int q3, int q4, int q5, int q6, int q7);
void ds_qec7dec(int q1, int q2, int q3, int q4, int q5, int q6, int q7);
int ds_qec7cor(int q1, int q2, int q3, int q4, int q5, int q6, int q7, int q8);

void ds_mesh(ds_Register reg, int nq_L);
void ds_swapreg(ds_Register reg, int nq_L);

int ds_measure(ds_Register reg, int nq2m, int *lq2m);

double ds_set_measure(ds_Register reg, int nq2m, int *lq2m, int val);

#endif
