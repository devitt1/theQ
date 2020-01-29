#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "sim.h"
#include "norm.h"

double ds_Pi, ds_Pio2, ds_Pio4, ds_Pio8, ds_root2_2;

ds_Register ds_create_register(int nq_L, double err_L, double sigma_L)
{
   ds_Register reg;
   FILE *out=fopen("mem_regs.txt","w");

   reg.nq = nq_L;
   reg.nc = pow(2, nq_L);
   reg.err = err_L;
   reg.sigma = sigma_L;
   reg.n_errors = (int*) calloc(3,sizeof(int));
   /*
   reg.n_errors[0] = 0;
   reg.n_errors[1] = 0;
   reg.n_errors[2] = 0;*/
   fprintf(out, "%d qubits!\n",nq_L); 

   reg.state = (ds_Complex *)calloc(reg.nc, sizeof(ds_Complex));
   if (reg.state == NULL) {
      fprintf(out, "Insufficient memory for statevector\n");
      exit(EXIT_FAILURE);
   }
   else {
      fprintf(out, "statevecotr allocated \n");
   }

   reg.steps = (int *)calloc(reg.nq, sizeof(int));
   if (reg.steps == NULL) {
      fprintf(out, "Insufficient memory for time steps\n");
      exit(EXIT_FAILURE);
   }
   else {
      fprintf(out, "time steps allocated \n");
   }
   fclose(out);

   return reg;
}

void ds_destroy_register(ds_Register reg)
{
   free(reg.steps);
   free(reg.state);
}

void ds_equate_registers(ds_Register reg1, ds_Register reg2)
{
   int i;

   for (i=0; i<reg2.nc; i++) reg1.state[i] = reg2.state[i];
   for (i=0; i<reg2.nq; i++) reg1.steps[i] = reg2.steps[i];
   reg1.nq = reg2.nq;
   reg1.nc = reg2.nc;
   reg1.err = reg2.err;
   reg1.sigma = reg2.sigma;
}

void ds_initialize_simulator(long seed)
{

   //ds_srand32(time(NULL));
   ds_srand32(seed);

   ds_Pi = acos(-1);
   ds_Pio2 = asin(1);
   ds_Pio4 = atan(1);
   ds_Pio8 = ds_Pio4 / 2;
   ds_root2_2 = 1 / sqrt(2);
}

void ds_clearreg(ds_Register reg)
{
   int i;

   for (i=0; i<reg.nc; i++) {
      reg.state[i].x = 0;
      reg.state[i].y = 0;
   }

   for (i=0; i<reg.nq; i++) {
      reg.steps[i] = 0;
   }
}
   

void ds_set_state(ds_Register reg, int n, double x, double y)
{
   reg.state[n].x = x;
   reg.state[n].y = y;
}

int ds_query_state(ds_Register reg, int n, double tol)
{
   if (reg.state[n].x >= 1-tol) return 1;
   return 0;
}

void dosmth(){
    FILE *out=fopen("out.dat","w");
    fprintf(out, "success from sim.c!");
    fclose(out);
}

void ds_print(ds_Register reg)
{
   int i;
   FILE *out=fopen("reg_out.dat","w");

   for (i=0; i<reg.nc; i++)
      if (reg.state[i].x !=0 || reg.state[i].y !=0 )
	 fprintf(out, "%d: %g: %g \n", i, reg.state[i].x, reg.state[i].y);
   fprintf(out, "\n");
   fclose(out);

}

void ds_change_err(ds_Register *regPtr, double err_L)
{
   regPtr->err = err_L;
}

void ds_change_sigma(ds_Register *regPtr, double sigma_L)
{
   regPtr->sigma = sigma_L;
}

void ds_update(ds_Register reg, int q1, int q2)
{
   while (reg.steps[q1] < reg.steps[q2]) {
      ds_lerr(reg, q1, 1);
      reg.steps[q1]++;
   }
   while (reg.steps[q2] < reg.steps[q1]) {
      ds_lerr(reg, q2, 1);
      reg.steps[q2]++;
   }
}

void ds_global_update(ds_Register reg)
{
   int q, max_steps = 0, i;

   for (q=0; q<reg.nq; q++)
      if (reg.steps[q] > max_steps)
	 max_steps = reg.steps[q];

   for (q=0; q<reg.nq; q++)
      for (i=reg.steps[q]; i<max_steps; i++)
	 ds_lerr(reg, q, 1);

   for (q=0; q<reg.nq; q++)
      reg.steps[q] = max_steps;
}

ds_Complex ds_eitheta(double theta)
{
   ds_Complex z;

   z.x = cos(theta);
   z.y = sin(theta);

   return z;
}

ds_Complex ds_add(ds_Complex z1, ds_Complex z2)
{
   ds_Complex z3;

   z3.x = z1.x + z2.x;
   z3.y = z1.y + z2.y;

   return z3;
}

ds_Complex ds_multiply(ds_Complex z1, ds_Complex z2)
{
   ds_Complex z3;

   z3.x = z1.x*z2.x - z1.y*z2.y;
   z3.y = z1.x*z2.y + z1.y*z2.x;

   return z3;
}

ds_Complex ds_zstarz(ds_Complex z1, ds_Complex z2)
{
   ds_Complex z3;

   z3.x = z1.x*z2.x + z1.y*z2.y;
   z3.y = z1.x*z2.y - z1.y*z2.x;

   return z3;
}

double ds_modsq(ds_Complex z)
{
   return (z.x)*(z.x)+(z.y)*(z.y);
}

double ds_inner_product(ds_Register reg1, ds_Register reg2)
{
   int i;
   ds_Complex temp;

   temp.x = temp.y = 0;
   for (i=0; i<reg1.nc; i++) temp = ds_add(temp, ds_zstarz(reg1.state[i], reg2.state[i]));

   return ds_modsq(temp);
}

void ds_esigi(ds_Complex *zPtr1, ds_Complex *zPtr2, double theta)
{
   double c = cos(theta/2), s = sin(theta/2);
   ds_Complex z1 = *zPtr1, z2 = *zPtr2;
   zPtr1->x = c*(z1.x) - s*(z2.y);
   zPtr1->y = s*(z1.x) + c*(z2.y);
   zPtr2->x = c*(z1.x) - s*(z2.y);
   zPtr2->y = s*(z1.x) + c*(z2.y);
}

void ds_esigx(ds_Complex *zPtr1, ds_Complex *zPtr2, double theta)
{
   double c = cos(theta/2), s = sin(theta/2);
   ds_Complex z1 = *zPtr1, z2 = *zPtr2;
   zPtr1->x =  c*(z1.x) -s*(z2.y);
   zPtr1->y =  c*(z1.y) +s*(z2.x);
   zPtr2->x = -s*(z1.y) +c*(z2.x);
   zPtr2->y =  s*(z1.x) +c*(z2.y);
}

void ds_esigy(ds_Complex *zPtr1, ds_Complex *zPtr2, double theta)
{
   double c = cos(theta/2), s = sin(theta/2);
   ds_Complex z1 = *zPtr1, z2 = *zPtr2;
   zPtr1->x =  c*(z1.x) +s*(z2.x);
   zPtr1->y =  c*(z1.y) +s*(z2.y);
   zPtr2->x = -s*(z1.x) +c*(z2.x);
   zPtr2->y = -s*(z1.y) +c*(z2.y);
}

void ds_u3_gate(ds_Complex *zPtr0, ds_Complex *zPtr1,
         double theta, double phi, double lambda)
{
   double c   = cos(theta/2),       s = sin(theta/2);
   double cp  = cos(phi),          sp = sin(phi);
   double cl  = cos(lambda),       sl = sin(lambda);
   double clp = cos(lambda+phi),  slp = sin(lambda+phi);

   ds_Complex z1 = *zPtr1, z0 = *zPtr0;

   zPtr1->x = s*(cp*z0.x - sp*z0.y) + c*(z1.x*clp - z1.y*slp);
   zPtr1->y = s*(cp*z0.y + sp*z0.x) + c*(z1.y*clp + z1.x*slp);
   zPtr0->x = c*z0.x - s*cl*z1.x + s*sl*z1.y;
   zPtr0->y = c*z0.y - s*cl*z1.y - s*sl*z1.x;
}

void ds_u2_gate(ds_Complex *zPtr0, ds_Complex *zPtr1,
         double phi, double lambda)
{
   ds_u3_gate(zPtr0, zPtr1, ds_Pi/2.0, phi, lambda);
}

void ds_u1_gate(ds_Complex *zPtr0, ds_Complex *zPtr1,
         double lambda)
{
   ds_u3_gate(zPtr0, zPtr1, 0, 0, lambda);
}

void ds_U3(ds_Register reg, int q, int time, double theta, double phi, double lambda)
{
   int i, j, k;

   for (i=0; i<reg.nc/2; i++) {
      ds_one_qubit_indices(i, q, &j, &k);
      ds_u3_gate(reg.state+j, reg.state+k, theta, phi, lambda);
   }

   ds_lerr(reg, q, time);
   reg.steps[q]+=time;
}


	//ds_U2(ds_reg, 0, 1, 0, ds_Pi);
void ds_U2(ds_Register reg, int q, int time, double phi, double lambda)
{
   int i, j, k;

   for (i=0; i<reg.nc/2; i++) {
      ds_one_qubit_indices(i, q, &j, &k);
      ds_u2_gate(reg.state+j, reg.state+k, phi, lambda);
   }

   ds_lerr(reg, q, time);
   reg.steps[q]+=time;
}

void ds_U1(ds_Register reg, int q, int time, double lambda)
{
   int i, j, k;

   for (i=0; i<reg.nc/2; i++) {
      ds_one_qubit_indices(i, q, &j, &k);
      ds_u1_gate(reg.state+j, reg.state+k, lambda);
   }

   ds_lerr(reg, q, time);
   reg.steps[q]+=time;
}

void ds_esigz(ds_Complex *zPtr1, ds_Complex *zPtr2, double theta)
{
   double c = cos(theta/2), s = sin(theta/2);
   ds_Complex z1 = *zPtr1, z2 = *zPtr2;
   zPtr1->x =  c*(z1.x) -s*(z1.y);
   zPtr1->y =  c*(z1.y) +s*(z1.x);
   zPtr2->x =  c*(z2.x) +s*(z2.y);
   zPtr2->y =  c*(z2.y) -s*(z2.x);
}

void ds_unitary(ds_Complex *zPtr1, ds_Complex *zPtr2,
             double alpha, double beta, double theta)
{
   double c = cos(theta/2), s = sin(theta/2);
   double c11 = cos(+alpha/2+beta/2);
   double s11 = sin(+alpha/2+beta/2);
   double c12 = cos(+alpha/2-beta/2);
   double s12 = sin(+alpha/2-beta/2);
   ds_Complex z1 = *zPtr1, z2 = *zPtr2;

   zPtr1->x =  c*(c11*z1.x-s11*z1.y) + s*(c12*z2.x-s12*z2.y);
   zPtr1->y =  c*(c11*z1.y+s11*z1.x) + s*(c12*z2.y+s12*z2.x);
   zPtr2->x = -s*(c12*z1.x+s12*z1.y) + c*(c11*z2.x+s11*z2.y);
   zPtr2->y = -s*(c12*z1.y-s12*z1.x) + c*(c11*z2.y-s11*z2.x);
}

/* Function used when performing a single qubit operation. n runs from 0
   to nc/2 - 1 and simply keeps track of how far through the array the
   calculation is. q runs from 0 to nq - 1 and denotes which qubit the
   operation is being performed on. The function outputs the next pair of
   indices to apply the desired transformation to.
   
     returns: nth pair of indices *iPtr and *jPtr, where in state *iPtr
              qubit q is not set and in state *jPtr qubit q is set
     e.g. for 3 qubits: n=0, q=1 -> *iPtr=0 (000) *jPtr=1 (001)
     e.g. for 3 qubits: n=0, q=2 -> *iPtr=0 (000) *jPtr=2 (010)
     e.g. for 3 qubits: n=1, q=3 -> *iPtr=1 (001) *jPtr=1 (101)
*/
void ds_one_qubit_indices(int n, int q, int *iPtr, int *jPtr)
{
   /* l is used as a mask to calculate nl --- the low part of n */
   int l, nl;

   /* calculate mask of 1s for section of n less significant than bit q */
   l = 1;
   l <<= q;
   l--;

   /* calculate section of n less significant than bit q */
   nl = n & l;

   /* remove section of n less significant than bit q and shift n higher for
      later insertion of bit q */
   n -= nl;
   n <<= 1;

   /* calculate value to be used to insert value of bit q */
   l++;

   *iPtr = n + 0 + nl;
   *jPtr = n + l + nl;
}

/* See comments for ds_one_qubit_indices(). */
void ds_controlled_indices(int n, int qcont, int qtarg, int *iPtr, int *jPtr)
{
   int m, l, nm, nl, qm, ql;

   if (qcont > qtarg) {
      qm = qcont;
      ql = qtarg;
   }
   else {
      qm = qtarg;
      ql = qcont;
   }

   /* calculate mask of 1s for section of n less significant than bit ql */
   l = 1;
   l <<= ql;
   l--;

   /* calculate section of n less significant than bit ql */
   nl = n & l;

   /* remove section of n less significant than bit ql and shift n higher for
      later insertion of bit q2 */
   n -= nl;
   n <<= 1;

   /* calculate mask of 1s for section of n less significant than bit qm but
      more significant than bit ql */
   m = 1;
   m <<= qm;
   m--;

   /* calculate section of n less significant than bit qm but more significant
      than bit ql */
   nm = n & m;

   /* remove section of n less significant than bit qm and shift n higher for
      later insertion of bit qm. Add nm and nl to complete preparation of n */
   n -= nm;
   n <<= 1;
   n += nm + nl;

   /* calculate values to be used to insert values of bits q1 and q2 */
   m++;
   l++;

   if (qcont > qtarg) {
      n += m;
      *iPtr = n;
      *jPtr = n + l;
   }
   else {
      n += l;
      *iPtr = n;
      *jPtr = n + m;
   }
}

/* See comments for ds_one_qubit_indices(). */
void ds_swap_indices(int n, int q1, int q2, int *iPtr, int *jPtr)
{
   int m, l, nm, nl, qm, ql;

   if (q1 > q2) {
      qm = q1;
      ql = q2;
   }
   else {
      qm = q2;
      ql = q1;
   }

   /* calculate mask of 1s for section of n less significant than bit ql */
   l = 1;
   l <<= ql;
   l--;

   /* calculate section of n less significant than bit ql */
   nl = n & l;

   /* remove section of n less significant than bit ql and shift n higher for
      later insertion of bit q2 */
   n -= nl;
   n <<= 1;

   /* calculate mask of 1s for section of n less significant than bit qm but
      more significant than bit ql */
   m = 1;
   m <<= qm;
   m--;

   /* calculate section of n less significant than bit qm but more significant
      than bit ql */
   nm = n & m;

   /* remove section of n less significant than bit qm and shift n higher for
      later insertion of bit qm. Add nm and nl to complete preparation of n */
   n -= nm;
   n <<= 1;
   n += nm + nl;

   /* calculate values to be used to insert values of bits q1 and q2 */
   m++;
   l++;

   *iPtr = n + l;
   *jPtr = n + m;
}

void ds_two_qubit_indices(int n, int q1, int q2, int *i1Ptr, int *i2Ptr, int *i3Ptr, int*i4Ptr)
{
   int m, l, nm, nl, qm, ql;

   if (q1 > q2) {
      qm = q1;
      ql = q2;
   }
   else {
      qm = q2;
      ql = q1;
   }

   /* calculate mask of 1s for section of n less significant than bit ql */
   l = 1;
   l <<= ql;
   l--;

   /* calculate section of n less significant than bit ql */
   nl = n & l;

   /* remove section of n less significant than bit ql and shift n higher for
      later insertion of bit q2 */
   n -= nl;
   n <<= 1;

   /* calculate mask of 1s for section of n less significant than bit qm but
      more significant than bit ql */
   m = 1;
   m <<= qm;
   m--;

   /* calculate section of n less significant than bit qm but more significant
      than bit ql */
   nm = n & m;

   /* remove section of n less significant than bit qm and shift n higher for
      later insertion of bit qm. Add nm and nl to complete preparation of n */
   n -= nm;
   n <<= 1;
   n += nm + nl;

   /* calculate values to be used to insert values of bits q1 and q2 */
   m++;
   l++;

   *i1Ptr = n;
   *i2Ptr = n + l;
   *i3Ptr = n + m;
   *i4Ptr = n + m + l;
}

void ds_irot(ds_Register reg, int q, double theta, int time)
{
   int i, j, k;

   for (i=0; i<reg.nc/2; i++) {
      ds_one_qubit_indices(i, q, &j, &k);
      ds_esigi(reg.state+j, reg.state+k, theta);
   }

   ds_lerr(reg, q, time);
   reg.steps[q]+=time;
}

void ds_xrot(ds_Register reg, int q, double theta, int time)
{
   int i, j, k;

   for (i=0; i<reg.nc/2; i++) {
      ds_one_qubit_indices(i, q, &j, &k);
      ds_esigx(reg.state+j, reg.state+k, theta);
   }

   ds_lerr(reg, q, time);
   reg.steps[q]+=time;
}

void ds_yrot(ds_Register reg, int q, double theta, int time)
{
   int i, j, k;

   for (i=0; i<reg.nc/2; i++) {
      ds_one_qubit_indices(i, q, &j, &k);
      ds_esigy(reg.state+j, reg.state+k, theta);
   }

   ds_lerr(reg, q, time);
   reg.steps[q]+=time;
}

void ds_zrot(ds_Register reg, int q, double theta, int time)
{
   int i, j, k;

   for (i=0; i<reg.nc/2; i++) {
      ds_one_qubit_indices(i, q, &j, &k);
      ds_esigz(reg.state+j, reg.state+k, theta);
   }

   ds_lerr(reg, q, time);
   reg.steps[q]+=time;
}

void ds_X_no_error(ds_Register reg, int q, int time)
{
   int i, j, k;
   ds_Complex z;

   for (i=0; i<reg.nc/2; i++) {
      ds_one_qubit_indices(i, q, &j, &k);
      z = reg.state[j];
      reg.state[j] = reg.state[k];
      reg.state[k] = z;
   }

   reg.steps[q]+=time;
}
void ds_X(ds_Register reg, int q, int time)
{
   int i, j, k;
   ds_Complex z;

   for (i=0; i<reg.nc/2; i++) {
      ds_one_qubit_indices(i, q, &j, &k);
      z = reg.state[j];
      reg.state[j] = reg.state[k];
      reg.state[k] = z;
   }

   ds_lerr(reg, q, time);
   reg.steps[q]+=time;
}

void ds_Z(ds_Register reg, int q, int time)
{
   int i, j, k;
   ds_Complex z;

   for (i=0; i<reg.nc/2; i++) {
      ds_one_qubit_indices(i, q, &j, &k);
      reg.state[k].x = -reg.state[k].x;
      reg.state[k].y = -reg.state[k].y;
   }

   ds_lerr(reg, q, time);

   reg.steps[q]+=time;
}

void ds_Z_no_error(ds_Register reg, int q, int time)
{
   int i, j, k;
   ds_Complex z;

   for (i=0; i<reg.nc/2; i++) {
      ds_one_qubit_indices(i, q, &j, &k);
      reg.state[k].x = -reg.state[k].x;
      reg.state[k].y = -reg.state[k].y;
   }

   reg.steps[q]+=time;
}

void ds_XZ_no_error(ds_Register reg, int q, int time){
   
   int i, j, k;
   ds_Complex z;

   for (i=0; i<reg.nc/2; i++) {
      ds_one_qubit_indices(i, q, &j, &k);
      z = reg.state[j];
      reg.state[j].x = -reg.state[k].x;
      reg.state[j].y = -reg.state[k].y;
      reg.state[k] = z;
   }

   reg.steps[q]+=time;
}

void ds_XZ(ds_Register reg, int q, int time)
{
   int i, j, k;
   ds_Complex z;

   for (i=0; i<reg.nc/2; i++) {
      ds_one_qubit_indices(i, q, &j, &k);
      z = reg.state[j];
      reg.state[j].x = -reg.state[k].x;
      reg.state[j].y = -reg.state[k].y;
      reg.state[k] = z;
   }

   ds_lerr(reg, q, time);
   reg.steps[q]+=time;
}

/* Change to discrete errors OR continuous */
void ds_lerr(ds_Register reg, int q, int time)
{
   int i, j, k;
   double p;
   double delta, alpha, beta, theta;
   ds_Complex z;

   if (time == 0) return;
    if (reg.err == 0) return;

   if (reg.err > 0) {
      p = ds_uniform();

      if (p < reg.err/3) {
	 ds_X_no_error(reg, q, 1);

	 //printf("error\n");
         reg.n_errors[0] += 1;
      }
      else if (p < 2*reg.err/3) {
	 ds_Z_no_error(reg, q, 1);

	 //printf("error\n");
         reg.n_errors[1] += 1;
      }
      else if (p < reg.err) {
	 ds_XZ_no_error(reg, q, 1);

	 //printf("error\n");
         reg.n_errors[2] += 1;
      }
   }

   if (reg.sigma > 0) {
      alpha = reg.sigma*ds_norm();
      beta = reg.sigma*ds_norm();
      theta = reg.sigma*ds_norm();
         for (i=0; i<reg.nc/2; i++) {
	 ds_one_qubit_indices(i, q, &j, &k);
	 ds_unitary(reg.state+j, reg.state+k, alpha, beta, theta);
	 }
   }
}

void ds_global_error(ds_Register reg)
{
   int q;

   for (q=0; q<reg.nq; q++) ds_lerr(reg, q, 1);
}

/* Kane compatible if |qcont-qtarg| = 1 */
void ds_cnot(ds_Register reg, int qcont, int qtarg, int time)
{
   int i, j, k;
   ds_Complex z;

   ds_update(reg, qcont, qtarg);

   for (i=0; i<reg.nc/4; i++) {
      ds_controlled_indices(i, qcont, qtarg, &j, &k);
      /* Note that the ds_toffoli gate does not work with pseudo ds_cnot gates
         ie : you cant use ds_esigy(reg.state+j, reg.state+k, Pi) */
      z = reg.state[j];
      reg.state[j] = reg.state[k];
      reg.state[k] = z;
   }

   ds_lerr(reg, qcont, time);
   reg.steps[qcont]+=time;
   ds_lerr(reg, qtarg, time);
   reg.steps[qtarg]+=time;
}

/* Kane compatible if |q1-q2| = 1 */
void ds_swap(ds_Register reg, int q1, int q2, int time)
{
   int i, j, k;
   ds_Complex z;

   ds_update(reg, q1, q2);

   for (i=0; i<reg.nc/4; i++) {
      ds_swap_indices(i, q1, q2, &j, &k);
      z = reg.state[j];
      reg.state[j] = reg.state[k];
      reg.state[k] = z;
   }

   ds_lerr(reg, q1, time);
   reg.steps[q1]+=time;
   ds_lerr(reg, q2, time);
   reg.steps[q2]+=time;
}

void ds_cnots(ds_Register reg, int qcont, int qtarg, int time)
{
   int i, i1, i2, i3, i4;
   ds_Complex z;

   ds_update(reg, qcont, qtarg);

   if (qcont > qtarg) {
      for (i=0; i<reg.nc/4; i++) {
	 ds_two_qubit_indices(i, qcont, qtarg, &i1, &i2, &i3, &i4);
	 z = reg.state[i2];
	 reg.state[i2] = reg.state[i4];
	 reg.state[i4] = reg.state[i3];
	 reg.state[i3] = z;
      }
   }
   else {
      for (i=0; i<reg.nc/4; i++) {
	 ds_two_qubit_indices(i, qcont, qtarg, &i1, &i2, &i3, &i4);
	 z = reg.state[i2];
	 reg.state[i2] = reg.state[i3];
	 reg.state[i3] = reg.state[i4];
	 reg.state[i4] = z;
      }
   }

   ds_lerr(reg, qcont, time);
   reg.steps[qcont]+=time;
   ds_lerr(reg, qtarg, time);
   reg.steps[qtarg]+=time;
}

/* Kane compatible if |q-x| = 1 */
void ds_phase(ds_Register reg, int q, double theta, int time)
{
   int i, j, k;

   for (i=0; i<reg.nc/2; i++) {
      ds_one_qubit_indices(i, q, &j, &k);
      reg.state[k] = ds_multiply(ds_eitheta(theta), reg.state[k]);
   }

   ds_lerr(reg, q, time);
   reg.steps[q]+=time;
}

void ds_cphase(ds_Register reg, int qcont, int qtarg, double theta, int time)
{
   int i, j, k;

   ds_update(reg, qcont, qtarg);

   for (i=0; i<reg.nc/4; i++) {
      ds_controlled_indices(i, qcont, qtarg, &j, &k);
      reg.state[k] = ds_multiply(ds_eitheta(theta), reg.state[k]);
   }

   ds_lerr(reg, qcont, time);
   reg.steps[qcont]+=time;
   ds_lerr(reg, qtarg, time);
   reg.steps[qtarg]+=time;
}

void ds_cphases(ds_Register reg, int qcont, int qtarg, double theta, int time)
{
   int i, j, k;
   ds_Complex z;

   ds_update(reg, qcont, qtarg);

   for (i=0; i<reg.nc/4; i++) {
      ds_controlled_indices(i, qcont, qtarg, &j, &k);
      reg.state[k] = ds_multiply(ds_eitheta(theta), reg.state[k]);
   }

   for (i=0; i<reg.nc/4; i++) {
      ds_swap_indices(i, qcont, qtarg, &j, &k);
      z = reg.state[j];
      reg.state[j] = reg.state[k];
      reg.state[k] = z;
   }

   ds_lerr(reg, qcont, time);
   reg.steps[qcont]+=time;
   ds_lerr(reg, qtarg, time);
   reg.steps[qtarg]+=time;
}

void ds_Hadamard(ds_Register reg, int q, int time)
{
   int i, j, k;
   ds_Complex z1, z2;

   for (i=0; i<reg.nc/2; i++) {
      ds_one_qubit_indices(i, q, &j, &k);
      z1 = reg.state[j];
      z2 = reg.state[k];
      reg.state[j].x = ds_root2_2*(z1.x+z2.x);
      reg.state[j].y = ds_root2_2*(z1.y+z2.y);
      reg.state[k].x = ds_root2_2*(z1.x-z2.x);
      reg.state[k].y = ds_root2_2*(z1.y-z2.y);
   }

   ds_lerr(reg, q, time);
   reg.steps[q]+=time;
}

void ds_Hcnot(ds_Register reg, int qcont, int qtarg, int time)
{
   int i, i1, i2, i3, i4;
   ds_Complex z1, z2, z3, z4;

   ds_update(reg, qcont, qtarg);

   if (qcont > qtarg) {
      for (i=0; i<reg.nc/4; i++) {
	 ds_two_qubit_indices(i, qcont, qtarg, &i1, &i2, &i3, &i4);
	 z1 = reg.state[i1];
	 z2 = reg.state[i2];
	 z3 = reg.state[i3];
	 z4 = reg.state[i4];
	 reg.state[i1].x = (z1.x+z3.x)*ds_root2_2;
	 reg.state[i1].y = (z1.y+z3.y)*ds_root2_2;
	 reg.state[i2].x = (z2.x+z4.x)*ds_root2_2;
	 reg.state[i2].y = (z2.y+z4.y)*ds_root2_2;
	 reg.state[i3].x = (z2.x-z4.x)*ds_root2_2;
	 reg.state[i3].y = (z2.y-z4.y)*ds_root2_2;
	 reg.state[i4].x = (z1.x-z3.x)*ds_root2_2;
	 reg.state[i4].y = (z1.y-z3.y)*ds_root2_2;
      }
   }
   else {
      for (i=0; i<reg.nc/4; i++) {
	 ds_two_qubit_indices(i, qcont, qtarg, &i1, &i2, &i3, &i4);
	 z1 = reg.state[i1];
	 z2 = reg.state[i2];
	 z3 = reg.state[i3];
	 z4 = reg.state[i4];
	 reg.state[i1].x = (z1.x+z2.x)*ds_root2_2;
	 reg.state[i1].y = (z1.y+z2.y)*ds_root2_2;
	 reg.state[i2].x = (z3.x-z4.x)*ds_root2_2;
	 reg.state[i2].y = (z3.y-z4.y)*ds_root2_2;
	 reg.state[i3].x = (z3.x+z4.x)*ds_root2_2;
	 reg.state[i3].y = (z3.y+z4.y)*ds_root2_2;
	 reg.state[i4].x = (z1.x-z2.x)*ds_root2_2;
	 reg.state[i4].y = (z1.y-z2.y)*ds_root2_2;
      }
   }

   ds_lerr(reg, qcont, time);
   reg.steps[qcont]+=time;
   ds_lerr(reg, qtarg, time);
   reg.steps[qtarg]+=time;
}

void ds_cnotH(ds_Register reg, int qcont, int qtarg, int time)
{
   int i, i1, i2, i3, i4;
   ds_Complex z1, z2, z3, z4;

   ds_update(reg, qcont, qtarg);

   if (qcont > qtarg) {
      for (i=0; i<reg.nc/4; i++) {
	 ds_two_qubit_indices(i, qcont, qtarg, &i1, &i2, &i3, &i4);
	 z1 = reg.state[i1];
	 z2 = reg.state[i2];
	 z3 = reg.state[i3];
	 z4 = reg.state[i4];
	 reg.state[i1].x = (z1.x+z4.x)*ds_root2_2;
	 reg.state[i1].y = (z1.y+z4.y)*ds_root2_2;
	 reg.state[i2].x = (z2.x+z3.x)*ds_root2_2;
	 reg.state[i2].y = (z2.y+z3.y)*ds_root2_2;
	 reg.state[i3].x = (z1.x-z4.x)*ds_root2_2;
	 reg.state[i3].y = (z1.y-z4.y)*ds_root2_2;
	 reg.state[i4].x = (z2.x-z3.x)*ds_root2_2;
	 reg.state[i4].y = (z2.y-z3.y)*ds_root2_2;
      }
   }
   else {
      for (i=0; i<reg.nc/4; i++) {
	 ds_two_qubit_indices(i, qcont, qtarg, &i1, &i2, &i3, &i4);
	 z1 = reg.state[i1];
	 z2 = reg.state[i2];
	 z3 = reg.state[i3];
	 z4 = reg.state[i4];
	 reg.state[i1].x = (z1.x+z4.x)*ds_root2_2;
	 reg.state[i1].y = (z1.y+z4.y)*ds_root2_2;
	 reg.state[i2].x = (z1.x-z4.x)*ds_root2_2;
	 reg.state[i2].y = (z1.y-z4.y)*ds_root2_2;
	 reg.state[i3].x = (z2.x+z3.x)*ds_root2_2;
	 reg.state[i3].y = (z2.y+z3.y)*ds_root2_2;
	 reg.state[i4].x = (-z2.x+z3.x)*ds_root2_2;
	 reg.state[i4].y = (-z2.y+z3.y)*ds_root2_2;
      }
   }

   ds_lerr(reg, qcont, time);
   reg.steps[qcont]+=time;
   ds_lerr(reg, qtarg, time);
   reg.steps[qtarg]+=time;
}

void ds_phases(ds_Register reg,int qcont, int q, double theta, int time)
{
   int i, j, k;
   ds_Complex z;

   for (i=0; i<reg.nc/2; i++) {
      ds_one_qubit_indices(i, q, &j, &k);
      reg.state[k] = ds_multiply(ds_eitheta(theta), reg.state[k]);
   }

   for (i=0; i<reg.nc/4; i++) {
      ds_swap_indices(i, qcont, q, &j, &k);
      z = reg.state[j];
      reg.state[j] = reg.state[k];
      reg.state[k] = z;
   }

   ds_lerr(reg, q, time);
   reg.steps[q]+=time;
   ds_lerr(reg, qcont, time);
   reg.steps[qcont]+=time;
}



void ds_Hcnots(ds_Register reg, int qcont, int qtarg, int time)
{
   int i, i1, i2, i3, i4;
   ds_Complex z1, z2, z3, z4;

   ds_update(reg, qcont, qtarg);

   if (qcont > qtarg) {
      for (i=0; i<reg.nc/4; i++) {
	 ds_two_qubit_indices(i, qcont, qtarg, &i1, &i2, &i3, &i4);
	 z1 = reg.state[i1];
	 z2 = reg.state[i2];
	 z3 = reg.state[i3];
	 z4 = reg.state[i4];
	 reg.state[i1].x = (z1.x+z3.x)*ds_root2_2;
	 reg.state[i1].y = (z1.y+z3.y)*ds_root2_2;
	 reg.state[i2].x = (z2.x-z4.x)*ds_root2_2;
	 reg.state[i2].y = (z2.y-z4.y)*ds_root2_2;
	 reg.state[i3].x = (z2.x+z4.x)*ds_root2_2;
	 reg.state[i3].y = (z2.y+z4.y)*ds_root2_2;
	 reg.state[i4].x = (z1.x-z3.x)*ds_root2_2;
	 reg.state[i4].y = (z1.y-z3.y)*ds_root2_2;
      }
   }
   else {
      for (i=0; i<reg.nc/4; i++) {
	 ds_two_qubit_indices(i, qcont, qtarg, &i1, &i2, &i3, &i4);
	 z1 = reg.state[i1];
	 z2 = reg.state[i2];
	 z3 = reg.state[i3];
	 z4 = reg.state[i4];
	 reg.state[i1].x = (z1.x+z2.x)*ds_root2_2;
	 reg.state[i1].y = (z1.y+z2.y)*ds_root2_2;
	 reg.state[i2].x = (z3.x+z4.x)*ds_root2_2;
	 reg.state[i2].y = (z3.y+z4.y)*ds_root2_2;
	 reg.state[i3].x = (z3.x-z4.x)*ds_root2_2;
	 reg.state[i3].y = (z3.y-z4.y)*ds_root2_2;
	 reg.state[i4].x = (z1.x-z2.x)*ds_root2_2;
	 reg.state[i4].y = (z1.y-z2.y)*ds_root2_2;
      }
   }

   ds_lerr(reg, qcont, time);
   reg.steps[qcont]+=time;
   ds_lerr(reg, qtarg, time);
   reg.steps[qtarg]+=time;
}

void ds_scnotH(ds_Register reg, int qcont, int qtarg, int time)
{
   int i, i1, i2, i3, i4;
   ds_Complex z1, z2, z3, z4;

   ds_update(reg, qcont, qtarg);

   if (qcont > qtarg) {
      for (i=0; i<reg.nc/4; i++) {
	 ds_two_qubit_indices(i, qcont, qtarg, &i1, &i2, &i3, &i4);
	 z1 = reg.state[i1];
	 z2 = reg.state[i2];
	 z3 = reg.state[i3];
	 z4 = reg.state[i4];
	 reg.state[i1].x = (z1.x+z4.x)*ds_root2_2;
	 reg.state[i1].y = (z1.y+z4.y)*ds_root2_2;
	 reg.state[i2].x = (z2.x+z3.x)*ds_root2_2;
	 reg.state[i2].y = (z2.y+z3.y)*ds_root2_2;
	 reg.state[i3].x = (z1.x-z4.x)*ds_root2_2;
	 reg.state[i3].y = (z1.y-z4.y)*ds_root2_2;
	 reg.state[i4].x = (-z2.x+z3.x)*ds_root2_2;
	 reg.state[i4].y = (-z2.y+z3.y)*ds_root2_2;
      }
   }
   else {
      for (i=0; i<reg.nc/4; i++) {
	 ds_two_qubit_indices(i, qcont, qtarg, &i1, &i2, &i3, &i4);
	 z1 = reg.state[i1];
	 z2 = reg.state[i2];
	 z3 = reg.state[i3];
	 z4 = reg.state[i4];
	 reg.state[i1].x = (z1.x+z4.x)*ds_root2_2;
	 reg.state[i1].y = (z1.y+z4.y)*ds_root2_2;
	 reg.state[i2].x = (z1.x-z4.x)*ds_root2_2;
	 reg.state[i2].y = (z1.y-z4.y)*ds_root2_2;
	 reg.state[i3].x = (z2.x+z3.x)*ds_root2_2;
	 reg.state[i3].y = (z2.y+z3.y)*ds_root2_2;
	 reg.state[i4].x = (z2.x-z3.x)*ds_root2_2;
	 reg.state[i4].y = (z2.y-z3.y)*ds_root2_2;
      }
   }

   ds_lerr(reg, qcont, time);
   reg.steps[qcont]+=time;
   ds_lerr(reg, qtarg, time);
   reg.steps[qtarg]+=time;
}

/* Top to bottom input : 43210
 * Top to bottom output : 21403 */
void ds_qec5en(ds_Register reg, int q4, int q3, int q2, int q1, int q0)
{
   ds_Hcnots(reg, 3, 2, 1);
   ds_cnots(reg, 4, 3, 1);
   ds_cnots(reg, 2, 1, 1);
   ds_swap(reg, 3, 2, 1);
   ds_Hcnot(reg, 1, 0, 1);
   ds_cnot(reg, 2, 1, 1);
   ds_Hcnot(reg, 2, 3, 1);
   ds_swap(reg, 1, 0, 1);
   ds_cnot(reg, 2, 1, 1);
}

/* Top to bottom input : 21403
 * Top to bottom output : 43210 */
void ds_qec5dec(ds_Register reg, int q4, int q3, int q2, int q1, int q0)
{
   ds_cnot(reg, 2, 1, 1);
   ds_cnotH(reg, 2, 3, 1);
   ds_swap(reg, 1, 0, 1);
   ds_cnot(reg, 2, 1, 1);
   ds_swap(reg, 3, 2, 1);
   ds_cnotH(reg, 1, 0, 1);
   ds_cnots(reg, 3, 4, 1);
   ds_cnots(reg, 1, 2, 1);
   ds_scnotH(reg, 3, 2, 1);
}

/* Note that ds_qec5dec must be called first */
int ds_qec5cor(ds_Register reg, int q4, int q3, int q2, int q1, int q0)
{
   int s, bit[4];

   bit[0] = q3;
   bit[1] = q2;
   bit[2] = q1;
   bit[3] = q0;

   s = ds_measure(reg, 4, bit);

   switch (s) {
      case 0:
	 break;
      case 1:
	 ds_X(reg, 0, 1);
	 break;
      case 2:
	 ds_X(reg, 1, 1);
	 break;
      case 3:
	 ds_Z(reg, 4, 1);
	 ds_X(reg, 1, 1);
	 ds_X(reg, 0, 1);
	 break;
      case 4:
	 ds_X(reg, 2, 1);
	 break;
      case 5:
	 ds_X(reg, 4, 1);
	 ds_X(reg, 2, 1);
	 ds_X(reg, 0, 1);
         break;
      case 6:
	 ds_Z(reg, 4, 1);
	 ds_X(reg, 2, 1);
	 ds_X(reg, 1, 1);
	 break;
      case 7:
	 ds_X(reg, 4, 1);
	 ds_X(reg, 2, 1);
	 ds_X(reg, 1, 1);
	 ds_X(reg, 0, 1);
         break;
      case 8:
	 ds_Z(reg, 4, 1);
	 ds_X(reg, 3, 1);
	 break;
      case 9:
	 ds_X(reg, 3, 1);
	 ds_X(reg, 0, 1);
	 break;
      case 10:
	 ds_X(reg, 4, 1);
	 ds_X(reg, 3, 1);
	 ds_X(reg, 1, 1);
	 break;
      case 11:
	 ds_X(reg, 4, 1);
	 ds_X(reg, 3, 1);
	 ds_X(reg, 1, 1);
	 ds_X(reg, 0, 1);
	 break;
      case 12:
	 ds_Z(reg, 4, 1);
	 ds_X(reg, 3, 1);
	 ds_X(reg, 2, 1);
	 break;
      case 13:
	 ds_X(reg, 4, 1);
	 ds_X(reg, 3, 1);
	 ds_X(reg, 2, 1);
	 ds_X(reg, 0, 1);
	 break;
      case 14:
	 ds_XZ(reg, 4, 1);
	 ds_X(reg, 3, 1);
	 ds_X(reg, 2, 1);
	 ds_X(reg, 1, 1);
	 break;
      case 15:
	 ds_Z(reg, 4, 1);
	 ds_X(reg, 3, 1);
	 ds_X(reg, 2, 1);
	 ds_X(reg, 1, 1);
	 ds_X(reg, 0, 1);
	 break;
   }

   ds_global_error(reg);

   return s;
}

void ds_qec7en0(int a1, int a2, int a3, int a4, int a5, int q, int a6)
{
}

int ds_qec7cor(int q1, int q2, int q3, int q4, int q5, int q6, int q7, int q8)
{
   return 0;
}

/* Kane compatable. nq must be even. */
void ds_mesh(ds_Register reg, int nq_L)
{
   int i, j, s;

   s = 1;
   for (i=nq_L/2; i<nq_L-1; i++) {
      for (j=0; j<s; j++)
         ds_swap(reg, i-2*j, i-2*j-1,1);
      s++;
   }
}

/* Kane compatable. nq_L must be even. */
void ds_swapreg(ds_Register reg, int nq_L)
{
   int i, j;

   for (i=0; i<nq_L/2; i++)
      for (j=0; j<=i; j++)
         ds_swap(reg, nq_L/2+i-2*j, nq_L/2+i-2*j-1,1);
   for (i=nq_L/2-2; i>=0; i--)
      for (j=0; j<=i; j++)
         ds_swap(reg, nq_L/2+i-2*j, nq_L/2+i-2*j-1,1);
}

/* nq2m : number of qubits to ds_measure
   lq2m : list of qubits to ds_measure (sorted from largest to smallest) */
int ds_measure(ds_Register reg, int nq2m, int *lq2m)
{
   int i, j, k, imeas, spread, mask, *masks;
   double r, p, *probabilities, norm;

   ds_global_update(reg);

   /* random value between 0 and 1 */
   r = ds_uniform();

   /* create array of single bit masks for reducing array index */
   masks = (int *)calloc(nq2m, sizeof(int));
   for (i=0; i<nq2m; i++) {
      masks[i] = 1;
      masks[i] <<= lq2m[i];
   }

   /* create array of possible measured values */
   probabilities = (double *)calloc(1 << nq2m, sizeof(double));

   for (i=0; i<reg.nc; i++) {
      /* reduce array index i */
      j = 0;
      for (k=0; k<nq2m; k++) {
         j <<= 1;
         if (i & masks[k]) {
            j++;
         }
      }
      probabilities[j] += ds_modsq(reg.state[i]);
   }

   i = 0;
   p = probabilities[0];
   while (p < r) {
      i++;
      p += probabilities[i];
   }
   
   imeas = i;
   norm = 1/sqrt(probabilities[imeas]);

   mask = 0;
   for (i=0; i<nq2m; i++) mask += masks[i];

   spread = 0;
   i = 1;
   for (j =0; j<nq2m; j++) {
      if (imeas & i) spread += masks[nq2m-1-j];
      i <<= 1;
   }

   for (i=0; i<reg.nc; i++) {
      if ((i & mask) == spread) {
         reg.state[i].x *= norm;
         reg.state[i].y *= norm;
      }
      else {
         reg.state[i].x = 0;
         reg.state[i].y = 0;
      }
   }

   free(masks);
   free(probabilities);

   return imeas;
}

void ds_scnot(ds_Register reg, int qcont, int qtarg, int time)
{
   int i, j, k;
   ds_Complex z;

   ds_update(reg, qcont, qtarg);

   for (i=0; i<reg.nc/4; i++) {
      ds_swap_indices(i, qcont, qtarg, &j, &k);
      z = reg.state[j];
      reg.state[j] = reg.state[k];
      reg.state[k] = z;
   }

   for (i=0; i<reg.nc/4; i++) {
      ds_controlled_indices(i, qcont, qtarg, &j, &k);
      /* Note that the ds_toffoli gate does not work with pseudo ds_cnot gates
         ie : you cant use ds_esigy(reg.state+j, reg.state+k, Pi) */
      z = reg.state[j];
      reg.state[j] = reg.state[k];
      reg.state[k] = z;
   }

   ds_lerr(reg, qcont, time);
   reg.steps[qcont]+=time;
   ds_lerr(reg, qtarg, time);
   reg.steps[qtarg]+=time;
}

/*------------set/ds_measure----------------------------*/

/* nq2m : number of qubits to measure
   lq2m : list of qubits to measure (sorted from largest to smallest)
   val  : value you wish to observe */

double ds_set_measure(ds_Register reg, int nq2m, int *lq2m, int val)
{
  int i, j, k, imeas, spread, mask, *masks;
  double p, *probabilities, norm;

  /* probably not necessary but simplifies things */
  ds_global_update(reg);

  /* create array of single bit masks for reducing array index */
  masks = (int *)calloc(nq2m, sizeof(int));
  for (i=0; i<nq2m; i++) {
    masks[i] = 1;
    masks[i] <<= lq2m[i];
  }

  /* create array of possible measured values */
  probabilities = (double *)calloc(1 << nq2m, sizeof(double));

  for (i=0; i<reg.nc; i++) {
    /* reduce array index i */
    j = 0;
    for (k=0; k<nq2m; k++) {
      j <<= 1;
      if (i & masks[k]) {
	j++;
      }
    }
    probabilities[j] += ds_modsq(reg.state[i]);
  }

  p = probabilities[val];
  imeas = val;

  if (probabilities[imeas] == 0){
    norm = 0;}
  else{
    norm = 1/sqrt(probabilities[imeas]);}

  mask = 0;
  for (i=0; i<nq2m; i++) mask += masks[i];

  spread = 0;
  i = 1;
  for (j =0; j<nq2m; j++) {
    if (imeas & i) spread += masks[nq2m-1-j];
    i <<= 1;
  }

  for (i=0; i<reg.nc; i++) {
    if ((i & mask) == spread) {
      reg.state[i].x *= norm;
      reg.state[i].y *= norm;
    }
    else {
      reg.state[i].x = 0;
      reg.state[i].y = 0;
    }
  }

  free(masks);
  free(probabilities);
  
  return p;
}

void ds_linearmesh(ds_Register reg,int tqmesh,int bqmesh,int inverse)
{
  
   int terms,i,j;

   terms = ((tqmesh-bqmesh)+1)/2;

  
   if (inverse == 0){
    for (i=0; i<terms-1; i++){
      for(j=0; j<terms-1-i; j++){
	ds_swap(reg,bqmesh+terms+j-i,bqmesh+terms-1+j-i,1);
      }
    }
   }
 

   else{
    for (i=terms-2; i>=0; i--){
      for(j=terms-2-i; j>=0; j--){
	ds_swap(reg,bqmesh+terms+j-i,bqmesh+terms-1+j-i,1);
      }
    }
   }  
  
}

void ds_swapupdown(ds_Register reg,int qstart,int number,int up)
{
   int i;

   for (i=0; i<number; i++) ds_swap(reg,qstart+up*i,qstart+up*i+up,1);

}

void ds_sqrtnot(ds_Register reg,int qcont,int qtarg,int dagger, int time)
{
   int i,j,k;
   ds_Complex angle1,angle2;
   ds_Complex z1,z2;

   ds_update(reg, qcont, qtarg);

   angle1 = (ds_eitheta((1-2*dagger)*ds_Pio4));
   angle2 = (ds_eitheta((2*dagger-1)*ds_Pio4));

   for (i=0; i<reg.nc/4; i++) {
      ds_controlled_indices(i, qcont, qtarg, &j, &k);

      z1 = reg.state[j];
      z2 = reg.state[k];

      reg.state[j].x = ds_root2_2*(ds_multiply(angle1, z1).x + ds_multiply(angle2, z2).x);
      reg.state[j].y = ds_root2_2*(ds_multiply(angle1, z1).y + ds_multiply(angle2, z2).y);
      reg.state[k].x = ds_root2_2*(ds_multiply(angle1, z2).x + ds_multiply(angle2, z1).x);
      reg.state[k].y = ds_root2_2*(ds_multiply(angle1, z2).y + ds_multiply(angle2, z1).y); 
   }

   ds_lerr(reg, qcont, time);
   reg.steps[qcont]+=time;
   ds_lerr(reg, qtarg, time);
   reg.steps[qtarg]+=time;
}


void sqrtnotswap(ds_Register reg,int qcont,int qtarg,int dagger, int time)
{
   int i,j,k;
   ds_Complex z1,z2,z,angle1,angle2;
   
   ds_update(reg, qcont, qtarg);

   angle1 = ds_eitheta((1-2*dagger)*ds_Pio4);
   angle2 = ds_eitheta((2*dagger-1)*ds_Pio4);		    

   for (i=0; i<reg.nc/4; i++) {
      ds_controlled_indices(i, qcont, qtarg, &j, &k);

      z1 = reg.state[j];
      z2 = reg.state[k];

      reg.state[j].x = ds_root2_2*(ds_multiply(angle1, z1).x + ds_multiply(angle2, z2).x);
      reg.state[j].y = ds_root2_2*(ds_multiply(angle1, z1).y + ds_multiply(angle2, z2).y);
      reg.state[k].x = ds_root2_2*(ds_multiply(angle1, z2).x + ds_multiply(angle2, z1).x);
      reg.state[k].y = ds_root2_2*(ds_multiply(angle1, z2).y + ds_multiply(angle2, z1).y); 
   }

   for (i=0; i<reg.nc/4; i++) {
      ds_swap_indices(i, qcont, qtarg, &j, &k);
      z = reg.state[j];
      reg.state[j] = reg.state[k];
      reg.state[k] = z;
   } 

   ds_lerr(reg, qcont, time);
   reg.steps[qcont]+=time;
   ds_lerr(reg, qtarg, time);
   reg.steps[qtarg]+=time;
}


void ds_conswaptoffoli(ds_Register reg,int qcont,int qtarg1,int qtarg2)
{

   ds_cnot(reg,qtarg2,qtarg1,1);
   ds_sqrtnot(reg,qtarg1,qtarg2,0,1);
   ds_cnot(reg,qcont,qtarg1,1);
   ds_sqrtnot(reg,qtarg1,qtarg2,1,1);
   ds_cnots(reg,qcont,qtarg1,1);
   sqrtnotswap(reg,qtarg1,qtarg2,0,1);
   ds_cnot(reg,qtarg1,qcont,1);

}

void ds_linearconswap(ds_Register reg,int N,int tq,int bq)
{

   int terms,i;

   terms = (tq-bq);

   ds_linearmesh(reg,bq+1+2*terms,bq+2,0);
  
   for (i=0; i<terms; i++) ds_conswaptoffoli(reg,N-2-2*i,N-3-2*i,N-4-2*i);

   ds_linearmesh(reg,bq+2+2*terms,bq+3,1);
   
}

/*---------Linear Shor's Algorithm addendum by Simon Devitt----------*/

void ds_linearqft(ds_Register reg,int N,int inverse,int tq,int bq,int maxangle)
{

   int j,i,terms,c,angle;
   terms = (tq-bq);

   if (inverse == 0) c = 1;
   else c = -1;
 
   ds_Hadamard(reg,tq,0);

   for (j=1; j<=terms; j++){

     if ((j % 2) == 0){
       ds_Hadamard(reg,tq,0);

       for (i=j; i>=2; i=i-2){
	 angle = pow(2,i);
	   if (angle > maxangle) ds_swap(reg,tq+1-i,tq-i,1);
	   else ds_cphases(reg,(tq+1-i),(tq-i),(c*ds_Pi/angle),1); 
       }
     }
     else{

       for (i=j; i>=1; i=i-2){
	 angle = pow(2,i);
	   if (angle > maxangle) ds_swap(reg,tq+1-i,tq-i,1);
	   else ds_cphases(reg,(tq+1-i),(tq-i),(c*ds_Pi/angle),1);
       }
     }
   }

   for (j=terms-1; j>=1; j--){

     if ((j % 2) == 0){
       ds_Hadamard(reg,tq,0);

       for(i=j; i>=2; i=i-2){
	 angle = pow(2,i);
	   if (angle > maxangle) ds_swap(reg,tq+1-i,tq-i,1);
	   else ds_cphases(reg,(tq+1-i),(tq-i),(c*ds_Pi/angle),1);
       }
     }
     else{

       for (i=j; i>=1; i=i-2){
	 angle = pow(2,i);
	   if (angle > maxangle) ds_swap(reg,tq+1-i,tq-i,1);
	   else ds_cphases(reg,(tq+1-i),(tq-i),(c*ds_Pi/angle),1);
       }
     }
   }
   
   ds_Hadamard(reg,tq,0);
   
}

void ds_adder(ds_Register reg,int initial,int N,int inverse,int tq,int bq,int carry, int time)
{
   int j,c,mask,factor;
   double angle;
   double totphase;

   if (inverse == 0) c = 1;
   else c = -1;

   angle = ds_Pi/pow(2,(tq-bq));
   mask = 1;
   mask <<= tq-bq+1;
   mask = mask - 1;

   for (j=tq-bq; j>=0; j--){

	factor = mask & initial;
	totphase = c*angle*factor; 

	if (carry == 1) ds_phases(reg,tq-j+1,tq-j,totphase,1);
	else ds_phase(reg,tq-j,totphase,time);
	printf("%d:%g\n",tq-j,totphase);
	initial <<= 1;
   }
	
}

void ds_contadder(ds_Register reg,int initial,int N,int inverse,int tq,int bq, int up)
{
   int j,c,mask,factor,mask2;
   double angle,totphase;
     
   if (inverse == 0) c = 1;
   else c = -1;

   angle = ds_Pi/pow(2,(tq-bq));
   mask = 1;
   mask <<= tq-bq+1;
   mask = mask - 1;

   if (up == 1){

     for (j=0; j<=tq-bq; j++){

	factor = mask & initial;
	totphase = c*angle*factor;

	ds_cphases(reg,bq+j-1,bq+j,totphase,1);
	
	initial <<= 1;
     }
   }

   else{

       for (j=tq-bq; j>=0; j--){

	mask2 = initial;
	mask2 <<= j;

	factor = mask & mask2;
	totphase = c*angle*factor;

	ds_cphases(reg,bq+j+1,bq+j,totphase,1);

       }
   }
}

void ds_addmod1(ds_Register reg,int N,int a,int mod,int tq,int bq,int maxangle)
{

   int terms, j;

   terms = tq;
  
   ds_contadder(reg,a,N,0,tq,bq,0);
   ds_adder(reg,mod,N,1,tq+1,bq+1,0,0);
   ds_linearqft(reg,N,1,tq+1,bq+1,maxangle);
   ds_swap(reg,tq+4,tq+3,1);
   ds_swap(reg,tq+2,tq+3,1);
   ds_cnots(reg,tq+1,tq+2,1);
   
   for (j=1; j<=terms; j++) ds_swap(reg,tq+2-j,tq+1-j,1);

   ds_linearqft(reg,N,0,tq+2,bq+2,maxangle);

   ds_contadder(reg,mod,N,0,tq+2,bq+2,1);

   ds_contadder(reg,a,N,1,tq+1,bq+1,1);

   ds_linearqft(reg,N,1,tq,bq,maxangle);

   for (j=0; j<=terms; j++) ds_swap(reg,tq+1-j,tq-j,1);

   ds_X(reg,tq+1,0);
   ds_cnot(reg,tq+1,tq+2,1);
   ds_X(reg,tq+1,0);

   ds_linearqft(reg,N,0,tq+1,bq+1,maxangle);

   ds_contadder(reg,a,N,0,tq+1,bq+1,1);

   ds_swap(reg,tq+2,tq+3,1);
   ds_swap(reg,tq+3,tq+4,1);
 

}

void ds_addmod1inv(ds_Register reg,int N,int a,int mod,int tq,int bq,int maxangle)
{
  
   int terms, j;

   terms = tq;

  
   ds_swap(reg,tq+3,tq+4,1);
   ds_swap(reg,tq+2,tq+3,1);

   ds_contadder(reg,a,N,1,tq,bq,0);
  
   ds_linearqft(reg,N,1,tq+1,bq+1,maxangle);

   /* make into one gate*/

   ds_X(reg,tq+1,0);
   ds_cnot(reg,tq+1,tq+2,1);
   ds_X(reg,tq+1,0);

   for (j=terms; j>=0; j--) ds_swap(reg,tq+1-j,tq-j,1);

   ds_linearqft(reg,N,0,tq,bq,maxangle);
   ds_contadder(reg,a,N,0,tq,bq,0);
   ds_contadder(reg,mod,N,1,tq+1,bq+1,0);
   ds_linearqft(reg,N,1,tq+2,bq+2,maxangle);
  
   for (j=terms; j>=0; j--) ds_swap(reg,tq+2-j,tq+1-j,1);

   ds_scnot(reg,tq+1,tq+2,1);
   ds_swap(reg,tq+2,tq+3,1);
   ds_swap(reg,tq+4,tq+3,1);
   ds_linearqft(reg,N,0,tq+1,bq+1,maxangle);
   ds_adder(reg,mod,N,0,tq+1,bq+1,0,0);
   ds_contadder(reg,a,N,1,tq+1,bq+1,1);
  
  
}

void ds_conmod1(ds_Register reg,int N,int a,int mod,int tq,int bq,int maxangle)
{

   int i,j,amod;

   for (i=0; i<tq; i++){

     amod = (pow(2.0,i)*a);
     amod = amod % mod;

     for (j=i; j>=0; j--) ds_swap(reg,tq+3+j,tq+4+j,1);

     if (i == 0) ds_sqrtnot(reg,tq+2,tq+1,0,1);
     ds_cnot(reg,tq+3,tq+2,1);
     ds_sqrtnot(reg,tq+2,tq+1,1,1);
     ds_cnots(reg,tq+3,tq+2,1);
     ds_sqrtnot(reg,tq+2,tq+1,0,1);

     ds_addmod1(reg,N,amod,mod,tq,bq,maxangle);
     
     
     ds_sqrtnot(reg,tq+2,tq+1,1,1);
     ds_scnot(reg,tq+3,tq+2,1);
     ds_sqrtnot(reg,tq+2,tq+1,0,1);
     ds_cnot(reg,tq+3,tq+2,1);
     if (i == tq-1) ds_sqrtnot(reg,tq+2,tq+1,1,1);
           
     for (j=0; j<=i; j++) ds_swap(reg,tq+3+j,tq+4+j,1);

   }
}

void conmod1inv(ds_Register reg,int N,int a,int mod,int tq,int bq,int maxangle)
{

   int i,j,amod;

   for (i=0; i<tq; i++){

     amod = (pow(2.0,i)*a);
     amod = amod % mod;

     for (j=i; j>=0; j--) ds_swap(reg,tq+3+j,tq+4+j,1);

     if (i == 0) ds_sqrtnot(reg,tq+2,tq+1,1,1);
     ds_scnot(reg,tq+3,tq+2,1);
     ds_sqrtnot(reg,tq+2,tq+1,0,1);
     ds_cnot(reg,tq+3,tq+2,1);
     ds_sqrtnot(reg,tq+2,tq+1,1,1);

     ds_addmod1inv(reg,N,amod,mod,tq,bq,maxangle);
     
     ds_sqrtnot(reg,tq+2,tq+1,0,1);
     ds_cnot(reg,tq+3,tq+2,1);
     ds_sqrtnot(reg,tq+2,tq+1,1,1);
     ds_cnots(reg,tq+3,tq+2,1);
     if (i == tq-1) ds_sqrtnot(reg,tq+2,tq+1,0,1);

     for (j=0; j<=i; j++) ds_swap(reg,tq+3+j,tq+4+j,1);
   }
}



void ds_ugate1(ds_Register reg,int N,int a,int mod,int ainv,int tq,int bq,int maxangle)
{
  
   ds_conmod1(reg,N,a,mod,tq,bq,maxangle);
   ds_linearqft(reg,N,1,tq,bq,maxangle); 

   ds_swapupdown(reg,tq, tq-bq,-1);
   ds_swapupdown(reg,tq+1,tq-bq,-1);
   ds_swapupdown(reg,tq+3,tq-bq,1);
   ds_swapupdown(reg,tq+2,tq-bq,1);
   ds_linearconswap(reg,N,tq,bq);
   ds_swapupdown(reg,N-1,tq-bq,-1);
   ds_swapupdown(reg,bq+2,tq-bq,1);
   ds_swapupdown(reg,bq+1,tq-bq,1);
   ds_swapupdown(reg,bq,tq-bq,1);

   ds_linearqft(reg,N,0,tq,bq,maxangle);
   conmod1inv(reg,N,ainv,mod,tq,bq,maxangle);
   
}

void ds_shortrick(ds_Register reg,int N,int a,int mod,int tq,int bq,int j,FILE *out,int maxangle){

   int i,y,k,m,d,amod,terms,mask,value,mask2;
   double totalphase;
   float prob;
   int measarr[1];

   terms = 2*tq;
   measarr[0] = tq+2;
   d = pow(2,terms-1);
   prob = 1;
   mask = 1;
   mask2 = 0;

   ds_linearqft(reg,N,0,tq,bq,maxangle);

   for (i=terms-1; i>=0; i--){

     value = mask & j;
     if (value != 0) value = 1;
     totalphase = ((j & mask2)*ds_Pi)/mask;

     ds_Hadamard(reg,tq+2,0);
     ds_countstate(reg);
 
     amod = (a % mod);
     for (y=1; y<d; y++) amod = ((a*amod) % mod);
     if (amod == 1) goto cont1;

     for (k=0; k<mod; k++){
       m = (1-amod*k);
      
       if (m % mod == 0){	
	 ds_ugate1(reg,N,amod,mod,k,tq,bq,maxangle);
	 goto cont1;
       }
     }

   cont1:

     ds_phase(reg,tq+2,totalphase,0);
     ds_Hadamard(reg,tq+2,0);
     ds_countstate(reg);
    
     prob = prob*ds_set_measure(reg,1,measarr,value);

     if (prob == 0) goto cont3;
     if (value == 1) ds_X(reg,tq+2,0);  
  
     mask2 = mask2 + mask;
     mask <<= 1; 
     d = d/2;
    
   }

   ds_linearqft(reg,N,1,tq,bq,maxangle);

 cont3:

   fprintf(out, "%d: %d: %d: %f \n",mod,a,j,prob);
    
}


/* ----------------General Case (i.e arbitrary Interactions) addendum by Simon Devitt---------*/

void ds_toffoli(ds_Register reg,int qcont1,int qcont2,int qtarg)
{
   ds_sqrtnot(reg,qcont2,qtarg,0,1);
   ds_global_update(reg);

   ds_cnot(reg,qcont1,qcont2,1);
   ds_global_update(reg);

   ds_sqrtnot(reg,qcont2,qtarg,1,1);
   ds_global_update(reg);

   ds_cnot(reg,qcont1,qcont2,1);
   ds_global_update(reg);

   ds_sqrtnot(reg,qcont1,qtarg,0,1);
   ds_global_update(reg);
}
 
void ds_generalqft(ds_Register reg,int N,int tq,int bq,int maxangle)
{

   int i;
   int j;
   double angle;

   ds_Hadamard(reg,tq,0);

   for (i=tq; i>bq; i= i-1) { 

     for (j=tq; j>i-1; j=j-1 ){

       angle = pow(2.0,j-i+1);

       if (angle < maxangle){ 
	 ds_cphase(reg,j,i-1,ds_Pi/angle,1);
	 ds_global_update(reg);
       }
     }
   
     ds_Hadamard(reg,i-1,0);

   }
}

void ds_generalqftinv(ds_Register reg,int N,int tq,int bq,int maxangle){

   int i;
   int j;
   double angle;

   ds_Hadamard(reg,bq,0);

   for (i=bq; i<tq; i++){

     for (j=i+1; j<=tq; j++){

	 angle = pow(2.0,j-i);

	 if (angle < maxangle){
	   ds_cphase(reg,j,i,-ds_Pi/angle,1);
	   ds_global_update(reg);
	 }

     }

     ds_Hadamard(reg,i+1,0);

   }
}
  
void ds_addergen(ds_Register reg,int N,int initial,int inverse,int tq,int bq)
{
   int j,c,mask,factor;
   double angle;
   double totphase;

   if (inverse == 0) c = 1;
   else c = -1;

   angle = ds_Pi/pow(2,(tq-bq));
   mask = 1;
   mask <<= tq-bq+1;
   mask = mask - 1;

   for (j=0; j<=tq-bq; j++){

	factor = mask & initial;
	totphase = c*angle*factor; 

        ds_phase(reg,tq-j,totphase,0);

	initial <<= 1;
   }

}

void ds_contaddergen(ds_Register reg,int N,int initial,int inverse,int tq,int bq,int control)
{
   int j,c,mask,factor,mask2;
   double angle,totphase;
     
   if (inverse == 0) c = 1;
   else c = -1;

   angle = ds_Pi/pow(2,(tq-bq));
   mask = 1;
   mask <<= tq-bq+1;
   mask = mask - 1;

      for (j=0; j<=tq-bq; j++){

	mask2 = initial;
	mask2 <<= tq-bq-j;

	factor = mask & mask2;
	totphase = c*angle*factor;

	ds_cphase(reg,control,bq+j,totphase,1);
	ds_global_update(reg);

      }
}

void ds_addmodgen(ds_Register reg,int N,int a,int mod,int inverse,int tq,int bq,int ancilla,int control,int maxangle)
{

   if (inverse == 0){

     ds_contaddergen(reg,N,a,0,tq,bq,control);
     ds_addergen(reg,N,mod,1,tq,bq);
     ds_generalqftinv(reg,N,tq,bq,maxangle);
     ds_cnot(reg,tq,ancilla,1);
     ds_global_update(reg);
     ds_generalqft(reg,N,tq,bq,maxangle);
     ds_contaddergen(reg,N,mod,0,tq,bq,ancilla);
     ds_contaddergen(reg,N,a,1,tq,bq,control);
     ds_generalqftinv(reg,N,tq,bq,maxangle);

     ds_X(reg,tq,0);
     ds_cnot(reg,tq,ancilla,1);
     ds_X(reg,tq,0);
     ds_global_update(reg);

     ds_generalqft(reg,N,tq,bq,maxangle);
     ds_contaddergen(reg,N,a,0,tq,bq,control);
   }

   else{

     ds_contaddergen(reg,N,a,1,tq,bq,control);
     ds_generalqftinv(reg,N,tq,bq,maxangle);

     ds_X(reg,tq,0);
     ds_cnot(reg,tq,ancilla,1);
     ds_X(reg,tq,0);
     ds_global_update(reg);

     ds_generalqft(reg,N,tq,bq,maxangle);
     ds_contaddergen(reg,N,a,0,tq,bq,control);
     ds_contaddergen(reg,N,mod,1,tq,bq,ancilla);
     ds_generalqftinv(reg,N,tq,bq,maxangle);
     ds_cnot(reg,tq,ancilla,1);
     ds_global_update(reg);
     ds_generalqft(reg,N,tq,bq,maxangle);
     ds_addergen(reg,N,mod,0,tq,bq);
     ds_contaddergen(reg,N,a,1,tq,bq,control);
   }
}

void ds_contmultgen(ds_Register reg,int N,int a,int mod,int inverse,int tq,int bq,int control,int maxangle)
{
   int i,factor;

   factor = (a % mod);

   for (i = 1; i<=tq; i++){

     ds_toffoli(reg,control,tq+3+i,tq+1);
       
     ds_addmodgen(reg,N,factor,mod,inverse,tq,bq,tq+3,tq+1,maxangle);

     ds_toffoli(reg,control,tq+3+i,tq+1);

     factor = (2*factor) % mod;
   }

}

void ds_cswap(ds_Register reg,int control,int qtarg1,int qtarg2)
{

   ds_cnot(reg,qtarg1,qtarg2,1);
   ds_global_update(reg);

   ds_toffoli(reg,control,qtarg2,qtarg1);

   ds_cnot(reg,qtarg1,qtarg2,1);
   ds_global_update(reg);

}

void ds_contswap(ds_Register reg,int tq,int bq,int control)
{

   int i;

   for (i = 0; i<tq; i++){

     ds_cswap(reg,control,bq+i,tq+4+i);
     ds_global_update(reg);
   
   }

}


void ds_ugategen(ds_Register reg,int N,int a,int mod,int ainv,int tq,int bq,int control,int maxangle)
{    
    ds_contmultgen(reg,N,a,mod,0,tq,bq,control,maxangle);
    ds_generalqftinv(reg,N,tq,bq,maxangle);
    ds_contswap(reg,tq,bq,control);
    ds_generalqft(reg,N,tq,bq,maxangle);
    ds_contmultgen(reg,N,ainv,mod,1,tq,bq,control,maxangle);
}

void ds_shortrickgen(ds_Register reg,int N,int a,int mod,int tq,int bq,int maxangle,int j,FILE *out)
{
   int i,y,k,m,d,amod,terms,mask,value,mask2;
   double totalphase;
   float prob;
   int measarr[1];

   terms = 2*tq;
   measarr[0] = tq+2;
   d = pow(2,terms-1);
   prob = 1;
   mask = 1;
   mask2 = 0;

   ds_generalqft(reg,N,tq,bq,maxangle);
  
   for (i=terms-1; i>=0; i--){

     value = mask & j;
     if (value != 0) value = 1;
     totalphase = ((j & mask2)*ds_Pi)/mask;

     ds_Hadamard(reg,tq+2,0);
 
     amod = (a % mod);
     for (y=1; y<d; y++) amod = ((a*amod) % mod);
     if (amod == 1) goto cont1;

     for (k=0; k<mod; k++){
       m = (1-amod*k);
      
       if (m % mod == 0){	
	 ds_ugategen(reg,N,amod,mod,k,tq,bq,tq+2,maxangle);
	 goto cont1;
       }
     }

   cont1:

     ds_phase(reg,tq+2,totalphase,0);

     ds_Hadamard(reg,tq+2,0);
    
     prob = prob*ds_set_measure(reg,1,measarr,value);

     if (prob == 0) goto cont3;
     if (value == 1) ds_X(reg,tq+2,1);  
  
     mask2 = mask2 + mask;
     mask <<= 1; 
     d = d/2;
    
   }

   ds_generalqftinv(reg,N,tq,bq,maxangle);

 cont3:

   fprintf(out, "%d: %d: %d: %f \n",mod,a,j,prob);
    
}

/*-----------------------------otherstuff---------------------------------*/

void ds_printmod(ds_Register reg,FILE *out)
{
   int i;

   for (i=0; i<reg.nc; i++){

      if (ds_modsq(reg.state[i]) > 1e-21 ) fprintf(out, "%d: %g \n", i, ds_modsq(reg.state[i]));

   }

   fprintf(out,"\n");
}

void ds_countstate(ds_Register reg)
{

  int i, count;

  count = 0;
  for (i=0; i<reg.nc; i++){

    if (ds_modsq(reg.state[i]) > 1e-21) count = count + 1;

  }

  printf("%d\n", count);
}
