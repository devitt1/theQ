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
	 ds_X(reg, q, 0);
      }
      else if (p < 2*reg.err/3) {
	 ds_Z(reg, q, 0);
      }
      else if (p < reg.err) {
	 ds_XZ(reg, q, 0);
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

/* Kane compatible if |q-x| = 1 */

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
