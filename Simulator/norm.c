#include <math.h>
#include "norm.h"
#include <stdlib.h>

#include <stdio.h>

//int ds_seed = 1, ds_b = 314159221;
//double ds_m = 2.147483648e9;

void ds_srand32(long s)
{
   //ds_seed = s;
   srand48(s);
}

/* Generates a pseudo-random integer from -(2^31) to 2^31-1 */
int ds_rand32()
{
   //ds_seed = ds_seed*ds_b+1;
   //return ds_seed;
   return lrand48();
}

double ds_uniform()
{
   //return ((double)ds_rand32()/ds_m+1)/2;

   //double u_rand = (double)random() / (double)RAND_MAX;
   double u_rand = drand48();

   //FILE *out=fopen("uniform.txt","a");

   //printf("%f\n",u_rand); 
   return u_rand;
}

double ds_norm()
{
   int i, j;
   double x1, x2, r2;

   i = ds_rand32();
   j = ds_rand32();
   x1 = (double)i/RAND_MAX;
   x2 = (double)j/RAND_MAX;
   r2 = x1*x1+x2*x2;
   while (r2 >= 1 || r2 == 0) {
      i = ds_rand32();
      j = ds_rand32();
      x1 = (double)i/RAND_MAX;
      x2 = (double)j/RAND_MAX;
      r2 = x1*x1+x2*x2;
   }

   return x1*sqrt(-2*log(r2)/r2);
}
