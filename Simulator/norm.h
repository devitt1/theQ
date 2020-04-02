#ifndef NORM_H
#define NORM_H

extern int ds_seed, ds_b;
extern double ds_m;

void ds_srand32(long s);
int ds_rand32();
double ds_uniform();
double ds_norm();

#endif
