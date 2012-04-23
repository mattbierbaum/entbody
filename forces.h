#ifndef __FORCES_H__
#define __FORCES_H__

void force_hertz(double *dx, double dist, double radi, double radj, int typei, int typej, double *f);
void force_morse(double *dx, double dist, double radi, double radj, int typei, int typej, double *f);

void force_damping(double *v, double *f);
void force_thermal(double T, double *f);
void force_kick(double *k, double *f);
void force_gravity(double *f);

#endif
