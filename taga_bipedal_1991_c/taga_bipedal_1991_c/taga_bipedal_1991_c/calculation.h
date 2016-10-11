#ifndef __CALCULATION_H_INCLUDED__
#define __CALCULATION_H_INCLUDED__

void Runge_Kutta(double (*f[])(double t, double *x), double t0, double *x, double tn, int div, int num);
#endif /* __CALCULATION_H_INCLUDED__ */
