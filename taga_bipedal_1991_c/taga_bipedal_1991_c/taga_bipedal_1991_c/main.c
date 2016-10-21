#include<stdio.h>
#include<math.h>
#include"calculation.h"

#define LOOP 10000

double f1(double t, double *x);
double f2(double t, double *x);
double f3(double t, double *x);

int main()
{
	int i;
	double t=0, tn=100;
	double x[] = {1., 1., 1.};
	double (*f[3])(double, double*);
	double h;

	f[0] = f1;
	f[1] = f2;
	f[2] = f3;

	h=(tn - t)/LOOP;

	printf("%f,%f,%f\n", x[0], x[1], x[2]);

	for(i=0; i<LOOP; i++){
		Runge_Kutta(f, t, x, t+h, 1, 3);
		printf("%f,%f,%f\n", x[0], x[1], x[2]);
		t+=h;
	}	

	return 0;
}

double f1(double t, double *x)
{
	return ((x[2] - 0.7)*x[0] - 3.5*x[1]);
}

double f2(double t, double *x)
{
	return (3.5*x[0] + (x[2] - 0.7)*x[1]);
}

double f3(double t, double *x)
{
	double a = .6;

	return (a + x[2] - x[2]*x[2]*x[2]/3 - (x[0]*x[0] + x[1]*x[1])*(1 + 0.25*x[2]));
}


