#include<stdio.h>
#include<math.h>
#include"calculation.h"

void Runge_Kutta(double (*f[])(double t, double *x), double t0, double *x, double tn, int div, int num)
{
	double h, t;
	double *k1, *k2, *k3, *k4, *temp;
	int i, j;

	if((k1 = (double*)malloc(sizeof(double)*num)) == NULL){
		fprintf(stderr, "Allocation Error\n");
		exit(1);
	}
	if((k2 = (double*)malloc(sizeof(double)*num)) == NULL){
		fprintf(stderr, "Allocation Error\n");
		exit(1);
	}
	if((k3 = (double*)malloc(sizeof(double)*num)) == NULL){
		fprintf(stderr, "Allocation Error\n");
		exit(1);
	}
	if((k4 = (double*)malloc(sizeof(double)*num)) == NULL){
		fprintf(stderr, "Allocation Error\n");
		exit(1);
	}
	if((temp = (double*)malloc(sizeof(double)*num)) == NULL){
		fprintf(stderr, "Allocation Error\n");
		exit(1);
	}

	h = (tn - t0);
	if(div) h /= div;

	t = t0;

	for(i=0; i<div; i++){
		for(j=0; j<num; j++){
			k1[j] = (*f[j])(t, x);
			temp[j] = x[j] + h*k1[j]/2;
		}
		for(j=0; j<num; j++){
			k2[j] = (*f[j])(t+h/2, temp);
		}
		for(j=0; j<num; j++){
			temp[j] = x[j] + h*k2[j]/2;
		}
		for(j=0; j<num; j++){
			k3[j] = (*f[j])(t+h/2, temp);
		}
		for(j=0; j<num; j++){
			temp[j] = x[j] + h*k3[j];
		}
		for(j=0; j<num; j++){
			k4[j] = (*f[j])(t+h, temp);
			x[j] += (k1[j] + 2*k2[j] + 2*k3[j] + k4[j])*h/6;
		}
		t += h;
	}

	free(k1);
	free(k2);
	free(k3);
	free(k4);
	free(temp);
}

