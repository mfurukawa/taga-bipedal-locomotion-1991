/* (c) Masahiro Furukawa
*
* m.furukawa@ist.osaka-u.ac.jp
*
* Created Oct 11, 2016
*
* Description
* - Runge-Kutta-Gill method (4th order)
*
* referece:
* http://www.edu.cc.uec.ac.jp/mce/NumAnalysis/pdf/NumericalAnalysis09.pdf
*/

#include <stdio.h>
#include "rkg4.h"

int rkg4(double t, double x, double v) {

	double a[3][GAUSS_JORDAN_MAXN + 10], b[3];
	int n, i, j;

	a[1][1] = 2;
	a[1][2] = 1;
	a[2][1] = 4;
	a[2][2] = 3;

	b[1] = 1;
	b[2] = 1;

	n = 2;

	printf("\nA \t\t\t\t| b\n");
	for (int i = 1; i <= n; i++) {
		for (int j = 1; j <= n; j++) {
			printf("%lf\t", a[i][j]);
		}
		printf("| %lf\n",b[i]);
	}

	// obtain solution and inverce matrix
	gauss_jordan(2, a, b);

	printf("\nsolution\n");
	for (int i = 1; i <= n; i++) {
		printf("x%d = %lf \n", i, b[i]);
	}
	
	printf("\ninverce matrix\n");
	for (int i = 1; i <= n; i++) {
		for (int j = 1; j <= n; j++) {
			printf("%lf\t", a[i][j]);
		}
		printf("\n");
	}
	return 0;
}

