/* (c) Masahiro Furukawa
*
* m.furukawa@ist.osaka-u.ac.jp
*
* Created Oct 11, 2016
*
* Description
* - gauss-jordan method test
*
* referece:
* http://www.yamamo10.jp/yamamoto/lecture/2006/5E/Linear_eauations/gaussj.pdf
*/

#include <stdio.h>
#include <memory.h>
#include "gauss_jordan.h"

int main() {

	double     a[4][GAUSS_JORDAN_MAXN + 10], b[4];
	double inv_a[4][GAUSS_JORDAN_MAXN + 10];
	int n, i, j;

	a[1][1] = 2.3;
	a[1][2] = 1;
	a[1][3] = 4;

	a[2][1] = 3.1;
	a[2][2] = 3;
	a[2][3] = 6;

	a[3][1] = 3;
	a[3][2] = 5;
	a[3][3] = 8;

	b[1] = 1;
	b[2] = 3;
	b[3] = 1;

	n = 3;

	printf("\nA \t\t\t\t\t\t| b\n");
	for (int i = 1; i <= n; i++) {
		for (int j = 1; j <= n; j++) {
			printf("%lf\t", a[i][j]);
		}
		printf("| %lf\n",b[i]);
	}

	memcpy(&inv_a[0][0], &a[0][0], sizeof(a));

	// obtain solution and inverce matrix
	gauss_jordan(3, inv_a, b);

	printf("\nsolution\n");
	for (int i = 1; i <= n; i++) {
		printf("\tx%d = %lf \n", i, b[i]);
	}
	
	printf("\ninverce matrix\n");
	for (int i = 1; i <= n; i++) {
		for (int j = 1; j <= n; j++) {
			printf("\t%lf", inv_a[i][j]);
		}
		printf("\n");
	}

	printf("\ninv_A * A\n");
	for (int i = 1; i <= n; i++) {
		for (int j = 1; j <= n; j++) {
		  double tmp = 0.0;
		  for (int k = 1; k <= n; k++) 
			tmp += a[i][k] * inv_a[k][j];
		  printf("\t%lf", tmp);
		}
		printf("\n");
	}
	return 0;
}

