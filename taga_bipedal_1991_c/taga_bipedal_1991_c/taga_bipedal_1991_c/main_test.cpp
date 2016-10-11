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
#include "gauss_jordan.h"

int main() {

	double a[3][GAUSS_JORDAN_MAXN + 10], b[3];
	int n, i, j;

	a[1][1] = 2;
	a[1][2] = 1;
	a[2][1] = 4;
	a[2][2] = 3;

	b[1] = 1;
	b[2] = 1;

	n = 2;

	gauss_jordan(2, a, b);

	printf("\nsolution\n");
	for (int i = 1; i <= n; i++) {
		printf("x%d = %lf \n", i, b[i]);
	}
	
	printf("\ninverce matrix\n");
	for (int i = 1; i <= n; i++) {
		for (int j = 1; j <= n; j++) {
			printf("%lf\t", i, a[i][j]);
		}
		printf("\n");
	}
	return 0;
}

