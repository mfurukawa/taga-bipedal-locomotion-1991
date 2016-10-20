/* (c) Masahiro Furukawa
*
* m.furukawa@ist.osaka-u.ac.jp
*
* Created Oct 11, 2016
*
* Description
* - gauss-jordan method to obtain inverce matrix and solution.
*
* referece:
* http://www.yamamo10.jp/yamamoto/lecture/2006/5E/Linear_eauations/gaussj.pdf
*/

#include <math.h>
#include "gauss_jordan.h"

int gauss_jordan(int n, long double a[][GAUSS_JORDAN_MAXN + 10], long double b[]) {

	int ipv, i, j;
	long double inv_pivot, temp;
	long double big;
	int pivot_row, row[GAUSS_JORDAN_MAXN + 10];
	
	for (ipv = 1; ipv <= n; ipv++) {
		
		// search maximum value
		big = 0.0;
		for (i = ipv; i <= n; i++) {
			if (lfabs(a[i][ipv]) > big) {
				big = lfabs(a[i][ipv]);
				pivot_row = i;
			}
		}
		if (big == 0.0) { return 0; } // singular matrix !!!
		row[ipv] = pivot_row;

		// swap row
		if (ipv != pivot_row) { 
			for (i = 1; i <= n; i++) {
				temp = a[ipv][i];
				a[ipv][i] = a[pivot_row][i];
				a[pivot_row][i] = temp;
			}
			temp = b[ipv];
			b[ipv] = b[pivot_row];
			b[pivot_row] = temp;
		}

		// regularize diagnal component
		inv_pivot = 1.0 / a[ipv][ipv];
		a[ipv][ipv] = 1.0;
		for (j = 1; j <= n; j++) {
			a[ipv][j] *= inv_pivot;
		}
		b[ipv] *= inv_pivot;

		// zero out pivot column except pivot row
		for (i = 1; i <= n; i++) {
			if (i != ipv) {
				temp = a[i][ipv];
				a[i][ipv] = 0.0;
				for (j = 1; j <= n; j++) {
					a[i][j] -= temp*a[ipv][j];
				}
				b[i] -= temp*b[ipv];

			}
		}

	}

	// swap column
	for (j = n; j >= 1; j--) {
		if (j != row[j]) {
			for (i = 1; i <= n; i++) {
				temp = a[i][j];
				a[i][j] = a[i][row[j]];
				a[i][row[j]] = temp;
			}
		}
	}

	return 1;
}

long double lfabs(long double x)
{
  if (x > 0.0) return x; 
  else return -x; 
}
