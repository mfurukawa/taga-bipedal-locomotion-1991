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
#include <stdlib.h>
#include "gauss_jordan.h"
#include "taga1991.h"

int main() {

	Taga1991 taga1991;
	double t = 0.0;

	for (int i = 0; i < 100000; i++) 
	{
		if(!taga1991.next() ){ 
			printf("\nsingular matrix !!!  [i] : %d\n",i); 
			exit(0); 
		}
		printf("% 1.4lf\t% 3.2lf\t% 3.2lf\t% 3.2lf\t% 3.2lf\t% 3.2lf\t% 3.2lf\t% 3.2lf\n", 
			   t, 
			   taga1991.u[1],
			   taga1991.u[2],
			   taga1991.u[3],
			   taga1991.u[4],
			   taga1991.u[5],
			   taga1991.v[1],
			   taga1991.v[2]);
		t += taga1991.dt;
	}
	return 0;
}

