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
		if(t>10){ 
		  exit(0); 
		}
		if(!taga1991.next() ){ 
		  printf("\nsingular matrix !!!  [i] : %d\n",i); 
		  exit(0); 
		}
		printf("% 1.4lf,% 3.10lf,% 3.10lf,% 3.10lf,% 3.10lf,% 3.10lf,% 3.10lf,% 3.10lf,% 3.10lf,% 3.10lf,% 3.10lf,% 3.10lf,% 3.10lf,% 3.10lf,% 3.10lf,% 3.10lf,% 3.10lf,% 3.10lf,% 3.10lf,% 3.10lf,% 3.10lf,% 3.10lf,% 3.10lf,% 3.10lf,% 3.10lf\n", 
			   t, 
			   taga1991.u[1],
			   taga1991.u[2],
			   taga1991.u[3],
			   taga1991.u[4],
			   taga1991.u[5],
			   taga1991.u[6],
			   taga1991.u[7],
			   taga1991.u[8],
			   taga1991.u[9],
			   taga1991.u[10],
			   taga1991.u[11],
			   taga1991.u[12],
			   taga1991.v[1],
			   taga1991.v[2],
			   taga1991.v[3],
			   taga1991.v[4],
			   taga1991.v[5],
			   taga1991.v[6],
			   taga1991.v[7],
			   taga1991.v[8],
			   taga1991.v[9],
			   taga1991.v[10],
			   taga1991.v[11],
			   taga1991.v[12]);
		t += taga1991.dt;
	}
	return 0;
}

