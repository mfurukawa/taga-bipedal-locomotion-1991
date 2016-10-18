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
		if(taga1991.x[2]<0){ 
		  fprintf(stderr, "\nfall down x[2] got less than ground level in cycle [i] : %d\n",i); 
		  exit(0); 
		}
		if(t>10.0){ 
		  fprintf(stderr, "\nsimulation time finished[i] : %d\n",i); 
		  exit(0); 
		}

		// time : 1
		printf("% 1.8lf,", t);

		// link position and orientation : 2-15
		for(int i=1; i<=14; i++) printf("% 3.10lf,", taga1991.x[i]);

		// ankle position : 16-23
		printf("% 3.10lf,% 3.10lf,% 3.10lf,% 3.10lf,% 3.10lf,% 3.10lf,% 3.10lf,% 3.10lf,",
			   taga1991.xr0,
			   taga1991.xr,
			   taga1991.yr0,
			   taga1991.yr,
			   taga1991.xl0,
			   taga1991.xl,
			   taga1991.yl0,
			   taga1991.yl);

		// neural rythm generator : 24-35
		for(int i=1; i<=12; i++) printf("% 3.10lf,", taga1991.u[i]);

		printf("\n");

		if(!taga1991.next() ){ 
		  fprintf(stderr, "\nsingular matrix !!!  [i] : %d\n",i); 
		  exit(0); 
		}

		t += taga1991.dt;
	}
	return 0;
}

