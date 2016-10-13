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

	for (int i = 0; i < 100000; i++) 
	{
		taga1991.dump();
	
		if(!taga1991.next() ){ 
			printf("\nsingular matrix !!!  [i] : %d\n",i); 
			exit(0); 
		}
	}
	return 0;
}

