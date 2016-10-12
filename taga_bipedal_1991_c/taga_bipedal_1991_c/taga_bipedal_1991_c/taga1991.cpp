// G. Taga, et al., 
// Self-organized control of bipedal locomotion by neural oscillators in unpredictable enviromnent, 
// Biological Cyberetics, 65, 147-150, 1991
//
// Masahiro Furukawa
// m.furukawa@ist.osaka-u.ac.jp
//
// rev 0.0, Created,  Oct 4, 2016-

#include "taga1991.h"


void draw() {
	// (6) derivered from Matsuoka 1985,1987 and Grillner 1981.

	for (int i = 1; i <= 12; i++) {
		y[i] = f(u[i]);

		print(u[i] + "\t");

		u_d[i] = -u[i] - beta*v[i] + u[0] + Feed[i];
		v_d[i] = -v[i] + y[i];

		for (int j = 1; j <= 12; j++)
			u_d[i] = u_d[i] + w[i][j] * y[j];

		u_d[i] = u_d[i] / tau[i];
		v_d[i] = u_d[i] / taud[i];
	}



}