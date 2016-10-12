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

	// A. The equations of motion for the bipedal musculo-skeletal system

	// Torques generated at each joint are given by:

	Tr[1] = p_he*y[2] - p_hf*y[1];
	Tr[2] = p_he*y[4] - p_hf*y[3];
	Tr[3] = p_ke*y[6] - p_kf*y[5];
	Tr[4] = p_ke*y[8] - p_kf*y[7];
	Tr[5] = (p_ae*y[10] - p_af*y[9]) *h(Fg[2]);
	Tr[6] = (p_ae*y[12] - p_af*y[11])*h(Fg[4]);
	
	// Feedback pathway

	Feed[1] = a[1] * (x[5] - PI / 2.0) - a[2] * (x[8] - PI / 2.0) + a[3] * (x[11] - PI / 2.0)*h(Fg[2]) + a[4] * h(Fg[4]);
	Feed[2] = a[1] * (PI / 2.0 - x[5]) - a[2] * (PI / 2.0 - x[8]) + a[3] * (PI / 2.0 - x[11])*h(Fg[2]) + a[4] * h(Fg[4]);
	Feed[3] = a[1] * (x[8] - PI / 2.0) - a[2] * (x[5] - PI / 2.0) + a[3] * (x[14] - PI / 2.0)*h(Fg[4]) + a[4] * h(Fg[2]);
	Feed[4] = a[1] * (PI / 2.0 - x[8]) - a[2] * (PI / 2.0 - x[5]) + a[3] * (PI / 2.0 - x[14])*h(Fg[4]) + a[4] * h(Fg[2]);
	Feed[5] = a[5] * (PI / 2.0 - x[14])*h(Fg[4]);
	Feed[6] = a[5] * (x[14] - PI / 2.0)*h(Fg[4]);
	Feed[7] = a[5] * (PI / 2.0 - x[11])*h(Fg[2]);
	Feed[8] = a[5] * (x[11] - PI / 2.0)*h(Fg[2]);
	Feed[9] = a[6] * (PI / 2.0 - x[11])*h(Fg[2]) - a[7] * (PI / 2.0 - x[14])*h(Fg[4]) + a[8] * xd[11] * h(Fg[2]);
	Feed[10] = a[6] * (x[11] - PI / 2.0)*h(Fg[2]) - a[7] * (x[14] - PI / 2.0)*h(Fg[4]) + a[8] * xd[11] * h(Fg[2]);
	Feed[11] = a[6] * (PI / 2.0 - x[14])*h(Fg[2]) - a[7] * (PI / 2.0 - x[11])*h(Fg[4]) + a[8] * xd[14] * h(Fg[4]);
	Feed[12] = a[6] * (x[14] - PI / 2.0)*h(Fg[2]) - a[7] * (x[11] - PI / 2.0)*h(Fg[4]) + a[8] * xd[14] * h(Fg[4]);

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

	// (xr,yr) and (xl, yl) represent the positions of the ankles, which are given by:

	xr = x[9] + (l[2] / 2.0)*cos(x[11]);
	yr = x[10] + (l[2] / 2.0)*sin(x[11]);
	xl = x[12] + (l[2] / 2.0)*cos(x[14]);
	yl = x[13] + (l[2] / 2.0)*sin(x[14]);


	// yg(x) is function which represents the terrain. When the ground is level, yg(x) = 0.
	// Horizontal and vertical forces on the ankles are given by:

	if (yr - yg(xr) < 0) {
		Fg[1] = -kg*(xr - xr0) - bg*xr_d;
		Fg[2] = -kg*(yr - yr0) - bg*f(-yr_d);
	}
	else {
		Fg[1] = 0;
		Fg[2] = 0;
	}

	if (yl - yg(xl) < 0) {
		Fg[3] = -kg*(xl - xl0) - bg*xl_d;
		Fg[4] = -kg*(yl - yl0) - bg*f(-yl_d);
	}
	else {
		Fg[3] = 0;
		Fg[4] = 0;
	}

}