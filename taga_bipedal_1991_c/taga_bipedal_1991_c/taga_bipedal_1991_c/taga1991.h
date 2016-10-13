#pragma once
// G. Taga, et al., 
// Self-organized control of bipedal locomotion by neural oscillators in unpredictable enviromnent, 
// Biological Cyberetics, 65, 147-150, 1991
//
// Masahiro Furukawa
// m.furukawa@ist.osaka-u.ac.jp
//
// rev 0.0, Created,  Oct 12, 2016-

#include <stdio.h>
#include <math.h>

#include "gauss_jordan.h"

#define a1 a[1]
#define a2 a[2]
#define a3 a[3]
#define a4 a[4]
#define a5 a[5]
#define a6 a[6]
#define a7 a[7]
#define a8 a[8]

#define Feed1 Feed[1]
#define Feed2 Feed[2]
#define Feed3 Feed[3]
#define Feed4 Feed[4]
#define Feed5 Feed[5]
#define Feed6 Feed[6]
#define Feed7 Feed[7]
#define Feed8 Feed[8]
#define Feed9 Feed[9]
#define Feed10 Feed[10]
#define Feed11 Feed[11]
#define Feed12 Feed[12]

#define x1 x[1]
#define x2 x[2]
#define x3 x[3]
#define x4 x[4]
#define x5 x[5]
#define x6 x[6]
#define x7 x[7]
#define x8 x[8]
#define x9 x[9]
#define x10 x[10]
#define x11 x[11]
#define x12 x[12]
#define x13 x[13]
#define x14 x[14]

#define xd5 xd[5]
#define xd8 xd[8]
#define xd11 xd[11]
#define xd14 xd[14]

#define xd52 xd[5]*xd[5]
#define xd82 xd[8]*xd[8]
#define xd112 xd[11]*xd[11]
#define xd142 xd[14]*xd[14]

#define s5  sin(x5)
#define s8  sin(x8)
#define s11 sin(x11)
#define s14 sin(x14)

#define c5  cos(x5)
#define c8  cos(x8)
#define c11 cos(x11)
#define c14 cos(x14)

class Taga1991 
{
private:

	// A. The equations of motion for the bipedal musculo-skeletal system

	double Fg1, Fg2, Fg3, Fg4;				// Horizontal and vertical forces on the ankles. See Fig. 12. 
	double Tr1, Tr2, Tr3, Tr4, Tr5, Tr6;	// Torque

	double x  [15];
	double xd [15];
	double xdd[15];
	double y  [15]; // the output of the i-th neuron. (6)

	// Newton-Eular method
	double P[15][9]; // position
	double Q[15];    // feedback etc
	double C[9][15]; // constraint
	double D[9];     // constraint
	double CQ[9]; 
	double DCQ[9]; 
	double CP[9][9];
	double Pinv_CP[15][9]; 
	double inv_CP[8][GAUSS_JORDAN_MAXN + 10], b[8];
	
	double xr, yr, xl, yl, xr0, yr0, xl0, yl0 = 0;
	double xr_d, yr_d, xl_d, yl_d = 0;

	double u   [15]; // the inner state of the i-th neuron. u0 is an external input with a constant rate
	double ud  [15];
	double v   [15]; // a variable represeinting the degree of the adaptation or self-inhibition effect of the i-th neuron 
	double vd  [15];
	double tau [13]; // time constants of the inner state
	double taud[13]; // the adaptation effect


	// C. Feedback pathway

	// Feedback signals from the musculo-skeletal system to the neural rhythm generetor are given by:

	double Feed[13];

	// D. Simlumation Parameters

	// Musculo-skeletal system

	double m1, m2;
	double l1, l2;
	double I1, I2;
	double b1, b2;
	double M, bk, kk, g, kg, bg, p_hf, p_he, p_kf, p_ke, p_af, p_ae = 75.0;
	
	// neural rhythm generator

	double beta;
	double w[13][13];  // a connecting weight
	double w_fe, w_rl, w_hka;

	// feedback

	double a[9];

	double f(double x) { if (x<0) return 0; else return x; } // max(0,x);
	double h(double x) { if (x<0) return 0; else return 1; }
	double yg(double x) { return 0; }

	// Runge Kutta Method coefficient
	double k1[15][4]; // [4] means twice diferentiated x+ u + v
	double k2[15][4];
	double k3[15][4];
	double k4[15][4];
	// escape current state 
	double u_esc[15], v_esc[15];
	double x_esc[15], xd_esc[15], xdd_esc[15];

	void update(void);

public:
	double dt;

	Taga1991();
	~Taga1991();

	void next(void);
};
