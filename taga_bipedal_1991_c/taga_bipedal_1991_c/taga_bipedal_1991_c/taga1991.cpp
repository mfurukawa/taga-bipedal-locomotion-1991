// G. Taga, et al., 
// Self-organized control of bipedal locomotion by neural oscillators in unpredictable enviromnent, 
// Biological Cyberetics, 65, 147-150, 1991
//
// Masahiro Furukawa
// m.furukawa@ist.osaka-u.ac.jp
//
// rev 0.0, Created,  Oct 13, 2016-
//
// reference: 
// http://www.birl.ethz.ch/education/open_resource

#include <stdio.h>
#include <memory.h>
#include <stdlib.h>
#include "taga1991.h"

//#define __DUMP_MATRIX__TAGA1991__

Taga1991::Taga1991()
{
	// D. Simlumation Parameters

	// Musculo-skeletal system

	M = 48.0;
	m1 = 7.0;   l1 = 0.5;  I1 = m1 * l1 * l1 / 12.0;
	m2 = 4.0;   l2 = 0.6;  I2 = m2 * l2 * l2 / 12.0;
	b1 = 10.0;  b2 = 10.0;  bk = 1000.0;
	kk = 10000.0; g = 9.8;
	kg = 10000.0; bg = 1000.0;
	p_hf = 15.0;  p_he = 85.0;  
	p_kf = 15.0;  p_ke = 15.0;  
	p_af = 10.0;  p_ae = 124.0;

	// dt is time division in second
	dt = 0.00001;

	// neural rhythm generator

	memset(&tau [0],0x00,sizeof(tau)); // time constants of the inner state
	memset(&taud[0],0x00,sizeof(taud)); // the adaptation effect

	tau [1] = tau [2] = tau [3] = tau [4] = 0.05;
	taud[1] = taud[2] = taud[3] = taud[4] = 0.60;
	tau [5] = tau [6] = tau [7] = tau [8] = tau [9] = tau [10] = tau [11] = tau [12] = tau[1]/2.0;
	taud[5] = taud[6] = taud[7] = taud[8] = taud[9] = taud[10] = taud[11] = taud[12] = taud[1]/2.0;
	beta = 2.5;

	// D. Simlumation Parameters
	// neural rhythm generator

	memset(w,0x00,sizeof(w));  // a connecting weight

	w[1][2] = w[2][1] = w[3][4] = w[4][3] = w[5][6] = w[6][5] = w[7][8] = w[8][7] = w[9][10] = w[10][9] = w[11][12] = w[12][11] = w_fe = -2.0;
	w[1][3] = w[3][1] = w[2][4] = w[4][2] = w_rl = -1.0;
	w[6][1] = w[6][2] = w[8][3] = w[8][4] = w[10][1] = w[10][2] = w[12][3] = w[12][4] = w_hka = -1.0;

	// feedback

	memset(&a,0x00,sizeof(a));

	a[1] = 1.5;  a[2] = 1.0;  a[3] = 1.5;  a[4] = 1.5;
	a[5] = 3.0;  a[6] = 1.5;  a[7] = 3.0;  a[8] = 1.5;

	u[0] = 5.5; // Fig 5A

	init();
}
void Taga1991::init()
{
  memset(&u[1] ,0x00,sizeof(u)-1); // the inner state of the i-th neuron. u0 is an external input with a constant rate
  memset(&ud[0],0x00,sizeof(ud));
  memset(&v[0] ,0x00,sizeof(v)); // a variable represeinting the degree of the adaptation or self-inhibition effect of the i-th neuron 
  memset(&vd[0],0x00,sizeof(vd));


  // E. Initial condotion

  memset(&  x[0], 0x00,sizeof(x));
  memset(& xd[0], 0x00,sizeof(xd));
  memset(&xdd[0], 0x00,sizeof(xdd));
  memset(&  y[0], 0x00,sizeof(y)); // the output of the i-th neuron. (6)

	x1 = 0.0;
	x2 = 1.09;
	x5 = x11 = 0.45*M_PI;
	x8 = x14 = 0.57*M_PI;
	x3 = x1 + (l1 / 2.0)*c5;
	x4 = x2 - (l1 / 2.0)*s5;
	x6 = x1 + (l1 / 2.0)*c8;
	x7 = x2 - (l1 / 2.0)*s8;
	x9 = l1 * c5 + (l2 / 2.0)*c11;
	x10 = x2 - l1 * s5 - (l2 / 2.0)*s11;
	x12 = l1 * c8 + (l2 / 2.0)*c14;
	x13 = x2 - l1 * s8 - (l2 / 2.0)*s14;
	
  xr = 0.0; yr = 0.0; xl = 0.0; yl = 0.0;
  xr_d = 0.0; yr_d = 0.0; xl_d = 0.0; yl_d = 0.0;

	xr = xr0 = x9 + (l2 / 2.0)*c11;
	yr = yr0 = x10 - (l2 / 2.0)*s11;
	xl = xl0 = x12 + (l2 / 2.0)*c14;
	yl = yl0 = x13 - (l2 / 2.0)*s14;
}

void Taga1991::y_vec(void)
{
	// yi is the output of ith neuron 
	// see also equation (6) 

	for (int i = 1; i <= 12; i++)
		y[i] = f(u[i]);
}
void Taga1991::xrl_vec(void)
{
	// A. The equations of motion for the bipedal musculo-skeletal system

	// (xr,yr) and (xl, yl) represent the positions of the ankles, which are given by:

	xr = x9  + (l2 / 2.0)*c11;
	yr = x10 - (l2 / 2.0)*s11;
	xl = x12 + (l2 / 2.0)*c14;
	yl = x13 - (l2 / 2.0)*s14;

	xr_d = xd9  - (l2 / 2.0)*s11*xd11;
	yr_d = xd10 - (l2 / 2.0)*c11*xd11; // debugged
	xl_d = xd12 - (l2 / 2.0)*s14*xd14;
	yl_d = xd13 - (l2 / 2.0)*c14*xd14; // debugged
}
void Taga1991:: fg_vec(void)
{
	// yg(x) is function which represents the terrain. When the ground is level, yg(x) = 0.
	// Horizontal and vertical forces on the ankles are given by:

	Fg1 = 0.0; Fg2 = 0.0; Fg3 = 0.0; Fg4 = 0.0;				        // Horizontal and vertical forces on the ankles. See Fig. 12.

	if (yr - yg(xr) < 0.0) {
	    Fg1 = -kg*(xr - xr0) - bg*xr_d;
		Fg2 = -kg*(yr - yr0) + bg*f(-yr_d);
	} else {
		Fg1 = 0.0;
		Fg2 = 0.0;
	}

	if (yl - yg(xl) < 0.0) {
		Fg3 = -kg*(xl - xl0) - bg*xl_d;
		Fg4 = -kg*(yl - yl0) + bg*f(-yl_d);
	} else {
		Fg3 = 0.0;
		Fg4 = 0.0;
	}
}
void Taga1991:: Tr_vec(void)
{
	// Torques generated at each joint are given by:

	Tr1 = 0.0; Tr2 = 0.0; Tr3 = 0.0; Tr4 = 0.0; Tr5 = 0.0; Tr6 = 0.0;	// Torque

	Tr1 = p_he*y[2] - p_hf*y[1];
	Tr2 = p_he*y[4] - p_hf*y[3];
	Tr3 = p_ke*y[6] - p_kf*y[5];
	Tr4 = p_ke*y[8] - p_kf*y[7];
	Tr5 = (p_ae*y[10] - p_af*y[9]) *h(Fg2);
	Tr6 = (p_ae*y[12] - p_af*y[11])*h(Fg4);
}
void Taga1991:: Feed_vec(void)
{
	// Feedback pathway

	memset(&Feed[0],0x00,sizeof(Feed));

	Feed1 = a1 * (x5 - M_PI_2) - a2 * (x8 - M_PI_2) + a3 * (x11 - M_PI_2)*h(Fg2) + a4 * h(Fg4);
	Feed2 = a1 * (M_PI_2 - x5) - a2 * (M_PI_2 - x8) + a3 * (M_PI_2 - x11)*h(Fg2) + a4 * h(Fg4);
	Feed3 = a1 * (x8 - M_PI_2) - a2 * (x5 - M_PI_2) + a3 * (x14 - M_PI_2)*h(Fg4) + a4 * h(Fg2);
	Feed4 = a1 * (M_PI_2 - x8) - a2 * (M_PI_2 - x5) + a3 * (M_PI_2 - x14)*h(Fg4) + a4 * h(Fg2);
	Feed5 = a5 * (M_PI_2 - x14)*h(Fg4);
	Feed6 = a5 * (x14 - M_PI_2)*h(Fg4);
	Feed7 = a5 * (M_PI_2 - x11)*h(Fg2);
	Feed8 = a5 * (x11 - M_PI_2)*h(Fg2);
	Feed9 =  a6 * (M_PI_2 - x11)*h(Fg2) + a7 * (M_PI_2 - x14)*h(Fg4) - a8 * xd11 * h(Fg2);
	Feed10 = a6 * (x11 - M_PI_2)*h(Fg2) + a7 * (x14 - M_PI_2)*h(Fg4) + a8 * xd11 * h(Fg2);
	Feed11 = a6 * (M_PI_2 - x14)*h(Fg4) + a7 * (M_PI_2 - x11)*h(Fg2) - a8 * xd14 * h(Fg4);
	Feed12 = a6 * (x14 - M_PI_2)*h(Fg4) + a7 * (x11 - M_PI_2)*h(Fg2) + a8 * xd14 * h(Fg4);
}
void Taga1991:: uv(void)
{
  double sum = 0.0;

  for (int i = 1; i <= 12; i++) {
	sum = 0.0;
	for (int j = 1; j <= 12; j++)
	  sum += w[i][j] * y[j];

	ud[i] = (-u[i] + sum - beta*v[i] + u[0] + Feed[i]) / tau[i];
	vd[i] = (-v[i] + y[i]) / taud[i];
  }
}
void Taga1991:: P_mat(void)
{
	// P[14][8]
	// The equations of motion of the bipedal musculo-skeletal system are derivered 
	// using the Newton-Euler method. All variables and conventions correspond to 
	// those shown in Fig.2 and Fig. 12.

	memset(&P,0x00,sizeof(P)); 

	P[1][1] = 1.0 / M;		P[1][3] = 1.0 / M;
	P[2][2] = 1.0 / M;		P[2][4] = 1.0 / M;
	P[3][1] = -1.0 / m1;	P[3][5] = 1.0 / m1;
	P[4][2] = -1.0 / m1;	P[4][6] = 1.0 / m1;

	P[6][3] = -1.0 / m1;	P[6][7] = 1.0 / m1;
	P[7][4] = -1.0 / m1;	P[7][8] = 1.0 / m1;

	P[9][5] = -1.0 / m2;
	P[10][6] = -1.0 / m2;

	P[12][7] = -1.0 / m2;
	P[13][8] = -1.0 / m2;

	P[5][1] = -l1*s5 / 2.0 / I1;
	P[5][2] = -l1*c5 / 2.0 / I1;
	P[5][5] = -l1*s5 / 2.0 / I1;
	P[5][6] = -l1*c5 / 2.0 / I1;

	P[8][3] = -l1*s8 / 2.0 / I1;
	P[8][4] = -l1*c8 / 2.0 / I1;
	P[8][7] = -l1*s8 / 2.0 / I1;
	P[8][8] = -l1*c8 / 2.0 / I1;

	P[11][5] = -l2*s11 / 2.0 / I2;
	P[11][6] = -l2*c11 / 2.0 / I2;

	P[14][7] = -l2*s14 / 2.0 / I2;
	P[14][8] = -l2*c14 / 2.0 / I2;
}
void Taga1991:: Q_mat(void)
{
	// Q[14]  feedback etc

	memset(&Q[0],   0x00,sizeof(Q));

	Q[2] = -g;
	Q[4] = -g;
	Q[5] = (Tr1 + Tr3 - b1*fabs(x5 - M_PI_2)*xd5 - (b2 + bk* f(x5 - x11))*(xd5 - xd11) - kk*h(x5 - x11)) / I1;
	Q[7] = -g; // debug
	Q[8] = (Tr2 + Tr4 - b1*fabs(x8 - M_PI_2)*xd8 - (b2 + bk* f(x8 - x14))*(xd8 - xd14) - kk*h(x8 - x14)) / I1;
	Q[9] =  Fg1 / m2;
	Q[10] = Fg2 / m2 - g;
	Q[11] = (-Tr3 - Tr5 -(b2 + bk*f(x5 - x11))*(xd11 - xd5) + kk*h(x5 - x11)) / I2;
	Q[12] = Fg3 / m2;
	Q[13] = Fg4 / m2 - g;
	Q[14] = (-Tr4 - Tr6 -(b2 + bk*f(x8 - x14))*(xd14 - xd8) + kk*h(x8 - x14)) / I2;
}
void Taga1991:: C_mat(void)
{
	// C[8][14]  constraint

	memset(&C,0x00,sizeof(C)); 

	C[1][1] = C[2][2] = C[3][1] = C[4][2] = C[5][3] = C[6][4]  = C[7][6]  = C[8][7]  =  1.0;
	C[1][3] = C[2][4] = C[3][6] = C[4][7] = C[5][9] = C[6][10] = C[7][12] = C[8][13] = -1.0;
	C[1][5] = -l1*s5 / 2.0;
	C[2][5] = -l1*c5 / 2.0;
	C[3][8] = -l1*s8 / 2.0;
	C[4][8] = -l1*c8 / 2.0;
	C[5][5] = -l1*s5 / 2.0;		C[5][11] = -l2*s11 / 2.0;
	C[6][5] = -l1*c5 / 2.0;		C[6][11] = -l2*c11 / 2.0;
	C[7][8] = -l1*s8 / 2.0;		C[7][14] = -l2*s14 / 2.0;
	C[8][8] = -l1*c8 / 2.0;		C[8][14] = -l2*c14 / 2.0;

}
void Taga1991:: D_mat(void)
{
	// D[8][1]

	memset(&D[0],         0x00,sizeof(D));     // constraint

	D[1] =   l1*c5*xd52 / 2.0;
	D[2] =  -l1*s5*xd52 / 2.0;
	D[3] =   l1*c8*xd82 / 2.0;
	D[4] =  -l1*s8*xd82 / 2.0;
	D[5] =  (l1*c5*xd52 + l2*c11*xd112) / 2.0;
	D[6] = -(l1*s5*xd52 + l2*s11*xd112) / 2.0;
	D[7] =  (l1*c8*xd82 + l2*c14*xd142) / 2.0;
	D[8] = -(l1*s8*xd82 + l2*s14*xd142) / 2.0;
}
int Taga1991::inv_CP_mat(void)
{
	// neuton-eular method - differential equations

	// CP[8][8] = C[8][14] * P[14][8] | product C(x)P(x) 
	
	memset(&CP,     0x00,sizeof(CP));
	memset(&inv_CP, 0x00,sizeof(inv_CP));

	for (int k = 1; k <= 8; k++) { 				// row idx for CP
		for (int j = 1; j <= 8; j++) {			// col idx for CP
			CP[k][j] = 0.0;
			for (int i = 1; i <= 14; i++)
				CP[k][j] += C[k][i] * P[i][j];

			inv_CP[k][j] = CP[k][j];			// prepare to calculate inverce matrix with gauss-jordan method

			b[j] = 1.0; 					// dummy (no use)
		}
	}

	// inv_CP[8][8] = CP[8][8]^-1 | calculate inverce matrix with gauss-jordan method

	if (!gauss_jordan(8, inv_CP, b)) return 0;

	return 1;
}
void Taga1991::xdd_vec(void)
{
	// CQ[8][1] = C[8][14] * Q[14][1] | product C(x)Q(x,xd,Tr(y),Fg(x,xd)) 
	
	memset(&CQ,        0x00,sizeof(CP)); 

	for (int j = 1; j <= 8; j++) {
		CQ[j] = 0.0;
		for (int i = 1; i <= 14; i++)
		  CQ[j] += C[j][i] * Q[i];
	}

	// DCQ[8][1] = D[8][1] - CQ[8][1] | subtruct {D(x,xd) - C(x)Q(x,xd,Tr(y),Fg(x,xd))}

	memset(&DCQ[0],       0x00,sizeof(DCQ)); 

	for (int i = 1; i <= 8; i++)
		DCQ[i] = D[i] - CQ[i];

	// (inv_CP)D[8][1] = inv_CP[8][8] , DCQ[8][1] 

	double inv_CP_D[9] = {};

	for (int k = 1; k <= 8; k++) {
	  inv_CP_D[k] = 0.0;
	  for (int i = 1; i <= 8; i++) 
		inv_CP_D[k] += inv_CP[k][i] * DCQ[i];
	}

	// XDD[14][1] = P[14][8] * (inv_CP)D[8][1] + Q[14][1] | product P(x){C(x)P(x)}^-1 {D(x,xd) - C(x)Q(x,xd,Tr(y),Fg(x,xd))} + Q(x,xd,Tr(y),Fg(x,xd))
	for (int j = 1; j <= 14; j++) {
		xdd[j] = 0.0;
		for (int i = 1; i <= 8; i++)
			xdd[j] += P[j][i] * inv_CP_D[i];
		xdd[j] += Q[j];
	}
}
int Taga1991::update(void)
{
  y_vec();
  xrl_vec();
  P_mat();
  C_mat();
  D_mat();
  fg_vec();
  Tr_vec();
  Q_mat();
  Feed_vec();
  if(!inv_CP_mat()) return 0;
  xdd_vec();
  uv();
  return 1;
}

int Taga1991::next(void)
{
  // Runge Kutta Method coefficient
  double k1[15][4];
  double k2[15][4];
  double k3[15][4];
  double k4[15][4];
  // escape current state 
  double u_esc[15], v_esc[15];
  double x_esc[15], xd_esc[15];

	// touch ground
	if (yr - yg(xr) > 0.0) { xr0 = xr; yr0 = yr; }
	if (yl - yg(xl) > 0.0) { xl0 = xl; yl0 = yl; }
  
	if(!update()) return 0;	// calcurate XDD, ud, vd	

	// Runge Kutta Nystrom Method (4th order) p.94

	// calc first coefficient k1 |  u[13][14] are dummy
	for (int i = 1; i <= 14; i++) {
		k1[i][0] = ud[i];
		k1[i][1] = vd[i];
		k1[i][2] = xd[i];
		k1[i][3] = xdd[i];
		// store current state
		 u_esc[i] =  u[i];
		 v_esc[i] =  v[i];
		 x_esc[i] =  x[i];
		xd_esc[i] = xd[i]; 
		// set next estimation state
		 u[i] += k1[i][0] / 2.0 * dt;
		 v[i] += k1[i][1] / 2.0 * dt;
		 x[i] += k1[i][2] / 2.0 * dt;
		xd[i] += k1[i][3] / 2.0 * dt;
	}

	if (!update()) return 0;	// re-calcurate XDD, ud, vd

	for (int i = 1; i <= 14; i++) {
		k2[i][0] = ud[i];
		k2[i][1] = vd[i];
		k2[i][2] = xd[i];
		k2[i][3] =xdd[i];

		 u[i] =  u_esc[i] + k2[i][0] / 2.0 * dt;
		 v[i] =  v_esc[i] + k2[i][1] / 2.0 * dt;
		 x[i] =  x_esc[i] + k2[i][2] / 2.0 * dt;
		xd[i] = xd_esc[i] + k2[i][3] / 2.0 * dt;
	}

	if (!update()) return 0;	// re-calcurate XDD, ud, vd

	for (int i = 1; i <= 14; i++) {
		k3[i][0] = ud[i];
		k3[i][1] = vd[i];
		k3[i][2] = xd[i];
		k3[i][3] =xdd[i];

 		 u[i] =  u_esc[i] + k3[i][0] * dt;
		 v[i] =  v_esc[i] + k3[i][1] * dt;
 		 x[i] =  x_esc[i] + k3[i][2] * dt;
		xd[i] = xd_esc[i] + k3[i][3] * dt;
	}

	if (!update()) return 0;	// re-calcurate XDD, ud, vd

	for (int i = 1; i <= 14; i++) {
		k4[i][0] = ud[i];
		k4[i][1] = vd[i];
		k4[i][2] = xd[i];
		k4[i][3] =xdd[i];

		// update state
		u[i] =  u_esc[i] + ((k1[i][0] + 2.0*k2[i][0] + 2.0*k3[i][0] + k4[i][0]) / 6.0 ) * dt; // u[13][14] are dummy
		v[i] =  v_esc[i] + ((k1[i][1] + 2.0*k2[i][1] + 2.0*k3[i][1] + k4[i][1]) / 6.0 ) * dt; // u[13][14] are dummy
		x[i] =  x_esc[i] + ((k1[i][2] + 2.0*k2[i][2] + 2.0*k3[i][2] + k4[i][2]) / 6.0 ) * dt; 
	   xd[i] = xd_esc[i] + ((k1[i][3] + 2.0*k2[i][3] + 2.0*k3[i][3] + k4[i][3]) / 6.0 ) * dt; 
	}

	return 1;
}
int Taga1991::dump(void) 
{
  printf("\n[DEBUG] int Taga1991::dump(void) \n");

  printf("\n\tx[14]\t\txd[14]\t\txdd[14]\n\n");	
  for (int i = 1; i <= 14; i++)
 	printf("\t % 1.1e\t % 1.1e\t % 1.1e\n", x[i], xd[i], xdd[i]);

  printf("\n\tu[12]\t\tud[12]\t\tv[12]\t\tvd[12]\n\n");	
  for (int i = 1; i <= 12; i++)
 	printf("\t% 1.1e\t% 1.1e\t % 1.1e\t% 1.1e\n", u[i], ud[i], v[i], vd[i]);
  return 1;
}

Taga1991::~Taga1991()
{
}

double Taga1991::f(double x)  
{
  return fmax(0,x);
} 

double Taga1991::h(double x) 
{
  if (x > 0.0) return 1.0; 
  else return 0.0; 
}

double Taga1991::yg(double x) 
{
  return 0.0;
}


