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
	p_af = 58.95;  p_ae = 171.7; // 200

	// dt is time division in second
	dt = 0.00005;

	// neural rhythm generator

	memset(&tau [0],0x00,sizeof(tau)); // time constants of the inner state
	memset(&taud[0],0x00,sizeof(taud)); // the adaptation effect

	tau [1] = tau [2] = tau [3] = tau [4] = 0.0895;
	taud[1] = taud[2] = taud[3] = taud[4] = 0.60;
	tau [5] = tau [6] = tau [7] = tau [8] = tau [9] = tau [10] = tau [11] = tau [12] = 0.0895/2.0;
	taud[5] = taud[6] = taud[7] = taud[8] = taud[9] = taud[10] = taud[11] = taud[12] = 0.30;
	beta = 2.5;

	// D. Simlumation Parameters
	// neural rhythm generator

	memset(&w[0][0],0x00,sizeof(w));  // a connecting weight

	w[1][2] = w[2][1] = w[3][4] = w[4][3] = w[5][6] = w[6][5] = w[7][8] = w[8][7] = w[9][10] = w[10][9] = w[11][12] = w[12][11] = w_fe = -2.0;
	w[1][3] = w[3][1] = w[2][4] = w[4][2] = w_rl = -1.0;
	w[6][1] = w[6][2] = w[8][3] = w[8][4] = w[10][1] = w[10][2] = w[12][3] = w[12][4] = w_hka = -1.0;

	// feedback

	memset(&a,0x00,sizeof(a));

	a[1] = 1.5;  a[2] = 1.0;  a[3] = 1.5;  a[4] = 1.5;
	a[5] = 3.0;  a[6] = 1.5;  a[7] = 3.0;  a[8] = 1.5;

	u[0] = 8.3; // Fig 5A

	init();
}
void Taga1991::init()
{
  // Newton-Eular method
  memset(&D[0],         0x00,sizeof(D));     // constraint
  memset(&CQ[0],        0x00,sizeof(CP)); 
  memset(&DCQ[0],       0x00,sizeof(DCQ)); 
  memset(&CP[0][0],     0x00,sizeof(CP));
  memset(&Pinv_CP[0][0],0x00,sizeof(Pinv_CP)); 
  memset(&inv_CP[0][0], 0x00,sizeof(inv_CP));
  
  memset(&u[1] ,0x00,sizeof(u)-1); // the inner state of the i-th neuron. u0 is an external input with a constant rate
  memset(&ud[0],0x00,sizeof(ud));
  memset(&v[0] ,0x00,sizeof(v)); // a variable represeinting the degree of the adaptation or self-inhibition effect of the i-th neuron 
  memset(&vd[0],0x00,sizeof(vd));

  // C. Feedback pathway

  // Feedback signals from the musculo-skeletal system to the neural rhythm generetor are given by:
  memset(&Feed[0],0x00,sizeof(Feed));

  // A. The equations of motion for the bipedal musculo-skeletal system
  Fg1 = 0.0; Fg2 = 0.0; Fg3 = 0.0; Fg4 = 0.0;				        // Horizontal and vertical forces on the ankles. See Fig. 12.
  Tr1 = 0.0; Tr2 = 0.0; Tr3 = 0.0; Tr4 = 0.0; Tr5 = 0.0; Tr6 = 0.0;	// Torque

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
	flag_r = 0;
	flag_l = 0;

	// init P[14][8]

	memset(&P[0][0],0x00,sizeof(P)); // position

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


	// init Q[14]

	memset(&Q[0],   0x00,sizeof(Q));    // feedback etc

	Q[2] = Q[4] = Q[7] = Q[10] = Q[13] = -g;

	// init C[8][14]

	memset(&C[0][0],0x00,sizeof(C)); // constraint

	C[1][1] = C[2][2] = C[3][1] = C[4][2] = C[5][3] = C[6][4]  = C[7][6]  = C[8][7]  =  1.0;
	C[1][3] = C[2][4] = C[3][6] = C[4][7] = C[5][9] = C[6][10] = C[7][12] = C[8][13] = -1.0;
}

int Taga1991::update(void)
{
#ifdef __DUMP_MATRIX__TAGA1991__
  printf("\n\n\n\n\n\n");
#endif

	// yi is the output of ith neuron 
	// see also equation (6) 

	for (int i = 1; i <= 12; i++)
		y[i] = f(u[i]);

	// A. The equations of motion for the bipedal musculo-skeletal system

	// (xr,yr) and (xl, yl) represent the positions of the ankles, which are given by:

	xr = x9  + (l2 / 2.0)*c11;
	yr = x10 - (l2 / 2.0)*s11;
	xl = x12 + (l2 / 2.0)*c14;
	yl = x13 - (l2 / 2.0)*s14;

	xr_d = xd9  - (l2 / 2.0)*s11*xd11;
	yr_d = xd10 + (l2 / 2.0)*c11*xd11;
	xl_d = xd12 - (l2 / 2.0)*s14*xd14;
	yl_d = xd13 + (l2 / 2.0)*c14*xd14;

	// yg(x) is function which represents the terrain. When the ground is level, yg(x) = 0.
	// Horizontal and vertical forces on the ankles are given by:

	if (yr - yg(xr) < 0.0) {
	  if(!flag_r) {
		flag_r = 1;
	  }
	    Fg1 = -kg*(xr - xr0) - bg*xr_d;
		Fg2 = -kg*(yr - yr0) + bg*f(-yr_d);
	}
	else {
		xr0 = xr;
		yr0 = yr;
		Fg1 = 0.0;
		Fg2 = 0.0;
		flag_r = 0;
	}

	if (yl - yg(xl) < 0.0) {
	  if(!flag_l) {
		flag_l = 1;
	  }
		Fg3 = -kg*(xl - xl0) - bg*xl_d;
		Fg4 = -kg*(yl - yl0) + bg*f(-yl_d);
	}
	else {
		xl0 = xl;	
		yl0 = yl; 
		Fg3 = 0.0;
		Fg4 = 0.0;
		flag_l = 0;
	}

	// Torques generated at each joint are given by:

	Tr1 = p_he*y[2] - p_hf*y[1];
	Tr2 = p_he*y[4] - p_hf*y[3];
	Tr3 = p_ke*y[6] - p_kf*y[5];
	Tr4 = p_ke*y[8] - p_kf*y[7];
	Tr5 = (p_ae*y[10] - p_af*y[9]) *h(Fg2);
	Tr6 = (p_ae*y[12] - p_af*y[11])*h(Fg4);

	// Feedback pathway

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
	
#ifdef __DUMP_MATRIX__TAGA1991__
	printf("\n [DEBUG] int Taga1991::update(void)  >  Feed[12]\n\n");
	for(int i=1; i<=12; i++) 
	  printf("\t% 4.2e", Feed[i]);
	printf("\n");
#endif

	// neural rhythm generator - differential equations

#ifdef __DUMP_MATRIX__TAGA1991__
	printf("\n [DEBUG] int Taga1991::update(void)  >  ud[12], vd[12]\n\n");
#endif

	for (int i = 1; i <= 12; i++) {
		ud[i] = -u[i];
		for (int j = 1; j <= 12; j++)
			ud[i] += w[i][j] * y[j];

		ud[i] = ud[i] - beta*v[i] + u[0] + Feed[i];
		ud[i] /= tau[i];

		vd[i] = (-v[i] + y[i]) / taud[i];
	}

#ifdef __DUMP_MATRIX__TAGA1991__
	printf("\tud[12]\t\tvd[12]\n\n");
	for(int i=1; i<=12; i++)   
	  printf("\t% 4.2e\t% 4.2e\n", ud[i],vd[i]);	
#endif


	// P[14][8]
	// The equations of motion of the bipedal musculo-skeletal system are derivered 
	// using the Newton-Euler method. All variables and conventions correspond to 
	// those shown in Fig.2 and Fig. 12.

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

	// Q[14]
	Q[5] = (Tr1 + Tr3 - b1*fabs(x5 - M_PI_2)*xd5 - (b2 + bk*f(x5 - x11))*(xd5 - xd11) - kk*h(x5 - x11)) / I1;
	Q[8] = (Tr2 + Tr4 - b1*fabs(x8 - M_PI_2)*xd8 - (b2 + bk*f(x8 - x14))*(xd8 - xd14) - kk*h(x8 - x14)) / I1;
	Q[9]  = Fg1 / m2;
	Q[10] = Fg2 / m2;
	Q[11] = (-Tr3 - Tr5 -(b2 + bk*f(x5 - x11))*(xd11 - xd5) + kk*h(x5 - x11)) / I2;
	Q[12] = Fg3 / m2;
	Q[13] = Fg4 / m2;
	Q[14] = (-Tr4 - Tr6 -(b2 + bk*f(x8 - x14))*(xd14 - xd8) + kk*h(x8 - x14)) / I2;

#ifdef __DUMP_MATRIX__TAGA1991__
	printf("\n\tQ[5] = (-b1*fabs(x5 - M_PI_2)*xd5 - (b2 + bk*f(x5 - x11))*(xd5 - xd11) - kk*h(x5 - x11) + Tr1 + Tr3) / I1;\n\n");
	printf("\tQ[5]      % 1.4e\n",Q[5]);
	printf("\tb1        % 1.4e\n",b1);
	printf("\tx5        % 1.4e\n",x5);
	printf("\t     x5 - M_PI_2         % 1.4e\n",x5 - M_PI_2);
	printf("\tfabs(x5 - M_PI_2)        % 1.4e\n",fabs(x5 - M_PI_2));
	printf("\txd5       % 1.4e\n",xd5);
	printf("\tb2        % 1.4e\n",b2);
	printf("\tbk        % 1.4e\n",bk);
	printf("\tf(x5-x11) % 1.4e\n",f(x5-x11));
	printf("\txd5-xd11  % 1.4e\n",(xd5-xd11));
	printf("\tx11       % 1.4e\n",x11);
	printf("\txd11      % 1.4e\n",xd11);
	printf("\th(x5-x11) % 1.4e\n",h(x5-x11));
	printf("\tkk        % 1.4e\n",kk);
	printf("\tTr1       % 1.4e\n",Tr1);
	printf("\tTr3       % 1.4e\n",Tr3);
	printf("\tI1        % 1.4e\n",I1);

	printf("\n [DEBUG] int Taga1991::update(void)  >  Q[14]\n\n");
	for (int k = 1; k <= 14; k++)		printf("\t% 4.2e\n", Q[k]);
#endif

	// C[8][14]

	C[1][5] = -l1*s5 / 2.0;
	C[2][5] = -l1*c5 / 2.0;
	C[3][8] = -l1*s8 / 2.0;
	C[4][8] = -l1*c8 / 2.0;
	C[5][5] = -l1*s5 / 2.0;		C[5][11] = -l2*s11 / 2.0;
	C[6][5] = -l1*c5 / 2.0;		C[6][11] = -l2*c11 / 2.0;
	C[7][8] = -l1*s8 / 2.0;		C[7][14] = -l2*s14 / 2.0;
	C[8][8] = -l1*c8 / 2.0;		C[8][14] = -l2*c14 / 2.0;

	// D[8][1]

	D[1] =   l1*c5*xd52 / 2.0;
	D[2] =  -l1*s5*xd52 / 2.0;
	D[3] =   l1*c8*xd82 / 2.0;
	D[4] =  -l1*s8*xd82 / 2.0;
	D[5] =  (l1*c5*xd52 + l2*c11*xd112) / 2.0;
	D[6] = -(l1*s5*xd52 + l2*s11*xd112) / 2.0;
	D[7] =  (l1*c8*xd82 + l2*c14*xd142) / 2.0;
	D[8] = -(l1*s8*xd82 + l2*s14*xd142) / 2.0;

	// neuton-eular method - differential equations

	// CP[8][8] = C[8][14] * P[14][8] | product C(x)P(x) 
	for (int k = 1; k <= 8; k++) { 				// row idx for CP
		for (int j = 1; j <= 8; j++) {			// col idx for CP
			CP[k][j] = 0.0;
			for (int i = 1; i <= 14; i++)
				CP[k][j] += C[k][i] * P[i][j];

			inv_CP[k][j] = CP[k][j];			// prepare to calculate inverce matrix with gauss-jordan method

			b[j] = 1.0; 					// dummy (no use)
		}
	}

#ifdef __DUMP_MATRIX__TAGA1991__
	printf("\n [DEBUG] int Taga1991::update(void) \n");
	printf("\n\tP[14][8]\n\n");	for (int k = 1; k <= 14; k++){	for (int j = 1; j <= 8; j++)	printf("\t% 1.0e", P[k][j]);		printf("\n");	}
	printf("\n\tC[8][14]\n\n");	for (int k = 1; k <= 8; k++){	for (int j = 1; j <= 14; j++)	printf("\t% 1.0e", C[k][j]);		printf("\n");	}
	printf("\n\tCP[8][8]\n\n");	for (int k = 1; k <= 8; k++){	for (int j = 1; j <= 8; j++)	printf("\t% 1.0e", CP[k][j]);	printf("\n");	}
	//	printf("\n\tinv_CP[8][8]\n\n");	for (int k = 1; k <= 8; k++){	for (int j = 1; j <= 8; j++)	printf("\t% 1.0e", inv_CP[k][j]);	printf("\n");	}
#endif

	// inv_CP[8][8] = CP[8][8]^-1 | calculate inverce matrix with gauss-jordan method
	if (!gauss_jordan(8, inv_CP, b)) return 0;

#ifdef __DUMP_MATRIX__TAGA1991__
	printf("\n [DEBUG] int Taga1991::update(void)  >  gauss_jordan(8, inv_CP, b)  >  inv_CP[i][j]\n\n");
	for (int k = 1; k <= 8; k++) {	 for (int j = 1; j <= 8; j++) 	printf("\t% 1.0e", inv_CP[k][j]); printf("\n"); }

	printf("\n [DEBUG] int Taga1991::update(void)  >  inv_CP[8][8] * CP[8][8]\n\n");
	for (int k = 1 ; k <= 8; k++){
	  for (int j = 1; j <= 8; j++){
		double tmp = 0.0;
		for (int i = 1; i <= 8; i++)
		  tmp += inv_CP[k][i] * CP[i][j];
		printf("\t% 1.3f", tmp);
	  }
	  printf("\n");
	}
#endif

	// CQ[8][1] = C[8][14] * Q[14][1] | product C(x)Q(x,xd,Tr(y),Fg(x,xd)) 
	for (int j = 1; j <= 8; j++) { 				// row idx for CP
		CQ[j] = 0.0;
		for (int i = 1; i <= 14; i++){
		  CQ[j] += C[j][i] * Q[i];

		}
	}
	// DCQ[8][1] = D[8][1] - CQ[8][1] | subtruct {D(x,xd) - C(x)Q(x,xd,Tr(y),Fg(x,xd))}
	for (int i = 1; i <= 8; i++)
		DCQ[i] = D[i] - CQ[i];
	// (inv_CP)D[8][1] = inv_CP[8][8] , DCQ[8][1] 
	double inv_CP_D[9];
	for (int k = 1; k <= 8; k++) {			// col idx for CP
	  inv_CP_D[k] = 0.0;
	  for (int i = 1; i <= 8; i++) {
		inv_CP_D[k] += inv_CP[k][i] * DCQ[i]; // NOTICE!! inv_CP's index number starts from '1'!
	  }
	}
	// XDD[14][1] = Pinv_CP[14][8] * DCQ[8][1] + Q[14][1] | product P(x){C(x)P(x)}^-1 {D(x,xd) - C(x)Q(x,xd,Tr(y),Fg(x,xd))} + Q(x,xd,Tr(y),Fg(x,xd))
	for (int j = 1; j <= 14; j++) { 				// row idx for CP
		xdd[j] = Q[j];
		for (int i = 1; i <= 8; i++){
			xdd[j] += P[j][i] * inv_CP_D[i];
		}
	}

	//exit(1);
	return 1;
}

int Taga1991::next(void)
{
  // Runge Kutta Method coefficient
  double k1[15][2]; // [4] means twice diferentiated x+ u + v
  double k2[15][2];
  double k3[15][2];
  double k4[15][2];
  double An[15], Bn[15], Cn[15], Dn[15], betan[15], deltan[15];
  // escape current state 
  double u_esc[15], v_esc[15];
  double x_esc[15], xd_esc[15];
  
	 
#ifdef __DUMP_MATRIX__TAGA1991__
	printf("\n[DEBUG] just entered next()");
	dump();
#endif

	if(!update()) return 0;	// calcurate XDD, ud, vd	

#ifdef __DUMP_MATRIX__TAGA1991__
	printf("\n[DEBUG] after update before k1");
	dump();
#endif

	// Runge Kutta Nystrom Method (4th order) p.94
	
#ifdef __DUMP_MATRIX__TAGA1991__
	printf("\n\n[DEBUG] int Taga1991::next(void)  >  k1[14][2] An betan\n\n\tdt*ud\t\tdt*vd\t\tAn\t\tbetan\n");
#endif

	for (int i = 1; i <= 14; i++) {
		// calc first coefficient k1 |  u[13][14] are dummy
		k1[i][0] = dt *  ud[i];
		k1[i][1] = dt *  vd[i];

		An[i] = xdd[i] * dt * 0.5; // An = kf(xn, yn, yn'), k=h/2, y'' = f(x, y, y') | h = dt
		betan[i] = (xd[i] + An[i] * 0.5) * 0.5 * dt;

#ifdef __DUMP_MATRIX__TAGA1991__
		printf("\t% 4.2e\t% 4.2e\t% 4.2e\t% 4.2e\n", k1[i][0], k1[i][1], An[i], betan[i]);
#endif
		// store current state
		 u_esc[i] =  u[i];
		 v_esc[i] =  v[i];

		 x_esc[i] =  x[i];
		xd_esc[i] = xd[i]; 
		// set next estimation state
		 u[i] += k1[i][0] / 2.0;
		 v[i] += k1[i][1] / 2.0;

		 x[i] =  x_esc[i] + betan[i]; // yn + betan : betan = k(yn' + An/2) | An = k1[i][3]
		xd[i] = xd_esc[i] + An[i];    // yn' + An : An = k1[i][3]
	}


	if (!update()) return 0;	// re-calcurate XDD, ud, vd

#ifdef __DUMP_MATRIX__TAGA1991__
	printf("\n[DEBUG] after k1");
	dump();

	printf("\n[DEBUG] int Taga1991::next(void)  >  k2[14][2] An Bn \n\n\tdt*ud\t\tdt*vd\t\tAn\t\tBn\n");
#endif

	for (int i = 1; i <= 14; i++) {
		// calc second coefficient k2 |  u[13][14] are dummy
		k2[i][0] = dt *  ud[i];
		k2[i][1] = dt *  vd[i];

		Bn[i] = xdd[i] * dt * 0.5; //  Bn = kf(xn + k, yn + betan, yn' + An)), xd[i] = yn' + An | Bn = k2[i][3]

#ifdef __DUMP_MATRIX__TAGA1991__
		printf("\t% 4.2e\t% 4.2e\t% 4.2e\t% 4.2e\n", k2[i][0], k2[i][1], An[i], Bn[i]);
#endif
		// set next estimation state
		 u[i] =  u_esc[i] + k2[i][0] / 2.0;
		 v[i] =  v_esc[i] + k2[i][1] / 2.0;

		 x[i] =  x_esc[i] + betan[i]; // yn + betan : betan = k(yn' + An/2) | An = k1[i][3], 
		xd[i] = xd_esc[i] + Bn[i];    // yn' + Bn : Bn = k2[i][3]
	}

	if (!update()) return 0;	// re-calcurate XDD, ud, vd

#ifdef __DUMP_MATRIX__TAGA1991__
	printf("\n[DEBUG] after k2");
	dump();

	printf("\n[DEBUG] int Taga1991::next(void)  >  k3[14][2] Cn deltan\n\n\tdt*ud\t\tdt*vd\t\tCn\t\tdeltan\n");
#endif

	for (int i = 1; i <= 14; i++) {
		// calc second coefficient k3 |  u[13][14] are dummy
		k3[i][0] = dt *  ud[i];
		k3[i][1] = dt *  vd[i];

		Cn[i] = xdd[i] * dt * 0.5; // Cn = kf(xn + k, yn + betan, yn' + Bn) | Cn = k3[i][3]
		deltan[i] = dt * (xd_esc[i] + Cn[i]);

#ifdef __DUMP_MATRIX__TAGA1991__
		printf("\t% 4.2e\t% 4.2e\t% 4.2e\t% 4.2e\n", k3[i][0], k3[i][1], Cn[i], deltan[i]);
#endif
		// set next estimation state
 		 u[i] =  u_esc[i] + k3[i][0];
		 v[i] =  v_esc[i] + k3[i][1];

		 x[i] =  x_esc[i] + deltan[i]; // yn + deltan : deltan = h(yn' + Cn) | Cn = k3[i][3]
		 xd[i] = xd_esc[i] + 2.0 * Cn[i]; // yn' + 2Cn : Cn = k3[i][3]
	}

	if (!update()) return 0;	// re-calcurate XDD, ud, vd

#ifdef __DUMP_MATRIX__TAGA1991__
	printf("\n[DEBUG] after k3");
	dump();

	printf("\n[DEBUG] int Taga1991::next(void)  >  k4[14][2] Cn Dn \n\n\tdt*ud\t\tdt*vd\t\tCn\t\tDn\n");
#endif

	for (int i = 1; i <= 14; i++) {
		// calc second coefficient k4 |  u[13][14] are dummy
		k4[i][0] = dt *  ud[i];
		k4[i][1] = dt *  vd[i];

		Dn[i] = xdd[i] * dt * 0.5; // Dn = kf(xn + h, yn + deltan, yn' + 2Cn), deltan = h(yn' + Cn) | Dn = k4[i][3]

#ifdef __DUMP_MATRIX__TAGA1991__
		printf("\t% 4.2e\t% 4.2e\t% 4.2e\t% 4.2e\n", k3[i][0], k3[i][1], Cn[i], Dn[i]);
#endif
		// update state
		  u[i] =  u_esc[i] + (k1[i][0] + 2.0*k2[i][0] + 2.0*k3[i][0] + k4[i][0]) / 6.0; // u[13][14] are dummy
 		  v[i] =  v_esc[i] + (k1[i][1] + 2.0*k2[i][1] + 2.0*k3[i][1] + k4[i][1]) / 6.0; // u[13][14] are dummy

		  x[i] =  x_esc[i] + dt * (xd_esc[i] + (An[i] + Bn[i] + Cn[i]) / 3.0 );         // yn+1 = yn + h[yn' + 1/3(An + Bn + Cn)]
		  xd[i]= xd_esc[i] + (An[i] + 2.0*Bn[i] + 2.0*Cn[i] + Dn[i]) / 3.0;
	}

#ifdef __DUMP_MATRIX__TAGA1991__
	printf("\n[DEBUG] after k4");
	dump();

	printf("\n\n\n\n\n\n");
#endif


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
  if (x<=0.0) return 0.0; // max(0,x);
  else return x; 
} 

double Taga1991::h(double x) 
{
  if (x<=0.0) return 0.0; 
  else return 1.0; 
}

double Taga1991::yg(double x) 
{
  return 0.0; 
}


