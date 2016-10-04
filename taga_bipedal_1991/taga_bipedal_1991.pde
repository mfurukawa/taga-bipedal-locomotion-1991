// G. Taga, et al., 
// Self-organized control of bipedal locomotion by neural oscillators in unpredictable enviromnent, 
// Biological Cyberetics, 65, 147-150, 1991
//
// Masahiro Furukawa
// m.furukawa@ist.osaka-u.ac.jp
//
// rev 0.0, Created,  Oct 4, 2016-


// C. Feedback pathway

double Feed1 (void) {
  return a[1]*(x[5] - PI/2.0) - a2*(x[8]-PI/2.0) + a[3]*(x[11] - PI/2.0)
} 


// D. Simlumation Parameters

// Musculo-skeletal system

double[] m        = new double[3];
double[] l        = new double[3];
double[] I        = new double[3];
double[] b        = new double[3];
double M;
M = 48.0;
m[1] = 7.0;   l[1] = 0.5;  I[1] = m[1] * l[1] * l[1] / 12.0;
m[2] = 4.0;   l[2] = 0.6;  I[2] = m[2] * l[2] * l[2] / 12.0;
b[1] = 10.0;  b[2] = 10.0;  bk = 1000.0;

double kk, kg, g, bg, p_hf, p_he, p_kf, p_ke, p_af, p_ae = 75.0;
kk = 10000.0; g = 9.8;
kg = 10000.0; bg = 1000.0;
p_hf = 15.0;  p_he = 85.0;  p_kf = 15.0;
p_ke = 15.0;  p_af = 100.0; p_ae = 75.0;

// neural rhythm generator

double[] tau      = new double[13];
double[] tau_dot  = new double[13];
tau    [1] = tau    [2] = tau    [3] = tau    [4] = 0.05;
tau_dot[1] = tau_dot[2] = tau_dot[3] = tau_dot[4] = 0.60;
tau    [5] = tau    [6] = tau    [7] = tau    [8] = tau    [9] = tau    [10] = tau    [11] = tau    [12] = 0.025;
tau_dot[5] = tau_dot[6] = tau_dot[7] = tau_dot[8] = tau_dot[9] = tau_dot[10] = tau_dot[11] = tau_dot[12] = 0.30;

double beta = 2.5;
double[][] w = new double[13][13];
double w_fe, w_rl, w_hka;
for(int i=0; i<13;i++) 
  for(int j=0; j<13;j++)
    w[i][j]=0.0;
w[1][2] = w[2][1] = w[3][4] = w[4][3] = w[5][6] = w[6][5] = w[7][8] = w[8][7] = w[9][10] = w[10][9] = w[11][12] = w[12][11] = w_fe = -2.0;
w[1][3] = w[3][1] = w[2][4] = w[4][2] = w_rl = -1.0;
w[6][5] = w[6][2] = w[8][3] = w[8][4] = w[10][1] = w[10][2] = w[12][3] = w[12][4] = w_hka = -1.0;

// feedback

double[] a = new double[9];
a[1] = 1.5;  a[2] = 1.0;  a[3] = 1.5;  a[4] = 1.5;
a[5] = 3.0;  a[6] = 1.5;  a[7] = 3.0;  a[8] = 1.5;


// E. Initial condotion

double[] x     = new double[13];
double[] x_dot = new double[15];
double[] u_dot = new double[13];
double[] v_dot = new double[13];
x[1] = 0.0;  x[2] = 1.09;  x[5] = x[11] = 0.45*PI;
x[8] = x[14] = 0.57*PI;
x[3] = x[1] + (l[1]/2.0)*cos(x[5]);  x[4] = x[2] - (l[1]/2.0)*sin(x[5]);
x[6] = x[1] + (l[1]/2.0)*cos(x[8]);  x[7] = x[2] - (l[1]/2.0)*sin(x[8]);
x[9]  =        l[1]*cos(x[5]) + (l[2]/2.0)*cos(x[11]);
x[10] = x[2] - l[1]*sin(x[5]) - (l[2]/2.0)*sin(x[11]);
x[12] =        l[1]*cos(x[8]) + (l[2]/2.0)*cos(x[14]);
x[13] = x[2] - l[1]*sin(x[8]) - (l[2]/2.0)*sin(x[14]);
for(int i=0; i<15; i++) x_dot[i] = 0.0;  
for(int i=0; i<13; i++) u_dot[i] = 0.0;  
for(int i=0; i<13; i++) v_dot[i] = 0.0;  

void setup() {
  
}

void draw() {
  ;
  
}