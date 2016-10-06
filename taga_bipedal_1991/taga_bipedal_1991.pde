// G. Taga, et al., 
// Self-organized control of bipedal locomotion by neural oscillators in unpredictable enviromnent, 
// Biological Cyberetics, 65, 147-150, 1991
//
// Masahiro Furukawa
// m.furukawa@ist.osaka-u.ac.jp
//
// rev 0.0, Created,  Oct 4, 2016-

float f(float x) { if (x<0) return 0; else return x; } // max(0,x);
float h(float x) { if (x<0) return 0; else return 1; }

  
// C. Feedback pathway

// Feedback signals from the musculo-skeletal system to the neural rhythm generetor are given by:

private float Fg2, Fg4;
private float[] Feed = new float[13];

// D. Simlumation Parameters

// Musculo-skeletal system

public float[] m        = new float[3];
public float[] l        = new float[3];
public float[] I        = new float[3];
public float[] b        = new float[3];
public float M, bk, kk, g, kg, bg, p_hf, p_he, p_kf, p_ke, p_af, p_ae = 75.0;

// neural rhythm generator

public float[] tau      = new float[13];
public float[] tau_dot  = new float[13];
public float beta;
public float[][] w = new float[13][13];
public float w_fe, w_rl, w_hka;

// feedback

private float[] a = new float[9];

// E. Initial condotion
  
public float[] x     = new float[15];
public float[] x_dot = new float[15];
public float[] u_dot = new float[13];
public float[] v_dot = new float[13];

void setup() {
  
  // D. Simlumation Parameters
  
  // Musculo-skeletal system

  M = 48.0;
  m[1] = 7.0;   l[1] = 0.5;  I[1] = m[1] * l[1] * l[1] / 12.0;
  m[2] = 4.0;   l[2] = 0.6;  I[2] = m[2] * l[2] * l[2] / 12.0;
  b[1] = 10.0;  b[2] = 10.0;  bk = 1000.0;
  kk = 10000.0; g = 9.8;
  kg = 10000.0; bg = 1000.0;
  p_hf = 15.0;  p_he = 85.0;  p_kf = 15.0;
  p_ke = 15.0;  p_af = 100.0; p_ae = 75.0;

  // neural rhythm generator
  
  tau    [1] = tau    [2] = tau    [3] = tau    [4] = 0.05;
  tau_dot[1] = tau_dot[2] = tau_dot[3] = tau_dot[4] = 0.60;
  tau    [5] = tau    [6] = tau    [7] = tau    [8] = tau    [9] = tau    [10] = tau    [11] = tau    [12] = 0.025;
  tau_dot[5] = tau_dot[6] = tau_dot[7] = tau_dot[8] = tau_dot[9] = tau_dot[10] = tau_dot[11] = tau_dot[12] = 0.30;
  beta = 2.5;
  
  for(int i=0; i<13;i++) 
  for(int j=0; j<13;j++)
      w[i][j]=0.0;
  w[1][2] = w[2][1] = w[3][4] = w[4][3] = w[5][6] = w[6][5] = w[7][8] = w[8][7] = w[9][10] = w[10][9] = w[11][12] = w[12][11] = w_fe = -2.0;
  w[1][3] = w[3][1] = w[2][4] = w[4][2] = w_rl = -1.0;
  w[6][5] = w[6][2] = w[8][3] = w[8][4] = w[10][1] = w[10][2] = w[12][3] = w[12][4] = w_hka = -1.0;
  
  // feedback
  
  a[1] = 1.5;  a[2] = 1.0;  a[3] = 1.5;  a[4] = 1.5;
  a[5] = 3.0;  a[6] = 1.5;  a[7] = 3.0;  a[8] = 1.5;
  
  // E. Initial condotion
  
  x[1] = 0.0;  
  x[2] = 1.09;  
  x[5] = x[11] = 0.45*PI;
  x[8] = x[14] = 0.57*PI;
  x[3] = x[1] + (l[1]/2.0)*cos(x[5]);
  x[4] = x[2] - (l[1]/2.0)*sin(x[5]);
  x[6] = x[1] + (l[1]/2.0)*cos(x[8]);
  x[7] = x[2] - (l[1]/2.0)*sin(x[8]);
  x[9]  =        l[1]*cos(x[5]) + (l[2]/2.0)*cos(x[11]);
  x[10] = x[2] - l[1]*sin(x[5]) - (l[2]/2.0)*sin(x[11]);
  x[12] =        l[1]*cos(x[8]) + (l[2]/2.0)*cos(x[14]);
  x[13] = x[2] - l[1]*sin(x[8]) - (l[2]/2.0)*sin(x[14]);
  for(int i=0; i<15; i++) x_dot[i] = 0.0;  
  for(int i=0; i<13; i++) u_dot[i] = 0.0;  
  for(int i=0; i<13; i++) v_dot[i] = 0.0;  
}

void draw() {
  
  Feed[1] = a[1]*(x[5]-PI/2.0) - a[2]*(x[8]-PI/2.0) + a[3]*(x[11]-PI/2.0)*h(Fg2) + a[4]*h(Fg4);  
  Feed[2] = a[1]*(PI/2.0-x[5]) - a[2]*(PI/2.0-x[8]) + a[3]*(PI/2.0-x[11])*h(Fg2) + a[4]*h(Fg4); 
  Feed[3] = a[1]*(x[8]-PI/2.0) - a[2]*(x[5]-PI/2.0) + a[3]*(x[14]-PI/2.0)*h(Fg4) + a[4]*h(Fg2);  
  Feed[4] = a[1]*(PI/2.0-x[8]) - a[2]*(PI/2.0-x[5]) + a[3]*(PI/2.0-x[14])*h(Fg4) + a[4]*h(Fg2);  
  Feed[5] = a[5]*(PI/2.0-x[14])*h(Fg4);
  Feed[6] = a[5]*(x[14]-PI/2.0)*h(Fg4);
  Feed[7] = a[5]*(PI/2.0-x[11])*h(Fg2);
  Feed[8] = a[5]*(x[11]-PI/2.0)*h(Fg2);
  Feed[9] = a[6]*(PI/2.0-x[11])*h(Fg2) - a[7]*(PI/2.0-x[14])*h(Fg4) + a[8]*x_dot[11]*h(Fg2); 
  Feed[10]= a[6]*(x[11]-PI/2.0)*h(Fg2) - a[7]*(x[14]-PI/2.0)*h(Fg4) + a[8]*x_dot[11]*h(Fg2);  
  Feed[11]= a[6]*(PI/2.0-x[14])*h(Fg2) - a[7]*(PI/2.0-x[11])*h(Fg4) + a[8]*x_dot[14]*h(Fg4);  
  Feed[12]= a[6]*(x[14]-PI/2.0)*h(Fg2) - a[7]*(x[11]-PI/2.0)*h(Fg4) + a[8]*x_dot[14]*h(Fg4);  
  

}