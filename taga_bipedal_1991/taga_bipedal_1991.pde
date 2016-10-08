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
float yg(float x){ return 0; }

// A. The equations of motion for the bipedal musculo-skeletal system
  
private float[] F   = new float[15]; // Force given on each link
private float[] Fg  = new float[7];  // Horizontal and vertical forces on the ankles. See Fig. 12. 
private float[] Tr  = new float[7];  // Torque

public float[] x        = new float[15];
public float[] x_dot    = new float[15];
public float[] x_dotdot = new float[15];
public float[] y        = new float[15]; // the output of the i-th neuron. (6)

private float xr, yr, xl, yl, xr0, yr0, xl0, yl0 = 0;
private float xr_dot, yr_dot, xl_dot, yl_dot = 0; 

public float[] u        = new float[13]; // the inner state of the i-th neuron. u0 is an external input with a constant rate
public float[] u_dot    = new float[13];
public float[] v        = new float[13]; // a variable represeinting the degree of the adaptation or self-inhibition effect of the i-th neuron 
public float[] v_dot    = new float[13];
public float[] tau      = new float[13];  // time constants of the inner state
public float[] tau_dash = new float[13];  // the adaptation effect
  

// C. Feedback pathway

// Feedback signals from the musculo-skeletal system to the neural rhythm generetor are given by:

private float[] Feed = new float[13];

// D. Simlumation Parameters

// Musculo-skeletal system

public float[] m        = new float[3];
public float[] l        = new float[3];
public float[] I        = new float[3];
public float[] b        = new float[3];
public float M, bk, kk, g, kg, bg, p_hf, p_he, p_kf, p_ke, p_af, p_ae = 75.0;

// neural rhythm generator

public float beta;
public float[][] w = new float[13][13];  // a connecting weight
public float w_fe, w_rl, w_hka;

// feedback

private float[] a = new float[9];


  
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
  
  tau     [1] = tau     [2] = tau     [3] = tau     [4] = 0.05;
  tau_dash[1] = tau_dash[2] = tau_dash[3] = tau_dash[4] = 0.60;
  tau     [5] = tau     [6] = tau     [7] = tau     [8] = tau     [9] = tau     [10] = tau     [11] = tau     [12] = 0.025;
  tau_dash[5] = tau_dash[6] = tau_dash[7] = tau_dash[8] = tau_dash[9] = tau_dash[10] = tau_dash[11] = tau_dash[12] = 0.30;
  beta = 2.5;
  
  for(int i=1; i<=12;i++) 
    for(int j=1; j<=12;j++)
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
  for(int i=1; i<=14; i++) x_dot[i] = 0.0;  
  for(int i=1; i<=12; i++) { u_dot[i] = v_dot[i] = u[i] = v[i] = 0.0; }   
  
  
  
  size(1024,700);
  background(0);
  stroke(255);
  line(0,500,1024,500);

}

void draw() {
  
  // A. The equations of motion for the bipedal musculo-skeletal system
  
  // Torques generated at each joint are given by:
  
  Tr[1] = p_he*y[2] - p_hf*y[1];
  Tr[2] = p_he*y[4] - p_hf*y[3];
  Tr[3] = p_ke*y[6] - p_kf*y[5];
  Tr[4] = p_ke*y[8] - p_kf*y[7];
  Tr[5] =(p_ae*y[10]- p_af*y[9]) *h(Fg[2]);
  Tr[6] =(p_ae*y[12]- p_af*y[11])*h(Fg[4]);
  
  // The equations of motion of the bipedal musculo-skeletal system are derivered 
  // using the Newton-Euler method. All variables and conventions correspond to 
  // those shown in Fig.2 and Fig. 12.
  
  x_dotdot[1] = ( F[1]+F[3])/M; 
  x_dotdot[2] = ( F[2]+F[4])/M-g; 
  x_dotdot[3] = (-F[1]+F[5])/m[1];
  x_dotdot[4] = (-F[2]+F[6])/m[1]-g;
  x_dotdot[5] = (-F[1]*(l[1]/2.0)*sin(x[5]) - F[2]*(l[1]/2.0)*cos(x[5])-F[5]*(l[1]/2.0)*sin(x[5]) - F[6]*(l[1]/2.0)*cos(x[5]) - b[1]*abs(x[5]-PI/2.0)*x_dot[5] - (b[2]+bk*f(x[5]-x[11]))*(x[5]-x[11])-kk*h(x[5]-x[11])+Tr[1]+Tr[3]) / I[1];
  x_dotdot[6] = (-F[3]+F[7])/m[1];
  x_dotdot[7] = (-F[4]+F[5]-m[1]*g)/m[1];
  x_dotdot[8] = (-F[3]*(l[1]/2.0)*sin(x[8]) - F[4]*(l[1]/2.0)*cos(x[8])-F[7]*(l[1]/2.0)*sin(x[8]) - F[8]*(l[1]/2.0)*cos(x[8]) - b[1]*abs(x[8]-PI/2.0)*x_dot[8] - (b[2]+bk*f(x[8]-x[14]))*(x[8]-x[14])-kk*h(x[8]-x[14])+Tr[2]+Tr[4]) / I[1];
  x_dotdot[9] = (-F[5]+Fg[1])/m[2];
  x_dotdot[10]= (-F[6]+Fg[2])/m[2]-g;
  x_dotdot[11]= (-F[5]*(l[2]/2.0)*sin(x[11]) - F[6]*(l[2]/2.0)*cos(x[11])-Fg[1]*(l[2]/2.0)*sin(x[11]) - Fg[2]*(l[2]/2.0)*cos(x[11]) - (b[2]+bk*f(x[5]-x[11]))*(x_dot[11]-x_dot[5])+kk*h(x[5]-x[11])-Tr[3]-Tr[5]) / I[2];
  x_dotdot[12]= (-F[7]+Fg[3])/m[2];
  x_dotdot[13]= (-F[8]+Fg[4])/m[2]-g; 
  x_dotdot[14]= (-F[7]*(l[2]/2.0)*sin(x[14]) - F[8]*(l[2]/2.0)*cos(x[14])-Fg[3]*(l[2]/2.0)*sin(x[14]) - Fg[4]*(l[2]/2.0)*cos(x[14]) - (b[2]+bk*f(x[8]-x[14]))*(x_dot[14]-x_dot[8])+kk*h(x[8]-x[14])-Tr[4]-Tr[6]) / I[2];
  
  
  // Feedback pathway
  
  Feed[1] = a[1]*(x[5]-PI/2.0) - a[2]*(x[8]-PI/2.0) + a[3]*(x[11]-PI/2.0)*h(Fg[2]) + a[4]*h(Fg[4]);  
  Feed[2] = a[1]*(PI/2.0-x[5]) - a[2]*(PI/2.0-x[8]) + a[3]*(PI/2.0-x[11])*h(Fg[2]) + a[4]*h(Fg[4]); 
  Feed[3] = a[1]*(x[8]-PI/2.0) - a[2]*(x[5]-PI/2.0) + a[3]*(x[14]-PI/2.0)*h(Fg[4]) + a[4]*h(Fg[2]);  
  Feed[4] = a[1]*(PI/2.0-x[8]) - a[2]*(PI/2.0-x[5]) + a[3]*(PI/2.0-x[14])*h(Fg[4]) + a[4]*h(Fg[2]);  
  Feed[5] = a[5]*(PI/2.0-x[14])*h(Fg[4]);
  Feed[6] = a[5]*(x[14]-PI/2.0)*h(Fg[4]);
  Feed[7] = a[5]*(PI/2.0-x[11])*h(Fg[2]);
  Feed[8] = a[5]*(x[11]-PI/2.0)*h(Fg[2]);
  Feed[9] = a[6]*(PI/2.0-x[11])*h(Fg[2]) - a[7]*(PI/2.0-x[14])*h(Fg[4]) + a[8]*x_dot[11]*h(Fg[2]); 
  Feed[10]= a[6]*(x[11]-PI/2.0)*h(Fg[2]) - a[7]*(x[14]-PI/2.0)*h(Fg[4]) + a[8]*x_dot[11]*h(Fg[2]);  
  Feed[11]= a[6]*(PI/2.0-x[14])*h(Fg[2]) - a[7]*(PI/2.0-x[11])*h(Fg[4]) + a[8]*x_dot[14]*h(Fg[4]);  
  Feed[12]= a[6]*(x[14]-PI/2.0)*h(Fg[2]) - a[7]*(x[11]-PI/2.0)*h(Fg[4]) + a[8]*x_dot[14]*h(Fg[4]);  
  
  // (6) derivered from Matsuoka 1985,1987 and Grillner 1981.
  
  for(int i=1; i<=12; i++) {
    y[i] = f(u[i]);
   
    print(u[i]+"\t");
   
    u_dot[i] = -u[i] - beta*v[i] + u[0] + Feed[i];
    v_dot[i] = -v[i] + y[i];
    
    for(int j=1; j<=12; j++) 
      u_dot[i] = u_dot[i] + w[i][j]*y[j];
    
    u_dot[i] = u_dot[i] / tau[i];
    v_dot[i] = u_dot[i] / tau_dash[i];
  }
  
  // (xr,yr) and (xl, yl) represent the positions of the ankles, which are given by:
  
  xr = x[9] +(l[2]/2.0)*cos(x[11]);
  yr = x[10]+(l[2]/2.0)*sin(x[11]);
  xl = x[12]+(l[2]/2.0)*cos(x[14]);
  yl = x[13]+(l[2]/2.0)*sin(x[14]);
  
  
  // yg(x) is function which represents the terrain. When the ground is level, yg(x) = 0.
  // Horizontal and vertical forces on the ankles are given by:

  if (yr-yg(xr) < 0)  {
    Fg[1] = -kg*(xr-xr0) - bg*xr_dot;
    Fg[2] = -kg*(yr-yr0) - bg*f(-yr_dot);
  }else{            
    Fg[1] = 0;
    Fg[2] = 0;
  }
  
  if (yl-yg(xl) < 0)  {
    Fg[3] = -kg*(xl-xl0) - bg*xl_dot;
    Fg[4] = -kg*(yl-yl0) - bg*f(-yl_dot);
  }else{            
    Fg[3] = 0;
    Fg[4] = 0;
  }
  
  // Numerical integration
  
  for(int i=1; i<=14; i++)  {
    x_dot[i] += x_dotdot[i]; 
    x[i]     += x_dot[i];
    
    /* print(x[i]+"\t"); */
  }

  for(int i=1; i<=12; i++)  {
    u[i] += u_dot[i];
    v[i] += v_dot[i];
  } 
  
  line(x[12],x[13]+500,xl,yl+500);
  
  println();
}