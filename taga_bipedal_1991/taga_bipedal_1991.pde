/**
 * Masahiro Furukawa
 * m.furukawa@ist.osaka-u.ac.jp
 * 
 * Oct 17, 2016 -
 *
 * taga 1991
 */

double speed_scale = 1.0;

int fps = 120;

String[] lines;
int index = 0;
int stp = 0;
boolean flag = false;
int[] idx = new int[1000];
int idx_cnt = 0;

void setup() {
  size(1400, 200);
  
  background(255);
  stroke(0);
  strokeWeight(2);
  frameRate(fps);         
  lines = loadStrings("/Users/Masahiro/taga-bipedal-locomotion-1991/taga_bipedal_1991_c/taga_bipedal_1991_c/taga_bipedal_1991_c/out.csv");

  String[] pieces;
  pieces = split(lines[0], ','); double t0 = float(pieces[0]); //println(t0);
  pieces = split(lines[1], ','); double t1 = float(pieces[0]); //println(t1);
  stp = (int)(1.0/(double)fps * speed_scale / (t1 - t0));

}
 
 
void draw() {

  //if(index == stp) delay(2000);
  
  background(255);
  
  stroke(0);
  strokeWeight(1);
  line(0,150,1400,150);
  
  strokeWeight(8);
  
  if (index < lines.length) {
    //println (lines[index]);
    String[] pieces = split(lines[index], ',');
    if (pieces.length >= 2) {
      //ellipse(float(pieces[1])*100+50, -float(pieces[2])*100+150,10,10);
      
  stroke(255,0,0);
      line(float(pieces[3])*100+50 + cos(float(pieces[5])   )*25, -float(pieces[4])*100+150 + sin(float(pieces[5])   )*25, 
           float(pieces[3])*100+50 + cos(float(pieces[5])+PI)*25, -float(pieces[4])*100+150 + sin(float(pieces[5])+PI)*25);
      
  stroke(0,255,0);
      line(float(pieces[6])*100+50 + cos(float(pieces[8])   )*25, -float(pieces[7])*100+150 + sin(float(pieces[8])   )*25, 
           float(pieces[6])*100+50 + cos(float(pieces[8])+PI)*25, -float(pieces[7])*100+150 + sin(float(pieces[8])+PI)*25);
      
  strokeWeight(8);
  
  stroke(0,0,255);
      line(float(pieces[9])*100+50 + cos(float(pieces[11])   )*30, -float(pieces[10])*100+150 + sin(float(pieces[11])   )*30, 
           float(pieces[9])*100+50 + cos(float(pieces[11])+PI)*30, -float(pieces[10])*100+150 + sin(float(pieces[11])+PI)*30);
      
  stroke(255,0,255);
      line(float(pieces[12])*100+50 + cos(float(pieces[14])   )*30, -float(pieces[13])*100+150 + sin(float(pieces[14])   )*30, 
           float(pieces[12])*100+50 + cos(float(pieces[14])+PI)*30, -float(pieces[13])*100+150 + sin(float(pieces[14])+PI)*30);
    
    }
    
    if (index % (stp*15) == 0) idx[idx_cnt++] = index;    
    
    // Go to the next line for the next run through draw()
    index = index + stp;

    // repeat
    if(float(pieces[3])*100+50 + cos(float(pieces[5])   )*25 > 1400 || index > lines.length){
       index = 0;
       idx_cnt = 0;
    }
    
  }

 // trail();
}

void trail() {
  for( int i = 0; i < idx_cnt ; i++) {
        
    if (idx[i] < lines.length) {
      //println (lines[idx[i]]);
      String[] pieces = split(lines[idx[i]], ',');
      
      if (pieces.length >= 2) 
      {
          strokeWeight(2);
          stroke(0);
          
              line(float(pieces[3])*100+50 + cos(float(pieces[5])   )*25, -float(pieces[4])*100+150 + sin(float(pieces[5])   )*25, 
                   float(pieces[3])*100+50 + cos(float(pieces[5])+PI)*25, -float(pieces[4])*100+150 + sin(float(pieces[5])+PI)*25);
              
              line(float(pieces[6])*100+50 + cos(float(pieces[8])   )*25, -float(pieces[7])*100+150 + sin(float(pieces[8])   )*25, 
                   float(pieces[6])*100+50 + cos(float(pieces[8])+PI)*25, -float(pieces[7])*100+150 + sin(float(pieces[8])+PI)*25);
              
              line(float(pieces[9])*100+50 + cos(float(pieces[11])   )*30, -float(pieces[10])*100+150 + sin(float(pieces[11])   )*30, 
                   float(pieces[9])*100+50 + cos(float(pieces[11])+PI)*30, -float(pieces[10])*100+150 + sin(float(pieces[11])+PI)*30);
              
              line(float(pieces[12])*100+50 + cos(float(pieces[14])   )*30, -float(pieces[13])*100+150 + sin(float(pieces[14])   )*30, 
                   float(pieces[12])*100+50 + cos(float(pieces[14])+PI)*30, -float(pieces[13])*100+150 + sin(float(pieces[14])+PI)*30);
            
      }
    }
  }
}

void keyPressed() {
 index = 0;
 idx_cnt = 0;
}
void mousePressed() {
 index = 0;
 idx_cnt = 0;
}