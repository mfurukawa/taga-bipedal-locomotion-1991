/**
 * LoadFile 1
 * 
 * Loads a text file that contains two numbers separated by a tab ('\t').
 * A new pair of numbers is loaded each frame and used to draw a point on the screen.
 */

String[] lines;
int index = 0;

void setup() {
  size(1000, 600);
  background(255);
  stroke(0);
  frameRate(120);
  lines = loadStrings("/Users/Masahiro/taga-bipedal-locomotion-1991/taga_bipedal_1991_c/taga_bipedal_1991_c/taga_bipedal_1991_c/out.csv");
}

void draw() {
  
  background(255);
  line(0,500,1200,500);
  
  if (index < lines.length) {
    //println (lines[index]);
    String[] pieces = split(lines[index], ',');
    if (pieces.length >= 2) {
      ellipse(float(pieces[1])*100+50, -float(pieces[2])*100+500,10,10);
      
      line(float(pieces[3])*100+50 + cos(float(pieces[5])   )*25, -float(pieces[4])*100+500 + sin(float(pieces[5])   )*25, 
           float(pieces[3])*100+50 + cos(float(pieces[5])+PI)*25, -float(pieces[4])*100+500 + sin(float(pieces[5])+PI)*25);
      
      line(float(pieces[6])*100+50 + cos(float(pieces[8])   )*25, -float(pieces[7])*100+500 + sin(float(pieces[8])   )*25, 
           float(pieces[6])*100+50 + cos(float(pieces[8])+PI)*25, -float(pieces[7])*100+500 + sin(float(pieces[8])+PI)*25);
      
      line(float(pieces[9])*100+50 + cos(float(pieces[11])   )*30, -float(pieces[10])*100+500 + sin(float(pieces[11])   )*30, 
           float(pieces[9])*100+50 + cos(float(pieces[11])+PI)*30, -float(pieces[10])*100+500 + sin(float(pieces[11])+PI)*30);
      
      line(float(pieces[12])*100+50 + cos(float(pieces[14])   )*30, -float(pieces[13])*100+500 + sin(float(pieces[14])   )*30, 
           float(pieces[12])*100+50 + cos(float(pieces[14])+PI)*30, -float(pieces[13])*100+500 + sin(float(pieces[14])+PI)*30);
    }
    // Go to the next line for the next run through draw()
    index = index + 50;    
  }
}