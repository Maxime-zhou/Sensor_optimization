
HBE=.1;   // height bottom electrode   
HALN=1.;      // height AlN
HTE=.1;    //  height top electrode   

Dinf=10;

WF=3.5;  // width finger
GF=2.75;  // gap finger
NF=2;   // number of fingers

lc0=HBE/1;   // 
lc1=HALN/1;   // 
lc2=4;

WT=(WF+GF)*NF+GF;  // largh totale
HT=HBE+HALN+HTE;  // alt totale

//
// bottom electrode
//

Point(1) = {0,0,0,lc0}; 
Point(2) = {WT,0,0,lc0}; 
Point(3) = {WT,HBE,0,lc0}; 
Point(4) = {0,HBE,0,lc0}; 

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

//
// boundary at infinity
//

Point(5) = {-Dinf,-Dinf,0,lc2}; 
Point(6) = {WT+Dinf,-Dinf,0,lc2}; 
Point(7) = {WT+Dinf,HT+Dinf,0,lc2}; 
Point(8) = {-Dinf,HT+Dinf,0,lc2}; 

Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,5};
Line Loop(9) = {5,6,7,8};    // infinity boundary

//
// upper electrode
//

y1=HBE+HALN;
y2=y1+HTE;
x=0;

sA[]={};  // lista boundary Air
sR[]={};  // lista boundary Resonator

pnumI=9;   // numero inziale per i punti
pnum=pnumI;  
lnum=10;    // numero inziale per le linee
snum=1;
 
For t In {1:NF+1}    // crea i gaps

 Point(pnum) = {x,y1,0,lc0}; 
 Point(pnum+1) = {x,y2,0,lc0}; 
 Point(pnum+2) = {x+GF,y1,0,lc0}; 
 Point(pnum+3) = {x+GF,y2,0,lc0}; 
 Line(lnum)={pnum,pnum+2};
 
 sA[] += {lnum};
 sR[] += {lnum};

 x=x+GF+WF;
 pnum=pnum+4;
 lnum=lnum+1;

EndFor

sEP[]={};  // lista positive elecrodes
sEG[]={};  // lista ground elecrodes

pnum=pnumI;
For t In {1:NF/2}    // creates positive elettrodi

  Line(lnum)={pnum+2,pnum+3};
  Line(lnum+1)={pnum+3,pnum+5};
  Line(lnum+2)={pnum+5,pnum+4};
  Line(lnum+3)={pnum+2,pnum+4};
  sA[] += {lnum,lnum+1,lnum+2};
  sR[] += {lnum+3};
  sEP[] += {lnum,lnum+1,lnum+2,lnum+3};
  pnum=pnum+8;
  lnum=lnum+4;

EndFor

pnum=pnumI+4;
For t In {1:NF/2}    // creates ground elettrodi

  Line(lnum)={pnum+2,pnum+3};
  Line(lnum+1)={pnum+3,pnum+5};
  Line(lnum+2)={pnum+5,pnum+4};
  Line(lnum+3)={pnum+2,pnum+4};
  sA[] += {lnum,lnum+1,lnum+2};
  sR[] += {lnum+3};
  sEG[] += {lnum,lnum+1,lnum+2,lnum+3};
  pnum=pnum+8;
  lnum=lnum+4;

EndFor

Line(lnum)={pnumI,4};
Line(lnum+1)={3,pnum-2};

Line Loop(lnum+2) = {-sA[],4,1,2,lnum,lnum+1};   // inner boundary Air
Plane Surface(snum) = {9,-(lnum+2)};  // surface for Air
nA=snum;

Line Loop(lnum+3) = {-sR[],-3,lnum,lnum+1};   // boundary resonator
Plane Surface(snum+1) = {lnum+3};  // surface for resonator
nR=snum+1;  

Physical Line (1) = {1,2,3,4};  // bottom electrode
Physical Line (2) = {sEP[]};  // positive electrode
Physical Line (3) = {sEG[]};  // ground
Physical Surface (4) = {nA};  // Air
Physical Surface (5) = {nR};  // Resonator



