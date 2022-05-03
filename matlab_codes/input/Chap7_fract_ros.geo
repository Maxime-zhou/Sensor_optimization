Mesh.Algorithm=2;  // 1 meshadapt algorithm
Mesh.ElementOrder=2;

a=1.;
rad=0.02;
V=10;
H=5;

pi=4*Atan(1.0);
num=8;
ang=pi/num;

lc1=0.015;
lc2=0.4;

Point(1) = {a,0.0,0.0,rad};
Point(2) = {a-rad,0.0,0.0,rad};
Line(1)={2,1};

i = 0 ;
ang=0;
For(1:num)
  i+=1;
  ang=ang+pi/num; 
  Point(2+i)={a-rad*Cos(ang),-rad*Sin(ang),0.0,rad};
  Circle(2*i) = {1+i,1,2+i};
  Line(2*i+1)={2+i,1};
  Line Loop(i) = {-(2*i-1),2*i,2*i+1};
  Plane Surface(i) = {i};
  lineloop[i-1] = -2*i ;
  surfloop[i-1] = i ;
EndFor

Point(num+3) = {0.0,0.0,0.0,rad};
Point(num+4) = {0.0,-V,0.0,lc2};
Point(num+5) = {H,-V,0.0,lc2};
Point(num+6) = {H,0.0,0.0,lc2};

Line(2*num+2)={2,num+3};
Line(2*num+3)={num+3,num+4};
Line(2*num+4)={num+4,num+5};
Line(2*num+5)={num+5,num+6};
Line(2*num+6)={num+6,num+2};

Line Loop(num+1) = {2*num+2,2*num+3,2*num+4,2*num+5,2*num+6, lineloop[]};
Plane Surface(num+1) = {num+1};


// physical entities

Physical Line(1) = {2*num+3};   // simmetria verticale
Physical Line(2) = {2*num+4};   // load inferiore
Physical Line(3) = {2*num+6, 2*num+1};   // simmetria hor
Physical Surface(4) = {num+1, surfloop[]}; 
