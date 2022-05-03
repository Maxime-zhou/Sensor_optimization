lc1=.125;
lc2=.25;

L=10.;
Ld2=L/2.;
H=5;

Point(1) = {0.0,0.0,0.0,lc1};
Point(2) = {H,0.0,0.0,lc1};
Point(3) = {H,Ld2,0.0,lc2};
Point(4) = {H,L,0.0,lc1};
Point(5) = {0.0,L,0.0,lc1};
Point(6) = {0.0,Ld2,0.0,lc2};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};
Line Loop(7) = {1,2,3,4,5,6};

Plane Surface(1) = {7};

Physical Line(1) = {1}; 
Physical Line(2) = {4}; 
Physical Line(3) = {5,6}; 
Physical Surface(4) = {1}; 

