lc1=.25;

L=10.;
H=1;

Point(1) = {0.0,0.0,0.0,lc1};
Point(2) = {0.0,-H,0.0,lc1};
Point(3) = {L,-H,0.0,lc1};
Point(4) = {L,0,0.0,lc1};
Point(5) = {L,H,0.0,lc1};
Point(6) = {0.0,H,0.0,lc1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line(5) = {4,5};
Line(6) = {5,6};
Line(7) = {6,1};

Line Loop(8) = {1,2,3,4};
Line Loop(9) = {-4,5,6,7};

Plane Surface(1) = {8};
Plane Surface(2) = {9};

Physical Line(1) = {1,7}; 
Physical Surface(2) = {1}; 
Physical Surface(3) = {2}; 

