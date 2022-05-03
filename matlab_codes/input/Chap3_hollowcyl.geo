// global dimensions

R1=1;
R2=2;

// mesh parameters

h=.05;

// point coordinates

Point(1) = {0,0,0,h};
Point(2) = {R1,0,0,h};
Point(3) = {R2,0,0,h};
Point(4) = {0,R2,0,h};
Point(5) = {0,R1,0,h};

// lines and line loops

Line(1) = {2,3};
Circle(2) = {3,1,4};
Line(3) = {4,5};
Circle(4) = {2,1,5};
Line Loop(5) = {1,2,3,-4};

// surface

Plane Surface(1) = {5};

// physical entities

Physical Line(1) = {1}; 
Physical Line(2) = {3}; 
Physical Line(3) = {4}; 
Physical Surface(4) = {1}; 

