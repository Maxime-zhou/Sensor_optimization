// global dimensions

R=1;

// mesh parameters
 
lcR=.3;
lc0=.1;

// point coordinates

Point(1) = {0,0,0,lc0};
Point(2) = {R,0,0,lcR};
Point(3) = {0,R,0,lcR};

// lines and line loops

Line(1) = {1,2};
Circle(2) = {2,1,3};
Line(3) = {3,1};
Line Loop(4) = {1,2,3};

// surface

Plane Surface(1) = {4};

// physical entities

Physical Line(1) = {2};      // to impose temp
Physical Surface(2) = {1}; 

