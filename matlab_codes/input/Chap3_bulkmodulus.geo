// global dimensions

H=10;
V=10;
R=5;

// mesh parameters

 
lc1=.125;
lc2=1.;


// point coordinates

Point(1) = {0,0,0,lc2};
Point(2) = {R,0,0,lc1};
Point(3) = {0,R,0,lc1};

Point(4) = {H,0,0,lc2};
Point(5) = {H,V,0,lc2};
Point(6) = {0,V,0,lc2};

// lines and line loops

Line(1) = {1,2};
Circle(2) = {2,1,3};
Line(3) = {3,1};
Line Loop(4) = {1,2,3};

Line(5) = {2,4};
Line(6) = {4,5};
Line(7) = {5,6};
Line(8) = {6,3};
Line Loop(9) = {-2,5,6,7,8};

// surface

Plane Surface(1) = {4};
Plane Surface(2) = {9};

// physical entities

Physical Line(1) = {1,5};  // lower disp
Physical Line(2) = {6,7};  // pressure
Physical Line(3) = {8,3};  // axi
Physical Surface(4) = {1}; // inclusion
Physical Surface(5) = {2}; // rest
 

