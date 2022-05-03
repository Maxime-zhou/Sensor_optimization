// global dimensions

H=10;
V=10;
R=2;

// mesh parameters

 
lc1=1;
lc2=0.5;


// point coordinates

Point(1) = {-H/2,-V/2,0,lc1};
Point(2) = {H/2,-V/2,0,lc1};
Point(3) = {H/2,V/2,0,lc1};
Point(4) = {-H/2,V/2,0,lc1};

Point(5) = {0,0,0,lc2};
Point(6) = {R,0,0,lc2};
Point(7) = {0,R,0,lc2};
Point(8) = {-R,0.,0,lc2};
Point(9) = {0,-R,0,lc2};

// lines and line loops

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {1,2,3,4};

Circle(6) = {6,5,7};
Circle(7) = {7,5,8};
Circle(8) = {8,5,9};
Circle(9) = {9,5,6};
Line Loop(10) = {-6,-7,-8,-9};

// surface

Plane Surface(1) = {5,10};

// physical entities

Physical Line(1) = {2}; 
Physical Line(2) = {4}; 
Physical Line(3) = {-6,-7,-8,-9}; 
Physical Surface(4) = {1}; 

