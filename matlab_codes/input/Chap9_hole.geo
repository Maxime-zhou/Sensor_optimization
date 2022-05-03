lc = .5;
rad=.3;

Point(1) = {0.0,0.0,0.0,lc};
Point(2) = {rad,0.0,0.0,.2*lc};
Point(3) = {5.0,0.0,0.0,lc};
Point(4) = {5.0,5.0,0.0,lc};
Point(5) = {0,5,0,lc};
Point(6) = {0,rad,0,.2*lc};

Line(1) = {2,3};
Line(2) = {3,4};
Line(3) = {4,5};
Line(4) = {5,6};
Circle(5) = {6,1,2};
Line Loop(6) = {1,2,3,4,5};

Plane Surface(1) = {6};

Physical Line(1) = {1}; 
Physical Line(2) = {2}; 
Physical Line(3) = {4}; 
Physical Surface(4) = {1}; 
