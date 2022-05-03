lc1=.3;

H=2.;
V=10;

Point(1) = {0.0,0.0,0.0,lc1};
Point(2) = {0,-H,0.0,lc1};
Point(3) = {V,-H,0.0,lc1};
Point(4) = {V,0,0.0,lc1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {1,2,3,4};

Plane Surface(1) = {5};

Physical Line(1) = {1}; 
Physical Line(2) = {3}; 
Physical Surface(3) = {1}; 

