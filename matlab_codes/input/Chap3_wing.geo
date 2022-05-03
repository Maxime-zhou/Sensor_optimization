
lc1 = 5.;

Point(1) = {0.0,0.0,0.0,lc1};
Point(2) = {48.0,44.0,0.0,lc1};
Point(3) = {48.0,60.0,0.0,lc1};
Point(4) = {0.0,44.0,0.0,lc1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {1,2,3,4};

Plane Surface(1) = {5};

Physical Line(1) = {2}; 
Physical Line(2) = {4}; 
Physical Surface(3) = {1}; 
