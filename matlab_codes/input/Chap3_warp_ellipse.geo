
lc1=.25;
a=2;
b=1;

Point(1) = {0,0,0,lc1};
Point(2) = {a,0,0,lc1};
Point(3) = {0,b,0,lc1};
Point(4) = {-a,0,0,lc1};
Point(5) = {0,-b,0,lc1};

Ellipse(1) = {2,1,1,3};
Ellipse(2) = {3,1,1,4};
Ellipse(3) = {4,1,1,5};
Ellipse(4) = {5,1,1,2};
Line Loop(5) = {1,2,3,4};

Plane Surface(1) = {5};

Physical Line(1) = {1,2,3,4}; 
Physical Surface(2) = {1}; 
