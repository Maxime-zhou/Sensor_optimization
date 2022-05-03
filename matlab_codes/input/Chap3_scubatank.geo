
lc1 = 0.15;
h=20;
r=10;
e=.5;

Point(1) = {r,0,0,lc1};
Point(2) = {r+e,0,0,lc1};
Point(3) = {r+e,h,0,lc1};
Point(4) = {0,r+h+e,0,lc1};
Point(5) = {0,r+h,0,lc1};
Point(6) = {r,h,0,lc1};
Point(7) = {0,h,0,lc1};

Line(1) = {1,2};
Line(2) = {2,3};
Circle(3) = {3,7,4};
Line(4) = {4,5};
Circle(5) = {5,7,6};
Line(6) = {6,1};
Line Loop(7) = {1,2,3,4,5,6};

Plane Surface (1) = {7};

Physical Line (1) = {1};
Physical Line (2) = {4};
Physical Line (3) = {5};
Physical Line (4) = {6};
Physical Surface (5) = {1};
