// global dimensions

H=10;
V=1;

lc1=.25;

Point(1) = {0, 0, 0, lc1};
Point(2) = {H, 0, 0, lc1};
Point(3) = {H, V, 0, lc1};
Point(4) = {0, V, 0, lc1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(5) = {1,2,3,4};

// surface

Plane Surface(1) = {5};

Physical Line (1) = {4};  // left line blocked
Physical Line (2) = {2};  // right line loaded
Physical Surface (3) = {1};  // surface
