// global dimensions

H=0.5;
L=1;

lc1=H/2;

Point(1) = {-L, -H, 0, lc1};
Point(2) = {L, -H, 0, lc1};
Point(3) = {L, H, 0, lc1};
Point(4) = {-L, H, 0, lc1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(5) = {1, 2, 3, 4};
Plane Surface(5) = {5};

Physical Line (1) = {2};     // left line
Physical Line (2) = {4};     // right line
Physical Surface (3) = {5};  // surface
