H=1;
L=10;

lc1=.25;

Point(1) = {H,H,0,lc1};
Point(2) = {-H,H,0,lc1};
Point(3) = {-H,-H,0,lc1};
Point(4) = {H,-H,0,lc1};

Point(5) = {H,H,L,lc1};
Point(6) = {-H,H,L,lc1};
Point(7) = {-H,-H,L,lc1};
Point(8) = {H,-H,L,lc1};

// lower surf

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(1) = {-1,-2,-3,-4};
Plane Surface(1) = {1};

// upper surf

Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,5};
Line Loop(2) = {5,6,7,8};
Plane Surface(2) = {2};

// lateral

Line(9) = {1,5};
Line(10) = {2,6};
Line(11) = {3,7};
Line(12) = {4,8};

Line Loop(3) = {9,5,-10,-1};
Plane Surface(3) = {3};
Line Loop(4) = {10,6,-11,-2};
Plane Surface(4) = {4};
Line Loop(5) = {11,7,-12,-3};
Plane Surface(5) = {5};
Line Loop(6) = {12,8,-9,-4};
Plane Surface(6) = {6};

// FEM volume

Surface Loop(1) = {1,2,3,4,5,6};
Volume(1) = {1};

Physical Surface(1) = {1};   // below
Physical Surface(2) = {2};   // upper
Physical Volume(3) = {1}; // volume


