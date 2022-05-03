lc1 = .5;
lc2 = .25;
lc3 = .1;


Point(1) = {-8.0,-3.0,0.0,lc1};
Point(2) = {-2.0,-3.0,0.0,lc2};
Point(3) = {0.0,-1.0,0.0,lc3};
Point(4) = {2.0,-3.0,0.0,lc2};
Point(5) = {8.0,-3.0,0.0,lc1};
Point(6) = {8.0,3.0,0.0,lc1};
Point(7) = {2.0,3.0,0.0,lc2};
Point(8) = {0.0,1.0,0.0,lc3};
Point(9) = {-2.0,3.0,0.0,lc2};
Point(10) = {-8.0,3.0,0.0,lc1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,9};
Line(9) = {9,10};
Line(10) = {10,1};
Line Loop(11) = {1,2,3,4,5,6,7,8,9,10};

Plane Surface(1) = {11};

Physical Line(1) = {5}; 
Physical Line(2) = {10}; 
Physical Surface(3) = {1}; 
