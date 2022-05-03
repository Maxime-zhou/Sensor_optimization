
L=5;
lc1=.25;

Point(1) = {0,0,0,lc1};
Point(2) = {L,0,0,lc1};
Point(3) = {L,L,0,lc1};
Point(4) = {0,L,0,lc1};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {1,2,3,4};

Plane Surface(1) = {5};


Physical Line(1) = {1,2,3,4}; 
Physical Surface(2) = {1}; 
