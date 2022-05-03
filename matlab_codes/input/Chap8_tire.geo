lc = .4;
rad1=2.;
rad2=3.;
ver=5;

Point(1) = {0.,ver,0.,lc};
Point(2) = {rad2,ver,0.,lc};
Point(3) = {0.,ver+rad2,0.,lc};
Point(4) = {-rad2,ver,0.,lc};
Point(5) = {0.,ver-rad2,0,lc};
Point(6) = {rad1,ver,0.,lc};
Point(7) = {0.,ver+rad1,0.,lc};
Point(8) = {-rad1,ver,0.,lc};
Point(9) = {0.,ver-rad1,0,lc};
Point(10) = {-4*rad2,-4*rad2,0,lc};
Point(11) = {rad2,rad2,0,lc};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};
Circle(5) = {6,1,7};
Circle(6) = {7,1,8};
Circle(7) = {8,1,9};
Circle(8) = {9,1,6};
Line(9) = {10,11};

Line Loop(10) = {1,2,3,4};
Line Loop(11) = {5,6,7,8};

Plane Surface(1) = {10,11};

Physical Surface(1) = {1}; 
Physical Line(2) = {9}; 