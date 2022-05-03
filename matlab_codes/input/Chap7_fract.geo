H=5;  //  semiwidth of plate
V=5;  //  semiheight of plate
a=1;  //  semilength of crack  (center crack is at x1=x2=0)

lc1=.5;
lc2=.2;
lc3=.025;  // element size at crack tip

Point(1) = {-H,-V,0.0,lc1};
Point(2) = {H,-V,0.0,lc1};
Point(3) = {H,0.0,0.0,lc2};
Point(4) = {H,V,0.0,lc1};
Point(5) = {-H,V,0.0,lc1};
Point(6) = {-H,0.0,0.0,lc2};
Point(7) = {-a,0.0,0.0,lc3};
Point(8) = {a,0.0,0.0,lc3};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,8};
Line(4) = {8,7};
Line(5) = {7,6};
Line(6) = {6,1};
Line(7) = {3,4};
Line(8) = {4,5};
Line(9) = {5,6};
Line(10) = {7,8};
Line Loop(11) = {1,2,3,4,5,6};    //  border of lower surface
Line Loop(12) = {7,8,9,-5,10,-3}; //  border of upper surface

Plane Surface(1) = {11};  // lower surface
Plane Surface(2) = {12};  // upper surface 

Physical Line(1) = {1};      // lower line for traction
Physical Line(2) = {8};      // upper line for traction
Physical Surface(3) = {1,2}; 
