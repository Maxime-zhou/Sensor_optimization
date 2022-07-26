h = 1;
w = h/2;
c = 0.5;

//Points
Point(1) = {0, 0, 0, c};
Point(2) = {w, 0, 0, c};
Point(3) = {w, h, 0, c};
Point(4) = {0, h, 0, c};


//Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

//Transfinite Curve{1} = 20;
//Transfinite Curve{2} = 38;
//Transfinite Curve{3} = 20;
//Transfinite Curve{4} = 38;


Transfinite Curve{1} = 41;
Transfinite Curve{2} = 81;
Transfinite Curve{3} = 41;
Transfinite Curve{4} = 81;
//Line Loops

Line Loop(1) = {1, 2, 3, 4};

//Surfaces 
Plane Surface(1) = {1};

Transfinite Surface(1);

//Generate Quadrangles mesh
Recombine Surface{1};

Physical Line(1) = {1}; 
Physical Line(2) = {3}; 
Physical Line(3) = {2,4}; 
Physical Surface(4) = {1}; 


