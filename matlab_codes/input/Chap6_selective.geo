// global dimensions

H=1;
V=1;
num=16;

lc1=V/num;

Point(1) = {0, 0, 0, lc1};
Point(2) = {0, V, 0, lc1};

Line(1) = {2, 1};

Extrude {H, 0, 0} { 
  Line{1};  Layers{H*num}; Recombine; 
}

Physical Line (1) = {4};  // lower line
Physical Line (2) = {2};  // right line
Physical Line (3) = {3};  // upper line
Physical Surface (4) = {5};  // surface
