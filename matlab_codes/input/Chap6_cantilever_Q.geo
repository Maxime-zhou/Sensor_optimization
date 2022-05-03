// global dimensions

H=10;
V=1;
num=2;

lc1=V/num;

Point(1) = {0, 0, 0, lc1};
Point(2) = {0, V, 0, lc1};

Line(1) = {2, 1};

Extrude {H, 0, 0} { 
  Line{1};  Layers{H*num}; Recombine; 
}

Physical Line (1) = {1};  // left line blocked
Physical Line (2) = {2};  // right line loaded
Physical Surface (3) = {5};  // surface
