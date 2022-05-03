
file='Chap6_wedge.msh';

% ANALYSIS TYPE
analysis.type = 4;

% MATERIAL: Young and Poisson
material = [1 .4999];

% SOLID: associates materials to sets
solid = [3 1];

% DBC: each row a bc. Physical set, direction, val
dbc = [1 1 .5;
       2 1 0];
  
% nodal DBC: each row a bc. node, direction, val
dbcn = [1 1 0;
        1 2 0;
        10 1 0];  

% TBC: each row a bc. Physical set, direction, val
tbc = [1 1 1;
       2 1 -1];


