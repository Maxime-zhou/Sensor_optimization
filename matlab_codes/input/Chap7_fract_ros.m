file='Chap7_fract_ros.msh';

% ANALYSIS TYPE
analysis.type = 4;

% MATERIAL: Young and Poisson
material = [1 1/3];

% SOLID: associates materials to sets
solid = [4 1];

% nodal DBC: each row a bc. node, direction, val
dbc = [1 1 0;
       3 2 0];  

% TBC: each row a bc. Physical set, direction, val
tbc = [2 2 -1];


