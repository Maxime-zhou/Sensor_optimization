file='Chap3_hollowcyl.msh';

% ANALYSIS TYPE
analysis.type = 4;

% MATERIAL: Young and Poisson
material = [1 0];

% SOLID: associates materials to sets
solid = [4 1];

% DBC: each row a bc. node, direction, val
dbc = [1 2 0;
       2 1 0];  

% TBC: each row a bc. Physical set, direction, val
tbc = [3 0 1];




