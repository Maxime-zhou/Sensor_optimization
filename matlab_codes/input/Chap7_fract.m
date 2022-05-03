file='Chap7_fract.msh';

% ANALYSIS TYPE
analysis.type = 4;

% MATERIAL: Young and Poisson
material = [1 1/3];

% SOLID: associates materials to sets
solid = [3 1];

% nodal DBC: each row a bc. node, direction, val
dbcn = [1 1 0;
        1 2 0;
        2 2 0];  

% TBC: each row a bc. Physical set, direction, val
tbc = [1 2 -1; 
       2 2 1];

