file='Chap3_scubatank.msh';

% ANALYSIS TYPE
analysis.type = 6;

% MATERIAL: Young and Poisson
material = [100 0.2];

% SOLID: associates materials to sets
solid = [5 1];

% DBC: each row a bc. Physical set, direction, val
dbc = [1 2 0;
       2 1 0];

% TBC: each row a bc. Physical set, direction, val
tbc = [3 0 -1;
       4 0 -1];


