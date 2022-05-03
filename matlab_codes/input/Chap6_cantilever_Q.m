
file='Chap6_cantilever_Q.msh';

% ANALYSIS TYPE
analysis.type=4;

% MATERIAL: Young and Poisson
material = [1000 0.0];

% SOLID: associates materials to sets
solid = [3 1];

% DBC: each row a bc. Physical set, direction, val
dbc = [1 1 0;
       1 2 0];

% TBC: each row a bc. Physical set, direction, val
tbc = [2 2 1];


