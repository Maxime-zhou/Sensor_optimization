file='Chap3_wing.msh';

% ANALYSIS TYPE
analysis.type=5;

% MATERIAL: Young and Poisson
material= [1 1/3];

% SOLID: associates materials to sets
solid = [3 1];

% DBC: each row a bc. Physical set, direction, val
dbc= [2 1 0;
      2 2 0];

% TBC: each row a bc. Physical set, direction, val
tbc = [1 2 0.0625];




