file='Chap6_selective.msh';

% ANALYSIS TYPE
analysis.type = 4;

% MATERIAL: Young and Poisson
material = [1 0.4999];

% SOLID: associates materials to sets
solid = [4 1];

% DBC: each row a bc. Physical set, direction, val
dbc = [1 2 0;
       2 1 0;
       3 1 0;
       3 2 -.1];
  
