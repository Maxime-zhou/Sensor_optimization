file='Appendix_example.msh';

% ANALYSIS TYPE: plane strain
analysis.type = 4;

% MATERIAL: Young and Poisson
material = [1 1/3];

% SOLID: associates materials to sets
solid = [4 1];

% DBC: each row a bc. Physical set, direction, val
dbc = [1 1 0.5;
       2 1 0];
  
% nodal DBC: each row a bc. node, direction, val
dbcn = [1 2 0];  

% TBC: each row a bc. Physical set, direction, val
tbc = [3 0 -2];




