file='Chap3_heat3D.msh';

% ANALYSIS TYPE: 3D potential problem
analysis.type = 3; 

% MATERIAL: conductivity
material = [1];

% SOLID: associates materials to sets
solid = [3 1];

% DBC: each row a bc. Physical set, direction, val
dbc = [1 1 0;
       2 1 10];
  




