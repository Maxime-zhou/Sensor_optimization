file='Chap3_bulkmodulus.msh';

% ANALYSIS TYPE
analysis.type = 6;

% MATERIAL: Young and Poisson
material = [.01 0.;
            10. 0.];

% SOLID: associates materials to sets
solid = [4 1;
         5 2];

% TBC: each row a bc. Physical set, direction, val
tbc = [2 0 -1];

% DBC: each row a bc. Physical set, direction, val
dbc = [1 2 0;
       3 1 0];
 
   
  


