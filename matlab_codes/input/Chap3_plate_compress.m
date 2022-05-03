file='Chap3_plate_compress.msh';

% ANALYSIS TYPE
analysis.type = 4;

% MATERIAL: Young and Poisson
material = [10 0.3];

% SOLID: associates materials to sets
solid = [4 1];

% DBC: each row a bc. Physical set, direction, val
dbc = [1 1 0;
       1 2 0;
       2 1 0;
       2 2 -0.1];
 
   
  


