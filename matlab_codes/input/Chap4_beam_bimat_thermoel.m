file='Chap4_beam_bimat.msh';

% ANALYSIS TYPE: 4 plane strain; 
analysis.type=4;

% MATERIAL: Young, Poisson and alpha
material= [1 0 2;
           1 0 1];

% SOLID: associates materials to sets
solid = [2 1;
         3 2];

% DBC: each row a bc. Physical set, direction, val
dbc= [1 1 0];
  
dbcn=[1 2 0];  





