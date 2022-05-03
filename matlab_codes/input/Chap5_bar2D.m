file='Chap5_bar2D.msh';

% ANALYSIS TYPE: plane strain
analysis.type = 4;

% MATERIAL: Young, Poisson, rho
material = [1 0. 1];

% SOLID: associates materials to sets
solid = [3 1];

% DBC on nodes
dbcn = [1 2 0];

% DBC: each row a bc. Physical set, direction, val
dbc = [1 1 0];

% TBC: each row a bc. Physical set, direction, val
tbc = [2 1 1];

% parameters for time history
tf =40.;
Dt = .02;
output_interval=.5;
   
  


