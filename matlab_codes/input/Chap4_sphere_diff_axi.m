file='Chap4_sphere_diff_axi.msh';

% ANALYSIS TYPE
analysis.type = 2;

% MATERIAL: k and rhoc
material = [1 10];

% SOLID: associates materials to sets
solid = [2 1];

% DBC: each row a bc. Physical set, direction, val
dbc = [1 1 1];

% parameters for time history
tf =1;
Dt = 0.001;
output_interval=.1;

  


