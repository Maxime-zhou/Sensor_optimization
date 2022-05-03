file='Chap3_warp_square.msh';

% ANALYSIS TYPE
analysis.type=1;

% MATERIAL: Young and Poisson
material= [1];

% SOLID: associates materials to sets
solid = [2 1];

% DBC: each row a bc. Physical set, direction, val
dbc= [1 1 0];

% body foces; elset and val
bf= 1;




