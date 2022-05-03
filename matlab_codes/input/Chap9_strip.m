file='Chap9_strip.msh';

% ANALYSIS TYPE
analysis.type=4;

% Time history
LambdaT = [0:0.0021:.0126];

% MATERIAL: Young and Poisson - Sigma0 - Hardenning
material = [1000 0.1 1 0];

% SOLID: associates materials to sets
solid = [3 1];

% DBC: each row a bc. Physical set, direction, val
dbc = [1 1 1;
       2 1 0];

% nBC: each row a bc. Physical set, direction, val
dbcn = [1 1 0;
        1 2 0];
