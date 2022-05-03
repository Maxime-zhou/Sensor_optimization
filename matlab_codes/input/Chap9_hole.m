file='Chap9_hole.msh'

% ANALYSIS TYPE
analysis.type = 4;

% Time history
LambdaT = [0. 1. 2. 3. 4. 5. 6. 7. 8.];

% MATERIAL: Young and Poisson - sigma0 - hardening 
material = [200000 0.3 140 0.];

% SOLID: associates materials to sets
solid = [4 1];

% DBC: each row a bc. Physical set, direction, val
dbc = [1 2 0;
       3 1 0;
	   2 1 .002];

% nodal DBC: each row a bc. node, direction, val
%dbcn = [];

% TBC: each row a bc. Physical set, direction, val
%tbc = [];