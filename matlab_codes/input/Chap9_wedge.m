file='Chap9_wedge.msh';

% ANALYSIS
analysis.type = 4;

% Time history
%LambdaT = [0. 0.08 0.16 0.24 0.32 0.4 0.46 0.52 0.57 0.59 0.6];
%LambdaT = [0. 0.08 0.16 0.24 0.32 0.4 0.48 0.56 0.64 0.72 0.8];

LambdaT = [0. 0.08 0.16 0.24 0.32 0.4 0.46];

% MATERIAL: Young, Poisson, Sigma0, Hardening, 
%material = [1 .3 0.88 0.];
material = [1 .3 0.88 0.05];

% SOLID: associates materials to sets
solid = [3 1];

% DBC: each row a bc. Physical set, direction, val
%dbc = [];
  
% nodal DBC: each row a bc. node, direction, val
dbcn = [1 1 0;
        1 2 0;
	   10 1 0];  

% TBC: each row a bc. Physical set, direction, val
tbc = [1 1  1.;
       2 1 -1.];

