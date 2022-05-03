file='Chap8_buckling.msh'

% ANALYSIS TYPE (PlaneStrain)
analysis.type = 4;

% load sequence
delta = 1/15;
LambdaT = [0.:delta:1.];

% MATERIAL: Young and Poisson
material = [1 .0];

% SOLID: associates materials to sets
solid = [4 1];

% DBC: each row a bc. Physical set, direction, val
dbc = [1 1  0.;
	   1 2  0.;
	   2 1  0.;
       2 2 -0.3];

% TBC: each row a bc. Physical set, direction, val
tbc = [3 1  0.00001];


