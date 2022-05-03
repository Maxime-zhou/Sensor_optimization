file='Chap8_tire.msh';

% ANALYSIS TYPE: plane strain
analysis.type = 4;

% MATERIAL: Young, Poisson, rho
material = [200 0.2 .3];

% SOLID: associates materials to sets
solid = [1 1];

% body forces val
bf = [0 -.25];

% parameters for contact
penalty = 20;   
f=.4;

% parameters for time history
tf = 10.;
Dt = .005;
output_interval=.2;

  
