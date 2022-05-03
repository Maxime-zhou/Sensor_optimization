file='Chap5_tire.msh';

% ANALYSIS TYPE: plane strain
analysis.type = 4;

% MATERIAL: Young, Poisson, rho
material = [100 0.2 .2];

% SOLID: associates materials to sets
solid = [1 1];

% body forces val
bf = [0 -.5];

% parameters for contact
penalty = 20;   

% parameters for time history
tf = 10.;
Dt = .01;
output_interval=.1;

  
