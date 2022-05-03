file='Chap3_capacitance.msh';

% ANALYSIS TYPE
analysis.type=1;

% MATERIAL: Young and Poisson
material= [1.;
           9.5];

% SOLID: associates materials to sets
solid = [4 1;
         5 2];

% DBC: each row a bc. Physical set, direction, val
dbc= [2 1 1;
      3 1 0];

constr = [1];
  



