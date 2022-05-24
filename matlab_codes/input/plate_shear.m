file='MyGeometry_rect_coarse.msh';

% ANALYSIS TYPE
analysis.type = 5;

% MATERIAL: Young and Poisson
material = [10 0.3];
% E1 = 130Gpa
% E2 = 10Gpa
% G12 = 5Gpa
% nu12 = 0.3

% SOLID: associates materials to sets
solid = [4 1];

% DBC: each row a bc. Physical set, direction, val

% dbc = [1 1 0;
%        1 2 0;
%        2 1 0;
%        2 2 -0.1];
dbc = [1 1 0;
       1 2 0];
	   
% TBC: each row a bc. Physical set, direction, val
tbc = [2 2 0.1];
% tbc = [2 1 0.1];
 
   
  


