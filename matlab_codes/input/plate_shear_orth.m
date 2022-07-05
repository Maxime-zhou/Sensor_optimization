file ='MyGeometry_rect_ddx000625.msh';
% file = 'MyGeometry_rect_fine.msh';

% ANALYSIS TYPE
% analysis.type = 5;
analysis.type = 8;
% MATERIAL: Young and Poisson
material = [130e3 10e3 0.3 5e3];
% material = [10 0.3];
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
% tbc = [2 2 0.1]; % y direction
tbc = [2 1 0.1]; % x direction 
 
   
  


