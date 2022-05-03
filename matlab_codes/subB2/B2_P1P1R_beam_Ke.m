
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% element stiffness matrix for linear 2 nodes beam element (u, phi)
% reduced integration scheme
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function Ke=B2_P1P1R_beam_Ke(X,EI,kmuS);

r1=X(1); r2=X(2);                            % radial coordinates of nodes
h=r2-r1;                                     % length of element
Ke=(EI/h)*[0,  0, 0,  0;...                  % element stiffness matrix 
     0,  1, 0, -1;...
     0,  0, 0,  0;...
     0, -1, 0,  1];
Ke=Ke+kmuS/(4*h)*[2*2,-2*h, -2*2, -2*h;...  
            -2*h, h*h,  2*h,  h*h;...
		     -2*2, 2*h,  2*2,  2*h;...
			 -2*h, h*h,  2*h,  h*h];
