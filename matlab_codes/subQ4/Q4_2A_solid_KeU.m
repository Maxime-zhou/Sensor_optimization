
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Linear elastic stiffness matrix: Q4
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Ke=Q4_2A_solid_KeU(X,mate)

E=mate(1);
nu=mate(2);
A=E/((1+nu)*(1-2*nu))*[1-nu,nu,0;
                       nu,1-nu,0;
                       0,0,(1-2*nu)/2];  

D=1/4*[-1 1 1 -1;                % w.r.t. a_1,a_2
       -1 -1 1 1]';
J=X'*D;                                    % jacobian matrix
detJ=J(1,1)*J(2,2)-J(1,2)*J(2,1);          % jacobian
invJ=1/detJ*[ J(2,2) -J(1,2); ...          % inverse jacobian matrix
             -J(2,1)  J(1,1)];
G=D*invJ;                                  % gradient of shape functions
B=[G(1,1) 0 G(2,1) 0 G(3,1) 0 G(4,1) 0 ;
   0 G(1,2) 0 G(2,2) 0 G(3,2) 0 G(4,2) ;
   G(1,2) G(1,1) G(2,2) G(2,1)...
   G(3,2) G(3,1) G(4,2) G(4,1)];

Ke=B'*A*B*detJ*4; % contrib. to stiffness matrix


