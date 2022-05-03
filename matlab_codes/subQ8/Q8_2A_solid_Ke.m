
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Linear elastic stiffness matrix: Q8
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Ke=Q8_2A_solid_Ke(X,mate)

E=mate(1);
nu=mate(2);
A=E/((1+nu)*(1-2*nu))*[1-nu,nu,0;
                       nu,1-nu,0;
                       0,0,(1-2*nu)/2];  
  
a_gauss=sqrt(3/5)*[-1 0 1];                  % Gauss abscissae
w_gauss=1/9*[5 8 5];                         % Gauss weights
Ke=zeros(16,16);
for g1=1:3,                                  % loop over Gauss points
 a1=a_gauss(g1);                             % param. coord. of gauss point
 for g2=1:3,                                 % loop over Gauss points
  a2=a_gauss(g2);                            % param. coord. for gauss point
  D=1/4*[(1-a2)*(2*a1+a2) -(2*a1-a2)*(-1+a2) (1+a2)*(2*a1+a2)  (1+a2)*(2*a1-a2) ...
          4*a1*(-1+a2) 2*(1-a2^2) -4*a1*(1+a2) -2*(1-a2^2);
         (1-a1)*(2*a2+a1) -(a1-2*a2)*(1+a1) (1+a1)*(a1+2*a2)  (-1+a1)*(a1-2*a2) ...
          -2*(1-a1^2) -4*(1+a1)*a2 2*(1-a1^2) 4*(-1+a1)*a2]';
  J=X'*D;                                    % jacobian matrix
  detJ=J(1,1)*J(2,2)-J(1,2)*J(2,1);          % jacobian
  invJ=1/detJ*[ J(2,2) -J(1,2); ...          % inverse jacobian matrix
               -J(2,1)  J(1,1)];
  G=D*invJ;                                  % gradient of shape functions
  B=[G(1,1) 0 G(2,1) 0 G(3,1) 0 G(4,1) 0 G(5,1) 0 G(6,1) 0 G(7,1) 0 G(8,1) 0 
     0 G(1,2) 0 G(2,2) 0 G(3,2) 0 G(4,2) 0 G(5,2) 0 G(6,2) 0 G(7,2) 0 G(8,2);
     G(1,2) G(1,1) G(2,2) G(2,1) G(3,2) G(3,1) G(4,2) G(4,1) ...
     G(5,2) G(5,1) G(6,2) G(6,1) G(7,2) G(7,1) G(8,2) G(8,1)];
  Ke=Ke+B'*A*B*detJ*w_gauss(g1)*w_gauss(g2); % contrib. to stiffness matrix
 end
end

