
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Linear elastic stiffness matrix: Q4
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Ke=Q4_2S_solid_Ke(X,mate)

E=mate(1);
nu=mate(2);
A=E/(1-nu^2)*[1,nu,0;
              nu,1,0;
              0,0,(1-nu)/2];   

a_gauss=1/sqrt(3)*[-1 1];                    % Gauss abscissae
w_gauss=[1 1];                               % Gauss weights
Ke=zeros(8,8);
for g1=1:2,                                  % loop over Gauss points
 a1=a_gauss(g1);                             % param. coord. of gauss point
 for g2=1:2,                                 % loop over Gauss points
  a2=a_gauss(g2);                            % param. coord. for gauss point
  D=1/4*[-(1-a2) (1-a2) (1+a2) -(1+a2) ;     % w.r.t. a_1,a_2
         -(1-a1) -(1+a1) (1+a1) (1-a1)]';
  J=X'*D;                                    % jacobian matrix
  detJ=J(1,1)*J(2,2)-J(1,2)*J(2,1);          % jacobian
  invJ=1/detJ*[ J(2,2) -J(1,2); ...          % inverse jacobian matrix
               -J(2,1)  J(1,1)];
  G=D*invJ;                                  % gradient of shape functions
  B=[G(1,1) 0 G(2,1) 0 G(3,1) 0 G(4,1) 0 ;
     0 G(1,2) 0 G(2,2) 0 G(3,2) 0 G(4,2) ;
     G(1,2) G(1,1) G(2,2) G(2,1)...
     G(3,2) G(3,1) G(4,2) G(4,1)];
  Ke=Ke+B'*A*B*detJ*w_gauss(g1)*w_gauss(g2); % contrib. to stiffness matrix
 end
end

