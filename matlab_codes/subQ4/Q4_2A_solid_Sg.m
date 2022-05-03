
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Stresses at gauss points: Q4 
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Sg=Q4_2A_solid_Sg(X,mate,Ue)

E=mate(1);
nu=mate(2);
A=E/((1+nu)*(1-2*nu))*[1-nu,nu,0;
                       nu,1-nu,0;
                       0,0,(1-2*nu)/2];  

a_gauss=1/sqrt(3)*[-1 1];                    % Gauss abscissae
Sg=zeros(3,3);
g=0;   
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
  B=[G(1,1) 0 G(2,1) 0 G(3,1) 0 G(4,1) 0;
     0 G(1,2) 0 G(2,2) 0 G(3,2) 0 G(4,2);
     G(1,2) G(1,1) G(2,2) G(2,1)...
     G(3,2) G(3,1) G(4,2) G(4,1)];
  g=g+1;
  Sg(g,:)=(A*B*Ue)';
 end 
end