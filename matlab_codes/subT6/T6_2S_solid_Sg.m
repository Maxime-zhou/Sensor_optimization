
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Stresses at gauss points: T6 
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Sg=T6_2A_solid_Sg(X,mate,Ue)

E=mate(1);
nu=mate(2);
A=E/(1-nu^2)*[1,nu,0;
              nu,1,0;
              0,0,(1-nu)/2];   
                   
a_gauss=1/6*[4 1 1; 1 4 1; 1 1 4]; 
Sg=zeros(3,3);
for g=1:3,
  a=a_gauss(g,:);                            % area coord. for gauss point
  D=[4*a(1)-1 0 -4*a(3)+1 4*a(2)...         % derivative of shape functions...
      -4*a(2) 4*(a(3)-a(1));                 % w.r.t. a_1,a_2
      0 4*a(2)-1 -4*a(3)+1 4*a(1) ...
      4*(a(3)-a(2)) -4*a(1)]';
  J=X'*D;                                   % jacobian matrix
  detJ=J(1,1)*J(2,2)-J(1,2)*J(2,1);          % jacobian
  invJ=1/detJ*[ J(2,2) -J(1,2); ...          % inverse jacobian matrix
               -J(2,1)  J(1,1)];
  G=D*invJ;                                % gradient of shape functions
  B=[G(1,1) 0 G(2,1) 0 G(3,1) 0 ...
     G(4,1) 0 G(5,1) 0 G(6,1) 0;
     0 G(1,2) 0 G(2,2) 0 G(3,2)...
     0 G(4,2) 0 G(5,2) 0 G(6,2);
     G(1,2) G(1,1) G(2,2) G(2,1)...
     G(3,2) G(3,1) G(4,2) G(4,1)...
     G(5,2) G(5,1) G(6,2) G(6,1)];
  Sg(g,:)=(A*B*Ue)';
end

