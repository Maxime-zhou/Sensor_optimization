
function Ke=Q8_2D_therm_Ke(X,mate)

k=mate(1);
a_gauss=sqrt(3/5)*[-1 0 1];                  % Gauss abscissae
w_gauss=1/9*[5 8 5];                         % Gauss weights
Ke=zeros(8,8);
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
  G=(D*invJ)';                               % gradient of shape functions
  Ke=Ke+k*G'*G*detJ*w_gauss(g1)*w_gauss(g2);  % contribution to elemental matrix
 end
end