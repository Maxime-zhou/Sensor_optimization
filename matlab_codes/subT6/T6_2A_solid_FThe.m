
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% generalised nodal forces due to temperature distribution
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Fe=T6_2A_solid_FTe(X,mate,Te)

E=mate(1);
nu=mate(2);
alpha=mate(3);

a_gauss=1/6*[4 1 1; 1 4 1; 1 1 4];           % Gauss abscissae
w_gauss=[1/6 1/6 1/6];                       % Gauss weights
Fe=zeros(12,1);
for g=1:3,                                   % loop over Gauss points
  a=a_gauss(g,:);                            % param. coordinates for gauss point
  N=[a(1)*(2*a(1)-1) a(2)*(2*a(2)-1) ...      % list of shape functions
   a(3)*(2*a(3)-1) 4*a(1)*a(2) ...
   4*a(2)*a(3) 4*a(1)*a(3)];
  T=N*Te';
  D=[4*a(1)-1 0 -4*a(3)+1 4*a(2)...          % derivative of shape functions...
     -4*a(2) 4*(a(3)-a(1));                  % w.r.t. a_1,a_2
     0 4*a(2)-1 -4*a(3)+1 4*a(1) ...
     4*(a(3)-a(2)) -4*a(1)]';
  J=X'*D;                                    % jacobian matrix
  detJ=J(1,1)*J(2,2)-J(1,2)*J(2,1);          % jacobian
  invJ=1/detJ*[ J(2,2) -J(1,2); ...          % inverse jacobian matrix
               -J(2,1)  J(1,1)];
  G=D*invJ;                                  % gradient of shape functions
  B=[G(1,1) G(1,2) G(2,1) G(2,2) G(3,1) G(3,2) ...
     G(4,1) G(4,2) G(5,1) G(5,2) G(6,1) G(6,2)];
  Fe=Fe+B'*alpha*E/(1-2*nu)*T*detJ*w_gauss(g);  % nodal forces
end

