
function Fe=T6_2D_therm_Fe(X,val)

a_gauss=1/6*[4 1 1; 1 4 1; 1 1 4];           % Gauss abscissae
w_gauss=[1/6 1/6 1/6];                       % Gauss weights
Fe=zeros(6,1);
for g=1:3,                                   % loop over Gauss points
 a=a_gauss(g,:);                            % param. coordinates for gauss point
 D=[4*a(1)-1 0 -4*a(3)+1 4*a(2)...          % derivative of shape functions...
    -4*a(2) 4*(a(3)-a(1));                  % w.r.t. a_1,a_2
    0 4*a(2)-1 -4*a(3)+1 4*a(1) ...
    4*(a(3)-a(2)) -4*a(1)]';
 J=X'*D;                                    % jacobian matrix
 detJ=J(1,1)*J(2,2)-J(1,2)*J(2,1);          % jacobian
 N=[a(1)*(2*a(1)-1) a(2)*(2*a(2)-1) ...      % list of shape functions
   a(3)*(2*a(3)-1) 4*a(1)*a(2) ...
   4*a(2)*a(3) 4*a(1)*a(3)];
 Fe=Fe+val*detJ*N'*w_gauss(g);              % contribution to stiffness matrix
end

