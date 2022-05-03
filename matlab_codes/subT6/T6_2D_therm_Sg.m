

function Sg=T6_2D_therm_Sg(X,mate,Ue)

k=mate(1);

a_gauss=1/6*[4 1 1; 1 4 1; 1 1 4];          % Gauss abscissae
Sg=zeros(3,2);
for g=1:3,                                  % loop over Gauss points
 a=a_gauss(g,:);                            % param. coordinates for gauss point
 D=[4*a(1)-1 0 -4*a(3)+1 4*a(2)...          % derivative of shape functions...
    -4*a(2) 4*(a(3)-a(1));                  % w.r.t. a_1,a_2
    0 4*a(2)-1 -4*a(3)+1 4*a(1) ...
    4*(a(3)-a(2)) -4*a(1)]';
 J=X'*D;                                    % jacobian matrix
 detJ=J(1,1)*J(2,2)-J(1,2)*J(2,1);          % jacobian
 invJ=1/detJ*[ J(2,2) -J(1,2); ...          % inverse jacobian matrix
              -J(2,1)  J(1,1)];
 G=(D*invJ)';                               % gradient of shape functions
 Sg(g,:)=(k*G*Ue)';
end

