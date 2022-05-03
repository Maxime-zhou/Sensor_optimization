
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% element stiffness matrix for 3 node line element and spherical symmetry
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Ke=B3_1S_solid_Ke(X,lambda,mu);

a_gauss=sqrt(3/5)*[-1 0 1];                  % 3 point Gauss abscissae
w_gauss=1/9*[5 8 5];                         % 3 point Gauss weights
Ke=zeros(3,3);
for g=1:length(a_gauss),                     % loop over Gauss points
 a=a_gauss(g);                               % param. coordinates for gauss point
 N=[.5*a*(a-1) (1-a^2) .5*a*(1+a)]';
 r=N'*X;
 D=[a-.5 -2*a a+.5]';
 J=X'*D;                                    % jacobian matrix
 G=D/J;                                     % gradient of shape functions
 Ke=Ke+((lambda+2*mu)*G*G'*r^2+...          % element stiffness matrix 
        2*lambda*(G*N'+N*G')*r+...
        4*(lambda+mu)*N*N')*J*w_gauss(g);  
end
