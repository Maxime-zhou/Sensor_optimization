
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Linear elastic stiffness matrix: T4
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Ke=Q4_2A_solid_KeS(X,mate)

E=mate(1);
nu=mate(2);
mu=E/(2*(1+nu));
Ke=zeros(8,8);

AI=[2*mu,0,0;                                 % first integrate mu eps:eps
   0,2*mu,0;                                 % contribution
   0,0,mu];  
a_gauss=1/sqrt(3)*[-1 1];                    % Gauss abscissae
w_gauss=[1 1];                               % Gauss weights
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
  Ke=Ke+B'*AI*B*detJ*w_gauss(g1)*w_gauss(g2); % contrib. to stiffness matrix
 end
end

lambda=2*mu*nu/(1-2*nu);                   % adds lambda treps contribution
AJ=lambda*[1,1,0;
          1,1,0;
          0,0,0];  
D=1/4*[-1 1 1 -1;                          % der. of shape functions in ..
       -1 -1 1 1]';                        % a-1=a_2=0
J=X'*D;                                    % jacobian matrix
detJ=J(1,1)*J(2,2)-J(1,2)*J(2,1);          % jacobian
invJ=1/detJ*[ J(2,2) -J(1,2); ...          % inverse jacobian matrix
             -J(2,1)  J(1,1)];
G=D*invJ;                                  % gradient of shape functions
B=[G(1,1) 0 G(2,1) 0 G(3,1) 0 G(4,1) 0 ;
   0 G(1,2) 0 G(2,2) 0 G(3,2) 0 G(4,2) ;
   G(1,2) G(1,1) G(2,2) G(2,1)...
   G(3,2) G(3,1) G(4,2) G(4,1)];
Ke=Ke+B'*AJ*B*detJ*4;                       % contrib. to stiffness matrix


