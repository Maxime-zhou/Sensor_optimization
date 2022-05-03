
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Element contributution to Gtheta: T6
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function GBe=T6_2A_solid_gtheta(X,mate,Ue,Ge)

E=mate(1);
nu=mate(2);
A=E/((1+nu)*(1-2*nu))*[1-nu,nu,0;
                       nu,1-nu,0;
                       0,0,(1-2*nu)/2]; 

GBe=0;
a_gauss=1/6*[4 1 1; 1 4 1; 1 1 4];           % Gauss abscissae
w_gauss=[1/6 1/6 1/6];                       % Gauss weights
for g=1:3,                                   % loop over Gauss points
  a=a_gauss(g,:);                            % param. coordinates for gauss point
  D=[4*a(1)-1 0 -4*a(3)+1 4*a(2)...          % derivative of shape functions...
     -4*a(2) 4*(a(3)-a(1));                  % w.r.t. a_1,a_2
     0 4*a(2)-1 -4*a(3)+1 4*a(1) ...
     4*(a(3)-a(2)) -4*a(1)]';
  J=X'*D;                                    % jacobian matrix
  detJ=J(1,1)*J(2,2)-J(1,2)*J(2,1);          % jacobian
  invJ=1/detJ*[ J(2,2) -J(1,2); ...          % inverse jacobian matrix
               -J(2,1)  J(1,1)];
  G=D*invJ;                                  % gradient of shape functions
  B=[G(1,1) 0 G(2,1) 0 G(3,1) 0 ...
     G(4,1) 0 G(5,1) 0 G(6,1) 0;
     0 G(1,2) 0 G(2,2) 0 G(3,2)...
     0 G(4,2) 0 G(5,2) 0 G(6,2);
     G(1,2) G(1,1) G(2,2) G(2,1)...
     G(3,2) G(3,1) G(4,2) G(4,1)...
     G(5,2) G(5,1) G(6,2) G(6,1)];
  sigma=(A*B*Ue)';
  GU=[G(1,1) 0 G(2,1) 0 G(3,1) 0 ...
      G(4,1) 0 G(5,1) 0 G(6,1) 0;
      G(1,2) 0 G(2,2) 0 G(3,2) 0 ...
      G(4,2) 0 G(5,2) 0 G(6,2) 0;
      0 G(1,1) 0 G(2,1) 0 G(3,1)...
      0 G(4,1) 0 G(5,1) 0 G(6,1);
      0 G(1,2) 0 G(2,2) 0 G(3,2)...
      0 G(4,2) 0 G(5,2) 0 G(6,2)];
  G_d=(GU*Ue);                                 % displacement gradient
  G_t=(GU*Ge);                                 % theta gradient
  prod_G=[G_d(1)*G_t(1)+G_d(2)*G_t(3) ...      % contraction of the two gradients
          G_d(1)*G_t(2)+G_d(2)*G_t(4) ...
          G_d(3)*G_t(1)+G_d(4)*G_t(3) ...
          G_d(3)*G_t(2)+G_d(4)*G_t(4)];
  div_t=G_t(1)+G_t(4);
  eps=[G_d(1) G_d(4) G_d(2)+G_d(3)];           % strain in engineering notation
  symprod_G=[prod_G(1) prod_G(4) ...           % symmetric part of prod_G 
             prod_G(2)+prod_G(3)];
  GBe=GBe+(symprod_G*sigma'- ...
       .5*eps*sigma'*div_t)*detJ*w_gauss(g);
end

