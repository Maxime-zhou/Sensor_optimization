%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Linear elastic stiffness matrix: T6 axi
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Ke=T6_AX_solid_Ke(X,mate)

E=mate(1);
nu=mate(2);
A=E/((1+nu)*(1-2*nu))*[1-nu,nu,nu,0;
                       nu,1-nu,nu,0;
                       nu,nu,1-nu,0;
                       0,0,0,(1-2*nu)/2]; 

a1=.3333333333333;                           % seven points Gauss rule
a2=.4701420641051;
a3=.1012865073235;
a_gauss=[a1 a1 a1; 
         a2 a2 1-2*a2; a2 1-2*a2 a2; 1-2*a2 a2 a2;
         a3 a3 1-2*a3; a3 1-2*a3 a3; 1-2*a3 a3 a3];       
w1=.1125;                   
w2=0.066197076394253;
w3=0.062969590272414;
w_gauss=[w1 w2 w2 w2 w3 w3 w3];              % Gauss weights 
                  
%a_gauss=1/6*[4 1 1; 1 4 1; 1 1 4];           % three points Gauss abscissae 
%w_gauss=[1/6 1/6 1/6];                       % three points Gauss weights 

Ke=zeros(12,12);

for g=1:length(w_gauss),                     % loop over Gauss points
 a=a_gauss(g,:);                             % param. coordinates for gauss point
 D=[4*a(1)-1 0 -4*a(3)+1 4*a(2)...           % derivative of shape functions...
     -4*a(2) 4*(a(3)-a(1));                  % w.r.t. a_1,a_2
     0 4*a(2)-1 -4*a(3)+1 4*a(1) ...
     4*(a(3)-a(2)) -4*a(1)]';
 J=X'*D;                                     % jacobian matrix
 detJ=J(1,1)*J(2,2)-J(1,2)*J(2,1);           % jacobian
 invJ=1/detJ*[ J(2,2) -J(1,2); ...           % inverse jacobian matrix
              -J(2,1)  J(1,1)];
 G=D*invJ;                                   % gradient of shape functions  
 N=[a(1)*(2*a(1)-1) a(2)*(2*a(2)-1) ...      % list of shape functions
   a(3)*(2*a(3)-1) 4*a(1)*a(2) ...
   4*a(2)*a(3) 4*a(1)*a(3)];
 r=N*X(:,1);                                 % radial coordinate
 B=[G(1,1) 0 G(2,1) 0 G(3,1) 0 ...           % first line as in T6
    G(4,1) 0 G(5,1) 0 G(6,1) 0;         
    0 G(1,2) 0 G(2,2) 0 G(3,2) ...           % second line as in T6
    0 G(4,2) 0 G(5,2) 0 G(6,2);
    N(1)/r 0 N(2)/r 0 N(3)/r 0 ...           % new line for eps_{tt}
    N(4)/r 0 N(5)/r 0 N(6)/r 0;
    G(1,2) G(1,1) G(2,2) G(2,1) G(3,2) G(3,1)... % last line as in T6
    G(4,2) G(4,1) G(5,2) G(5,1) G(6,2) G(6,1)];
 Ke=Ke+B'*A*B*detJ*w_gauss(g)*r;             % contrib to stiffness matrix
end
