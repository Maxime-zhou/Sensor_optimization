
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% T6 element : Plain strain isotropic Saint-Venant-Kirchhoff model
% Tangent stiffness matrix
% Internal force vector
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function [KTe,Finte]=T6_2A_solid_KTeNL(X,mate,Ue)

E=mate(1);
nu=mate(2);
A=E/((1+nu)*(1-2*nu))*[1-nu,nu,0;
                       nu,1-nu,0;
                       0,0,(1-2*nu)/2];  
a_gauss=1/6*[4 1 1; 1 4 1; 1 1 4];           % Gauss abscissae
w_gauss=[1/6 1/6 1/6];                       % Gauss weights
KTe=zeros(12,12);
Finte=zeros(12,1);
for g=1:3,                                   % loop over Gauss points
 a=a_gauss(g,:);                             % coordinates of gauss point
 D=[4*a(1)-1 0 -4*a(3)+1 4*a(2)...           % derivative of shape functions...
    -4*a(2) 4*(a(3)-a(1));                   % w.r.t. a_1,a_2
    0 4*a(2)-1 -4*a(3)+1 4*a(1) ...
    4*(a(3)-a(2)) -4*a(1)]';
 J=X'*D;                                     % jacobian matrix
 detJ=J(1,1)*J(2,2)-J(1,2)*J(2,1);           % jacobian
 invJ=1/detJ*[ J(2,2) -J(1,2); ...           % inverse jacobian matrix
             -J(2,1)  J(1,1)];
 G=D*invJ;                                   % gradient of shape functions
 GV=[G(1,1) 0 G(2,1) 0 G(3,1) 0 G(4,1) 0 G(5,1) 0 G(6,1) 0;
     G(1,2) 0 G(2,2) 0 G(3,2) 0 G(4,2) 0 G(5,2) 0 G(6,2) 0;
     0 G(1,1) 0 G(2,1) 0 G(3,1) 0 G(4,1) 0 G(5,1) 0 G(6,1);
     0 G(1,2) 0 G(2,2) 0 G(3,2) 0 G(4,2) 0 G(5,2) 0 G(6,2)];
 gradU=GV*Ue;                                % gradient of u
 e=[gradU(1)+.5*gradU(1)^2+.5*gradU(3)^2, ... % Green-Lagrange tensor e^{(k)} 
    gradU(4)+.5*gradU(2)^2+.5*gradU(4)^2, ...
    gradU(2)+gradU(3)+gradU(1)*gradU(2)+ ...     
    gradU(3)*gradU(4)]';
 SK=A*e;                                     % Piola tensor S^{[k]}
 B=[GV(1,:); GV(4,:); GV(2,:)+GV(3,:)];
 Bu=[gradU(1)*GV(1,:)+gradU(3)*GV(3,:);
     gradU(2)*GV(2,:)+gradU(4)*GV(4,:);
     gradU(1)*GV(2,:)+gradU(3)*GV(4,:)+ ...
     gradU(2)*GV(1,:)+gradU(4)*GV(3,:)];
 Bnl=B+Bu;
 Bg=[SK(1)*GV(1,:)+SK(3)*GV(2,:);  
     SK(3)*GV(1,:)+SK(2)*GV(2,:);
     SK(1)*GV(3,:)+SK(3)*GV(4,:);
     SK(3)*GV(3,:)+SK(2)*GV(4,:)];
 Kg=GV'*Bg;
 Finte=Finte-Bnl'*SK*detJ*w_gauss(g);              % internal force contribution
 KTe=KTe+(Bnl'*A*Bnl+Kg)*detJ*w_gauss(g);         % tangent stiff matrix contribution
end
