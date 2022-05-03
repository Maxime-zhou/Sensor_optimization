
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% T3 element : Plain strain isotropic Saint-Venant-Kirchhoff model
% Tangent stiffness matrix
% Internal force vector
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function [KTe, Finte]=T3_2A_solid_KTeNL(X,mate,Ue)

E=mate(1);
nu=mate(2);
A=E/((1+nu)*(1-2*nu))*[1-nu,nu,0;
                       nu,1-nu,0;
                       0,0,(1-2*nu)/2];  

x11=X(1,1); x21=X(2,1); x31=X(3,1);          % nodal coordinates 
x12=X(1,2); x22=X(2,2); x32=X(3,2);
S=0.5*((x21-x11)*(x32-x12)-...
       (x31-x11)*(x22-x12));                 % element area
GV=1/(2*S)*[x22-x32,0,x32-x12,0,x12-x22,0;    % gradient matrix 
            x31-x21,0,x11-x31,0,x21-x11,0;
            0,x22-x32,0,x32-x12,0,x12-x22;
            0,x31-x21,0,x11-x31,0,x21-x11];
gradU=GV*Ue;                                   % gradient of u^{(k)}
e=[gradU(1)+.5*gradU(1)^2+.5*gradU(3)^2, ...  % Green-Lagrange tensor e^{(k)} 
   gradU(4)+.5*gradU(2)^2+.5*gradU(4)^2, ...
   gradU(2)+gradU(3)+gradU(1)*gradU(2)+ ...     
   gradU(3)*gradU(4)]';
SK=A*e;                                        % Piola tensor S^{(k)}
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
Finte=-Bnl'*SK*S;
KTe=(Bnl'*A*Bnl+Kg)*S;
