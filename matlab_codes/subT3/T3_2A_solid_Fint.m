%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%  Internal nodal forces, large displacements: T3
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function Fe=T3_2A_solid_Fint(X,mate,Ue)

E=mate(1);
nu=mate(2);
A=E/(1-nu^2)*[1,nu,0;
              nu,1,0;
              0,0,(1-nu)/2];
          
x11=X(1,1); x21=X(2,1); x31=X(3,1);          % nodal coordinates 
x12=X(1,2); x22=X(2,2); x32=X(3,2);
S=((x21-x11)*(x32-x12)-...
   (x31-x11)*(x22-x12))/2;
GV=1/(2*S)*[x22-x32,0,x32-x12,0,x12-x22,0;    % gradient matrix 
            x31-x21,0,x11-x31,0,x21-x11,0;
            0,x22-x32,0,x32-x12,0,x12-x22;
            0,x31-x21,0,x11-x31,0,x21-x11];
grad=GV*Ue;                                   % gradient of u^{(k)}
e=[grad(1)+.5*grad(1)^2+.5*grad(3)^2, ...    % Green-Lagrange tensor e^{(k)} 
   grad(4)+.5*grad(2)^2+.5*grad(4)^2, ...
   grad(2)+grad(3)+grad(1)*grad(2)+ ...     
   grad(3)*grad(4)]';
PK=A*e;                                      % Piola Kirkhhoff tensor pi^{(k)}
B=[GV(1,:); GV(4,:); GV(2,:)+GV(3,:)];
Bu=[grad(1)*GV(1,:)+grad(3)*GV(3,:);
    grad(2)*GV(2,:)+grad(4)*GV(4,:);
    grad(1)*GV(2,:)+grad(3)*GV(4,:)+ ...
    grad(2)*GV(1,:)+grad(4)*GV(3,:)];
Bt=B+Bu;
Fe=Bt'*PK*S;

