
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Stresses in one T3 element due to displacements Ve (constant): linear
% elastic constitutive law
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Sg=T3_2S_solid_Sg(X,mate,Ue)

E=mate(1);
nu=mate(2);
A=E/(1-nu^2)*[1,nu,0;
              nu,1,0;
              0,0,(1-nu)/2];   

x11=X(1,1); x21=X(2,1); x31=X(3,1);          % nodal coordinates
x12=X(1,2); x22=X(2,2); x32=X(3,2);
S=.5*((x21-x11)*(x32-x12)-...
      (x31-x11)*(x22-x12));                  % element area
B=[x22-x32,0,x32-x12,0,x12-x22,0;
   0,x31-x21,0,x11-x31,0,x21-x11;
   x31-x21,x22-x32,x11-x31, ...
   x32-x12,x21-x11,x12-x22]/(2*S);
Sg=(A*B*Ue)';
