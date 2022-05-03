
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Stresses in one T3 element due to displacements Ve (constant): linear
% elastic constitutive law
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Sg=T3_2D_therm_Sg(T,mate,Ue)

k=mate(1);
x11=T(1,1); x21=T(2,1); x31=T(3,1);          % nodal coordinates
x12=T(1,2); x22=T(2,2); x32=T(3,2);
S=.5*((x21-x11)*(x32-x12)-...
      (x31-x11)*(x22-x12));                  % element area
G=[x22-x32,x32-x12,x12-x22;
    x31-x21,x11-x31,x21-x11]/(2*S);
Sg=(k*G*Ue)';
