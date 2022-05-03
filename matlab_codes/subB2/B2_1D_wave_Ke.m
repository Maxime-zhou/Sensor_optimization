
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Element stiffness matrix for 2node line element
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Ke=B2_1D_wave_Ke(X,E)

x1=X(1); x2=X(2);
h=x2-x1;                                     % element length
Ke=E/h*[1,-1;-1,1];                          % element stiffness matrix


