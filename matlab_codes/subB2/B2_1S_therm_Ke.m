
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% computes element conductivity matrix for B2 element with spherical symmetry
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Ke=B2_1S_therm_Ke(X,k)

r1=X(1); r2=X(2);                            % node coordinates
h=r2-r1;                                     % element length 
Ke=k/(3*h)*(r1^2+r1*r2+r2^2)*[1,-1; -1,1];   % element conductivity matrix

