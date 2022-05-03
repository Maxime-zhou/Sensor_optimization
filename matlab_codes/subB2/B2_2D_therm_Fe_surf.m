
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Nodal forces due to thermal flux: B2
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Fe=B2_2D_thermo_Fe_surf(X,val,idir)

x11=X(1,1); x12=X(1,2);                      % first node
x21=X(2,1); x22=X(2,2);                      % second node
L=sqrt((x21-x11)^2+(x22-x12)^2);             % element length
NL=L/2*[1 1]';
Fe=NL*val;                                   % nodal forces due to fluxes

