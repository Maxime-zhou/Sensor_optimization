
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Nodal forces due to surface tractions: T3
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Fe=T3_2D_therm_Fe_surf(X,val,idir)

x11=X(1,1); x21=X(2,1); x31=X(3,1);          % nodal coordinates
x12=X(1,2); x22=X(2,2); x32=X(3,2);
S=.5*((x21-x11)*(x32-x12)-...
      (x31-x11)*(x22-x12));                  % element area
NL=S/3*[1 1 1]';
Fe=NL*val;                                   % nodal forces due to fluxes

end

