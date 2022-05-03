
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Nodal forces due to surface fluxes: B3
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Fe=B3_AX_therm_Fe_surf(X,val,idir)

a_gauss=[-1/sqrt(3) 1/sqrt(3)];             % gauss-legendre two point rule
w_gauss=[1 1];
Fe=zeros(6,1);
for g=1:2                                   % loop over gauss points
 a=a_gauss(g);
 D=[a-.5  a+.5 -2*a];                       % derivative of shape functions
 J=D*X;
 detJ=sqrt(J(1)^2+J(2)^2);                  % jacobian
 N=[.5*a*(a-1) .5*a*(1+a) 1-a^2];           % shape functions
 Fe=Fe+N'*val*detJ*w_gauss(g);
end

