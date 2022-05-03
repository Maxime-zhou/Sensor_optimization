
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Nodal forces due to surface tractions: B3
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Fe=B3_2A_solid_Fe_surf(X,val,idir)

a_gauss=[-1/sqrt(3) 1/sqrt(3)];             % gauss-legendre two point rule
w_gauss=[1 1];
Fe=zeros(6,1);
for g=1:2                                   % loop over gauss points
 a=a_gauss(g);
 D=[a-.5  a+.5 -2*a];                      % derivative of shape functions
 J=D*X;
 detJ=sqrt(J(1)^2+J(2)^2);                  % jacobian
 TD=zeros(2,1);
 if idir>0,                                 % sets traction vector
  TD(idir)=val;
 else
  TD=val/detJ*[J(2) -J(1)]';                % traction along normal vector
 end
 NL=[.5*a*(a-1) .5*a*(1+a) 1-a^2];          % shape functions
 N=[NL(1) 0 NL(2) 0 NL(3) 0;
    0 NL(1) 0 NL(2) 0 NL(3)];
 Fe=Fe+N'*TD*detJ*w_gauss(g);
end

