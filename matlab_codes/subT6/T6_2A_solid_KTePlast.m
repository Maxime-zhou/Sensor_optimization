
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Local plastic correction for T6
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function [KEPe, Finte, sigma_n, Dp]=T6_2A_solid_KTePlast(T, mate, DUe, sigma, p)

E  = mate(1);
nu = mate(2);
A  = E/((1+nu)*(1-2*nu))*[1-nu,nu,0;         % Plane-strain linear
                          nu,1-nu,0;         % isotropic elastic matrix
                          0,0,(1-2*nu)/2];   %

a_gauss=1/6*[4 1 1; 1 4 1; 1 1 4];           % Gauss abscissae
w_gauss=[1/6 1/6 1/6];                       % Gauss weights

KEPe=zeros(12,12);
Finte=zeros(12,1);
for g=1:3,                                   % loop over Gauss points
  a=a_gauss(g,:);                            % param. coordinates for gauss point
  DN=[4*a(1)-1 0 -4*a(3)+1 4*a(2)...         % derivative of shape functions...
      -4*a(2) 4*(a(3)-a(1));                 % w.r.t. a_1,a_2
      0 4*a(2)-1 -4*a(3)+1 4*a(1) ...
      4*(a(3)-a(2)) -4*a(1)]';
  J=T'*DN;                                   % jacobian matrix
  detJ=J(1,1)*J(2,2)-J(1,2)*J(2,1);          % jacobian
  invJ=1/detJ*[ J(2,2) -J(1,2); ...          % inverse jacobian matrix
               -J(2,1)  J(1,1)];
  GN=DN*invJ;                                % gradient of shape functions
  Be=[GN(1,1) 0 GN(2,1) 0 GN(3,1) 0 ...
      GN(4,1) 0 GN(5,1) 0 GN(6,1) 0;
      0 GN(1,2) 0 GN(2,2) 0 GN(3,2)...
      0 GN(4,2) 0 GN(5,2) 0 GN(6,2);
      GN(1,2) GN(1,1) GN(2,2) GN(2,1)...
      GN(3,2) GN(3,1) GN(4,2) GN(4,1)...
      GN(5,2) GN(5,1) GN(6,2) GN(6,1)];
  Deps=Be*DUe;                               % computes elastic increment of stresses
  [AEP,Dp(g,:),sigma_n(g,:)]=...
    RRTM_VonMises_2A_R(mate,...              % radial return algorithm
                    sigma(g,1:4)',p(g),Deps);
  KEPe  = KEPe+Be'*AEP*Be*detJ*w_gauss(g);   % tangent matrix
  Finte = Finte - ...                        % internal nodal forces
    Be'*sigma_n(g,[1:2, 4])'*detJ*w_gauss(g);
end

