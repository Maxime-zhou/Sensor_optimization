
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Radial Return for Von Mises Linear Isotropic Hardening (Plane strain analysis)
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function [Dp,sigma_hat]=RR_VonMises_2A_R(mate,sigma,p,Deps)

E=mate(1);                                   % Young modulus
nu=mate(2);                                  % Poisson coefficient
sigma0=mate(3);                              % Yield limit
h=mate(4);                                   % hardenning coefficient
mu=E/(2*(1+nu));                             % Lame coeff
kappa=E/(3*(1-2*nu));                        % bulk modulus

vl=[1 1 1 0]';                               % identity tensor
trDeps=Deps(1)+Deps(2);                      % trace of increment of strain
De=[Deps(1:2); 0; Deps(3)/2]-1/3*trDeps*vl;  % Deviatoric tensor (4 components)

sigma_elas=sigma+kappa*trDeps*vl+2*mu*De;    % elastic prediction of stresses
trsigma=sum(sigma_elas(1:3));                % volumetric part (stress)
s_elas=sigma_elas-1/3*trsigma*vl;            % deviatoric elastic prediction (stress)
sigeq_elas=sqrt(1.5*...                      % equivalent deviatoric stress 
            (s_elas(1)^2+s_elas(2)^2+ ...
             s_elas(3)^2+2*s_elas(4)^2));    

f_elas=sigeq_elas-h*p-sigma0;
if(f_elas>0)                                 % if plastic process
  Dp=f_elas/(3*mu+h);                        % increment of plastic eq. strain
  Depsp=3/2*Dp*s_elas/sigeq_elas;            % increment over step of plastic strain
  sigma_hat=sigma_elas-2*mu*Depsp;           % new total stress
else                                         % elseif elastic process
  Dp=0;   
  sigma_hat=sigma_elas;                      % new total stress
end


