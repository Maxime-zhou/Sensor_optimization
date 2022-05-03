
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Radial Return for Von Mises Linear Isotropic Hardening (Plane strain analysis)
%               with tangent matrix
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function [AEP,Dp,sigma_hat]=RRTM_VonMises_2A_R(mate,sigma,p,Deps)

E=mate(1);                                 % Young modulus
nu=mate(2);                                % Poisson coefficient
sigma0=mate(3);                            % Yield limit
h=mate(4);                                 % hardenning coefficient
mu=E/(2*(1+nu));                           % Lame coeff
kappa=E/(3*(1-2*nu));                      % bulk modulus
A  = E/((1+nu)*(1-2*nu))*[1-nu,nu,0;       % elastic tensor (plane strain)
                          nu,1-nu,0;
                          0,0,(1-2*nu)/2];  
MK=1/3*[2 -1 -1 0; -1 2 -1 0; ...          % deviatoric extractor
        -1 -1 2 0; 0  0  0 3/2];
vl=[1 1 1 0]';                             % identity tensor
trDeps=Deps(1)+Deps(2);                    % trace of increment of strain
De=[Deps(1:2); 0; Deps(3)/2]-1/3*trDeps*vl;% Deviatoric tensor (4 comp.) 

sigma_elas=sigma+kappa*trDeps*vl+2*mu*De ; % elastic prediction of stresses 
trsigma=sum(sigma_elas(1:3));              % volumetric part (stress)
s_elas=sigma_elas-1/3*trsigma*vl;          % Deviatoric elastic prediction
sigeq_elas=sqrt(1.5*...                    % equivalent deviatoric stress
            (s_elas(1)^2+s_elas(2)^2+ ...
             s_elas(3)^2+2*s_elas(4)^2));     

f_elas=sigeq_elas-h*p-sigma0;
if(f_elas>0)                               % if plastic process
  Dp=f_elas/(3*mu+h);                      % incr acc. plastic strain
  n_elas=s_elas/sigeq_elas;
  Depsp=3/2*Dp*n_elas;                     % incr. of plastic strain
  sigma_hat=sigma_elas-2*mu*Depsp;         % new total stress

  beta=3*mu*Dp/sigeq_elas;                 % coefficients gamma and beta
  gamma=3*mu/(3*mu+h);
  D=3*mu*(gamma-beta)*n_elas*n_elas'+...   % D matrix
    2*mu*beta*MK;
  AEP=A-[D(1:2,1:2) D(1:2,4); ...          % % A-D (select components)
           D(4,1:2) D(4,4)]; 
else                                       % elseif elastic process
  Dp=0;   
  sigma_hat=sigma_elas;                    % new total stress
  AEP=A;
end
