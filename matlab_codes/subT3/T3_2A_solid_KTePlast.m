
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%   Tangent plastic matrix: T3 (element level)
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function [KEPe, Finte, sigma_hat, Dp]=T3_2A_solid_KTePlast(T, mate, DUe, sigma, p)

E  = mate(1);
nu = mate(2);
A=E/((1+nu)*(1-2*nu))*[1-nu,nu,0;       % Plane-strain linear
                       nu,1-nu,0;       % isotropic elastic matrix
                       0,0,(1-2*nu)/2];  
   
x11=T(1,1); x21=T(2,1); x31=T(3,1);     % nodal coordinates
x12=T(1,2); x22=T(2,2); x32=T(3,2);
S=.5*((x21-x11)*(x32-x12)-...
      (x31-x11)*(x22-x12));             % element area
Be=[x22-x32,0,x32-x12,0,x12-x22,0;
    0,x31-x21,0,x11-x31,0,x21-x11;
    x31-x21,x22-x32,x11-x31, ...
    x32-x12,x21-x11,x12-x22]/(2*S);
Deps=Be*DUe;                            % computes increment of strain
[AEP,Dp,sigma_hat(1,:)]=...             % radial return algorithm
    RRTM_VonMises_2A_R... 
      (mate,sigma(1,1:4)',p,Deps);
KEPe  =  S*Be'*AEP*Be;                  % element stiffness matrix

%[Dp,sigma_hat(1,:)]=RR_VonMises2A...    % radial return algorithm
%      (mate,sigma(1,1:4)',p,Deps);
%KEPe  =  S*Be'*A*Be;                    % element stiffness matrix

Finte = -S*Be'*sigma_hat([1:2 4])';     % nodal forces equivalent to stresses
