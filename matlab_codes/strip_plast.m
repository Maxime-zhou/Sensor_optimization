
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% FILE: strip_plast.m
% Plane strain strip under uniform tension
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

addpath(genpath('.'))

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L=1;                                         % halflength of strip
E=1000;                                       % Young modulus
nu=.1;                                       % Poisson coefficient
sigma0=1;                                    % yield limit     
H=0;                                         % hardening parameter
mate= [E nu sigma0 H];                       % material parameters: [E, nu, sigma0, H]

%q=[0:.005:.12];                                % history of imposed displacements
q0=2*sigma0*L*(1-nu^2)/(E*sqrt(1-nu+nu^2))   % elastic limit
step=q0/10;
q=[0,q0:step:6*q0];                           % history of imposed displacements
numstep=length(q);                           % number of steps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre-processor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p=0;                                         % ini. c. plastic strains
sigma=zeros(4,1);                            % ini. stresses
toll=1.d-4;                                  % user fixed tolerance
output=zeros(numstep,2);                     % save solution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% History analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for istep=2:numstep,                         % loop over all load steps
  Dq=q(istep)-q(istep-1);                    % delta of displacement applied in step
  iter=0;
  resid=1;
  Dp=0;                                      % init. of plastic mult. increment in step
  Deps=zeros(3,1);

  while resid > toll,                        % solution of one load step
    iter=iter+1;   
    if iter==1                               % imposes equilibrium: sigma_22=0 
      Deps=Dq/(2*L)*[1 -nu/(1-nu) 0]';       % increment of total eps_11
    else 
      deps2=-sigma_hat(2)*(1+nu)* ...
                       (1-2*nu)/(E*(1-nu));           
      Deps=Deps+[0 deps2 0]';
    end
    [Dp,sigma_hat]=RR_VonMises_2A_R(mate,sigma,p,Deps);
    resid=abs(sigma_hat(2));                 % abs value of sig_22 taken as residuum
  end                                        % end of one load step

  p=p+Dp;                                    % increments plastic multiplier
  sigma=sigma_hat;
  output(istep,:)=[sigma(1) sigma(3)];       % save solution
end                                          % end loop over time steps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post-processing: comparison with exact results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

history=q/(2*L);
thetae=asin(sqrt(3)/(2*sqrt(1-nu+nu^2)))-pi/6;
sigmae=sigma0/sqrt(1-nu+nu^2);
epse=(1-nu^2)/E*sigmae;
const=epse-2*sigma0/(sqrt(3)*E)*...
	     (3/4*log(abs(tan(thetae/2+pi/3)))+...
      (1-2*nu)*sin(thetae-pi/6));
vtheta=[thetae:.005:pi/3];
vsigmaxx=2*sigma0/sqrt(3)*sin(vtheta+pi/6);
vsigmazz=2*sigma0/sqrt(3)*sin(vtheta-pi/6);
veps=const+2*sigma0/(sqrt(3)*E)*...
	      (3/4*log(abs(tan(vtheta/2+pi/3)))+...
       (1-2*nu)*sin(vtheta-pi/6));
veps=[0 veps];
vsigmaxx=[0 vsigmaxx];
vsigmazz=[0 vsigmazz];

figure('Color','w')
plot(veps,vsigmaxx,'k-')
hold on
plot(history,output(:,1),'ko')
hold on
plot(veps,vsigmazz,'k-.')
hold on
plot(history,output(:,2),'k*')
axis([min(veps),1.1*max(veps),...
	     min(vsigmazz),1.1*max(vsigmaxx)])

xlabel('\epsilon_{11}','Fontsize',12)
ylabel('stress','Fontsize',12)
legend('Exact \sigma_{11}',...
       'Numerical \sigma_{11} (50 load steps)',...
	      'Exact \sigma_{33}', ...
	      'Numerical \sigma_{33} (50 load steps)',...
	      'Location','SouthEast') 
      
      