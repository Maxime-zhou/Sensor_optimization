
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Heat diffusion: sphere with imposed temperature at external boundary
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

addpath(genpath('.'))

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R=1;                                         % outer radius
NE=20;                                       % number of elements
rhoc=10.;                                    % capacity
k=1.;                                        % conductivity
TD=1;                                        % value of temperature in r=R
tf=5;                                        % final time
Dt=.001;                                     % time step
implicit=0;                                  % implicit or not? 
if implicit==1, lumped=0;            
else lumped=1;                               % if explicit lumped cap matrix
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pre-processor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=R/NE;                                      % size of one element
NN=NE+1;                                     % number of nodes
coor=[0:h:R]';                               % nodal coordinates: uniform discretization 
neq=NN-1;                                    % number of equations

U=zeros(NN,1);                               % initialisation of temperatures
U(NN)=TD;                                    % boundary conditions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assemblage of matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K=zeros(neq,neq);                            % allocates conductivity matrix
M=zeros(neq,neq);                            % allocates capacity matrix
F=zeros(neq,1);                              % allocates rhs side for diffusion part

for e=1:NE,                                  % loop over elements to assemble matrices
 Xe=[coor(e) coor(e+1)];                     % creates segment
 Me=B2_1S_therm_Me(Xe,rhoc,lumped);          % element capacity 
 Ke=B2_1S_therm_Ke(Xe,k);                    % element conductivity
 Ie=[e e+1];                                 % sets nodal degrees of freedom
 if e<NE                                     % if not last element 
  M(Ie,Ie)=M(Ie,Ie)+Me;                      % the whole Me goes into M 
  K(Ie,Ie)=K(Ie,Ie)+Ke;                      % the whole Ke goes into K 
 else
  M(neq,neq)=M(neq,neq)+Me(1,1);             % else only the first coefficient
  K(neq,neq)=K(neq,neq)+Ke(1,1);             % else only the first coefficient
  F(neq)=-Ke(1,2)*U(NN);                     % and Ce also contributes to rhs
 end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time marching solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if implicit==0,                              % computes stable  Dt for explicit option
 uu=eigs(K,M,1);
 stable_Dt=2/uu
end  
%Dt=.99*stable_Dt;

if implicit==1,                              % prepares Cholesky decomposition
 K=K+M/Dt;                                 
 KC=chol(K);
 clear K
end

T=zeros(neq,1);
nstep=floor(tf/Dt);                          % number of time steps assuming constant dt
vec_time=zeros(nstep,1);                     % vector for time values (for post-process.)
vec_temperature=zeros(nstep,1);              % vector for cavity temperatures

%-----------------------------------------------------------------------------------------
disp(['Number of elements: ' num2str(NE)])
percold=0;
%-----------------------------------------------------------------------------------------

for step=1:nstep,                            % loop over time steps

%-----------------------------------------------------------------------------------------
 perc=10*floor(10*step/nstep);
 if perc~=percold,
  percold=perc;
  disp(['Time history: ' num2str(percold) '%'])
 end
%-----------------------------------------------------------------------------------------

 if implicit==1,                             
  Fg=F+1/Dt*M*T;                             % implicit case
  T=KC\(KC'\Fg);
 else
  Fg=F-K*T;                                  % explicit case
  T=T+Dt*Fg./diag(M);                        
 end 
 U(1:NN-1)=T;                                % temps over all nodes
 vec_time(step)=step*Dt;
 vec_temperature(step)=U(1);
end

clear RK M


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots solution and compares with exact one
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bkgcol=1.*[1 1 1];
figure('Color',bkgcol)                       % plots numerical and exact radial displ.
%axis square
plot(coor,U,'k')
title('Temperature along the radius at final time')
xlabel('Radial coordinate')
ylabel('Temperature')
set(gca,'Color','w')                         % sets color of current axes

bkgcol=1.*[1 1 1];
figure('Color',bkgcol)                       % plots numerical and exact radial displ.
%axis square
axis on
plot(vec_time,vec_temperature,'k-.')
lambda=k/rhoc;
nn=100;
for step=1:nstep,
  time=vec_time(step);
  an=(1:nn);
  vn=exp(-((an*pi/R).^2)*lambda*time)*2;     % The exact solution is an infinite series...
  an=(-1).^an;                               % here truncated to nn terms 
  vec2(step,1)=an*vn'+1;
end
hold on
plot(vec_time,vec2,'k')
title('Temperature history at sphere center')
xlabel('Time','FontSize',12)
ylabel('Temperature','FontSize',12)
legend('FEM','exact','Location','Southeast') 
set(gca,'Color','w')                         % sets color of current axes
%axis([0 tf 0 1.4])
axis([min(vec_time) max(vec_time) ...
      min(vec2) 1.05*max(vec2)]);


