
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Bar with imposed traction: direct integration with Newmark scheme
% Heaviside o blast loading
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

addpath(genpath('.'))

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L=1;                                         % bar length
NE=100;                                      % number of elements
rho=1.;                                      % density
E=1.;                                        % Young modulus
TD=1;                                        % value of traction x=L
tf=10;                                       % final time
Dt=.01;                                      % time step
gamma=.55;                                   % parameters for Newmark algorithm...
beta=.25;                                    % beta=0; gamma=.5;  central differences...
beta=(gamma/2+1/4)^2;                        % highest diss for high freq for givn gamma

blast=0;                                     % blast loading o heaviside?
T=.1;                                        % ramp duration for blast loading

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pre-processor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=L/NE;                                      % size of one element
NN=NE+1;                                     % number of nodes
coor=[0:h:L]';                               % nodal coordinates: uniform discretization 
neq=NN-1;                                    % number of equations
dof=[0 1:1:NN-1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assemblage of matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M=zeros(neq,neq);                            % allocates mass matrix
K=zeros(neq,neq);                            % allocates rigidity matrix
F=zeros(neq,1);                              % allocates rhs side

for e=1:NE,                                  % loop over elements to assemble matrices
 X=[coor(e) coor(e+1)];                      % creates segment
 dofe=dof(e:e+1);
 pe=find(dofe>0);
 Ie=dofe(pe);                                % gets value of associated DOFs 
 Me=B2_1D_wave_Me(X,rho);                    % element mass consistent matrix
 M(Ie,Ie)=M(Ie,Ie)+Me(pe,pe);                % assemblage of mass matrix
 Ke=B2_1D_wave_Ke(X,E);                      % element stiffness matrix
 K(Ie,Ie)=K(Ie,Ie)+Ke(pe,pe);                % assemblage of stiffness matrix
end
F(neq)=F(neq)+TD;                            % rhs 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial conditions and acceleration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

omj2=eigs(K,M,1);                            % computes stability interval for...
stable_time_increment=2/sqrt(omj2);          % explicit option 

U=zeros(neq,1);                              % initial value of displacements
Ud=zeros(neq,1);                             % initial value of velocities
Udd=M\F;                                     % initial value of accelerations
S=M+beta*Dt^2*K;                          
CS=chol(S);                                  % Cholesky decomposition of S matrix 
clear S M  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time marching solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nstep=floor(tf/Dt);                          % number of time steps with constant dt
history=zeros(nstep,3);

time=0;
for step=1:nstep,                            % loop over time steps
 time=time+Dt;
 U=U+Dt*Ud+0.5*Dt^2*(1-2*beta)*Udd;          % prediction phase for displacements
 Ud=Ud+Dt*(1-gamma)*Udd;                     % prediction phase for velocities
 
 % loading definition
 fun=1;                                      % defines loading for heaviside
 if blast==1                                 % or for blast 
  if time < T,
   fun=(1/T)^2*time;                         % ascending ramp
  elseif time < 2*T, 
   fun=(1/T)^2*(2*T-time);                   % descending ramp
  else
   fun=0;                                    % zero
  end
 end
 
 % resumes Newmark scheme
 Fg=fun*F-K*U;      
 Udd=CS\(CS'\Fg);                            % new accelerations 
 U=U+Dt^2*beta*Udd;                          % actualisation of displacements
 Ud=Ud+Dt*gamma*Udd;                         % actualisation of velocities
 history(step,3)=Dt*step;                    % saves hisotry at two points for post-proc.
 history(step,1)=U(neq);
 history(step,2)=U(round(neq/2));
end

clear CS K

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bkgcol=1.*[1 1 1];
figure('Color',bkgcol)                       % plots numerical and exact radial displ.
xlabel('Time')
ylabel('Displacement')
plot(history(:,3),history(:,1),'k-') 
hold on
plot(history(:,3),history(:,2),'k-.') 

axis([0 tf 0 2])                             % heaviside loading
if blast==1
 axis([0 tf -1.2 1.2])                       % blast loading
end 

xlabel('Time','FontSize',12)
ylabel('Displacement','FontSize',12)
legend('Tip point','Mid point') 
set(gca,'Color','w')                         % sets color of current axes

