
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% semi-infinite bar with imposed traction and direct integration: Absorbing Boundary
% Conditions
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

addpath(genpath('.'))

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L=1;                                         % bar length
NE=20;                                       % number of elements
rho=1.;                                      % density
E=1.;                                        % Young modulus
TD=1;                                        % value of traction x=L
tf=20;                                       % final time
tinv=5;
omega=1;
%Dt=.01445;                                  % time step
Dt=.001;                                     % time step
beta=.0;                                     % parameters for Newmark algorithm...
gamma=.5;                                    % beta=0; gamma=.5;  central differences...
                                             % beta=.25; gamma=.5;  incond. stable                                           
c=sqrt(E/rho);
                                             
                                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pre-processor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=L/NE;                                      % size of one element
NN=NE+1;                                     % number of nodes
coor=[0:h:L]';                               % nodal coordinates: uniform discretization 
neq=NN;                                    % number of equations
dof=[1:1:NN];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assemblage of matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M=zeros(neq,neq);                            % allocates mass matrix
K=zeros(neq,neq);                            % allocates rigidity matrix
F=zeros(neq,1);                              % allocates rhs side

for e=1:NE,                                 % loop over elements to assemble matrices
 X=[coor(e) coor(e+1)];                     % creates segment
 dofe=dof(e:e+1);
 pe=find(dofe>0);
 Ie=dofe(pe);                               % gets value of associated DOFs 
 Me=B2_1D_wave_Me(X,rho);                   % element mass consistent matrix
 M(Ie,Ie)=M(Ie,Ie)+Me(pe,pe);               % assemblage of mass matrix
 Ke=B2_1D_wave_Ke(X,E);                     % element stiffness matrix
 K(Ie,Ie)=K(Ie,Ie)+Ke(pe,pe);               % assemblage of stiffness matrix
end
F(neq)=TD;                                  % rhs 
M(1,1)=M(1,1)+E/c*gamma*Dt;                 % contrib from ABC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial conditions and acceleration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

omj2=eigs(K,M,1);                            % computes stability interval for...
stable_time_increment=2/sqrt(omj2)           % explicit option 

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
  U=U+Dt*Ud+0.5*Dt^2*(1-2*beta)*Udd;         % prediction phase for displacements
  Ud=Ud+Dt*(1-gamma)*Udd;                    % prediction phase for velocities
  F(1)=-E/c*Ud(1);
  Fg=F-K*U;      
  Udd=CS\(CS'\Fg);                           % new accelerations 
  U=U+Dt^2*beta*Udd;                         % actualisation of displacements
  Ud=Ud+Dt*gamma*Udd;                        % actualisation of velocities
  history(step,3)=Dt*step;                   % saves history at two points for post-proc.
  history(step,1)=U(neq);                    % tip point 
  history(step,2)=U(1);                      % x=0 point
end

clear CS K

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bkgcol=.93*[1 1 1];
figure('Color',bkgcol)                       % plots numerical and exact radial displ.
xlabel('Time')
ylabel('Displacement')
plot(history(:,3),history(:,1),'k-') 
hold on
plot(history(:,3),history(:,2),'k--') 
axis([0 tf min(history(:,1)) ...
	                max(history(:,1))]);
xlabel('Time','FontSize',12)
ylabel('Displacement','FontSize',12)
legend('Tip point','Base point') 
set(gca,'Color','w')                         % sets color of current axes

