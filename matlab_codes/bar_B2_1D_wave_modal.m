
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Bar with imposed traction: modal analysis and static correction
% Heaviside o blast loading
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
tf=10;                                       % final time
Dt=.01;

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

bN=4;
[V,D]=eigs(K,M,bN,'SM');                     % modal base

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time marching solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f=V'*F;
resid=F-M*V*f;
Du=K\resid;                                  % static correction when g(t)=1

nstep=floor(tf/Dt);                          % number of time stepsdt
history=zeros(nstep,3);
history(:,3)=[Dt:Dt:tf];

D=diag(D);
omega=sqrt(D);
U=zeros(bN,1);

for step=1:nstep,                            % loop over time steps
 time=Dt*step; 
 if blast==1,                                % if impulse
  U=-f./(T^2*omega.^3).*(sin(omega*time)- ... 
              2*sin(omega*(time-T))+ ...     % hat approx of impulse
              sin(omega*(time-2*T))); 
  if time > 2*T,
   fun=0;                                    % after impulse
  elseif time > T, 
   fun=(1/T)^2*(2*T-time);                   % descending ramp
  else
   fun=(1/T)^2*time;                         % ascending ramp
  end
 else                                        % heaviside loading
  U=f./D.*(1-cos(omega*time));   
  fun=1;                                     % multiplier for static correction
 end
 history(step,3)=Dt*step;                    % saves history for post-proc.
 history(step,1)=V(neq,:)*U+fun*Du(neq);
 history(step,2)=V(round(neq/2),:)*U+fun*Du(round(neq/2));
end

clear CS K

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bkgcol=1*[1 1 1];
figure('Color',bkgcol)                       % plots numerical and exact radial displ.
xlabel('Time')
ylabel('Displacement')
plot(history(:,3),history(:,1),'k-') 
hold on
plot(history(:,3),history(:,2),'k-.') 
if blast==1
axis([0 tf min(history(:,1)) ...
	                max(history(:,1))]);
else
 axis([0 tf -.1 2])
end    
xlabel('Time','FontSize',12)
ylabel('Displacement','FontSize',12)
legend('Tip point','Mid point') 
set(gca,'Color','w')                         % sets color of current axes

