
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Elastic sphere under imposed pressure and displacements
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

addpath(genpath('.'))

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R1=1; R2=2;                                  % inner and outer radii
NE=5;                                        % number of elements
E=1.;                                        % Young modulus
nu=.3;                                       % Poisson coefficient
p=1;                                         % value of inner pressure
d=-.1;                                       % value of imposed displ at r=r2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pre-processor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=(R2-R1)/NE;                                % size of one element
NN=NE+1;                                     % number of nodes
coor=[R1:h:R2]';                             % nodal coordinates: uniform discretization 
mu=E/(2*(1+nu));                             % definition of Lame constants
lambda=2*nu*mu/(1-2*nu);

neq=NN-1;                                    % number of unknown nodal values
U=zeros(NN,1);                               % initialisation of displacements
U(NN)=d;                                     % puts imposed displacement on last node

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assemblage of stiffness matrix and rhs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K=zeros(neq,neq);                            % allocates stiffness matrix
F=zeros(neq,1);                              % allocates rhs side
for e=1:NE                                   % assemblage of stiffness matrix
 X=[coor(e) coor(e+1)];                      % creates segment
 Ke=B2_1S_solid_Ke(X,lambda,mu);             % element stiffness matrix
 Ie=[e e+1];                                 % sets nodal degrees of freedom
 if e<NE                                     % if not last element 
  K(Ie,Ie)=K(Ie,Ie)+Ke;                      % the whole Ke goes into K 
 else
  K(neq,neq)=K(neq,neq)+Ke(1,1);             % else only the first coefficient
  F(neq)=-Ke(1,2)*U(NN);                     % and Ke also contributes to rhs
 end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assemblage of rhs due to external loads
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F(1)=F(1)+coor(1)^2*p;                       % contribution to rhs from inner pressure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solution phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U(1:NN-1)=K\F;                               % solution of linear system

clear K F                                    % clears memory

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% post phase: computes radial stresses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stress=zeros(NE,1);                          % vector for saving radial stresses 
for e=1:NE
 X=[coor(e) coor(e+1)];                      % creates segment
 Ue=U(e:e+1);                                % gets nodal displacements
 stress(e)=B2_1S_solid_Sg(X,lambda,mu,Ue);   % stresses in the middle of element
end
centers=.5*(coor(1:NN-1)+coor(2:NN));        % vector with center-points of elements

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% post phase: evaluates exact solution on a finer grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coor2=[R1:(R2-R1)/100:R2];                   % finer grid
a=(-(1+nu)*p*R1^3+2*d*E*R2^2)/...            
         ((1+nu)*R1^3+2*(1-2*nu)*R2^3);
b=R1^3*R2^2*(d*E+(1-2*nu)*p*R2)/...
         ((1+nu)*R1^3+2*(1-2*nu)*R2^3);
pa=a*(1-2*nu)/E;
pb=b*(1+nu)/E;
U_ex=pa.*coor2+pb./coor2.^2;                 % exact displacements 
stress_ex=a-2*b./(coor2.^3);                 % exact radial stresses

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% post phase: plots comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bkgcol=1.*[1 1 1];
figure('Color',bkgcol)                       % plots numerical and exact radial displ.
plot(coor,U,'ko','MarkerSize',8)
%plot(coor,U,'k:')
hold on
plot(coor2,U_ex,'k')
xlabel('Radial coordinate','FontSize',12)
ylabel('Displacement u_r','FontSize',12)
legend('FEM','exact') 
set(gca,'Color','w')                         % sets color of current axes

figure('Color',bkgcol)                       % plots numerical and exact sigma_rr
plot(centers,stress,'ko','MarkerSize',8)
hold on
plot(coor2,stress_ex,'k')
xlabel('Radial coordinate','FontSize',12)
ylabel('Stress \sigma_{rr}','FontSize',12)
legend('FEM','exact','Location','Southeast') 
set(gca,'Color','w')                         % sets color of current axes

