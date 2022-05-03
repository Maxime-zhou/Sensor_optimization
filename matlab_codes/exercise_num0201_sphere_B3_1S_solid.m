
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Exercise 2.1
% Elastic sphere under imposed pressure and displacements: B3 elements
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

addpath(genpath('.'))

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R1=1; R2=2;                                  % inner and outer radii
NE=4;                                        % number of elements
E=1.;                                        % Young modulus
nu=0.3;                                      % Poisson coefficient
p=1;                                         % value of inner pressure
d=-0.1;                                      % value of imposed displacement at r=r2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pre-processor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=(R2-R1)/(2*NE);                            % size of one element
NN=2*NE+1;                                   % number of nodes
coor=[R1:h:R2]';                             % nodal coordinates: uniform discretization 
alpha=0.0;                                   % distorsion parameter
for n=2:2:NN-1
 coor(n)=coor(n)+alpha/2*(coor(n+1)-coor(n-1));    
end
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
 X=[coor(2*e-1) coor(2*e) coor(2*e+1)]';     % creates segment
 Ke=B3_1S_solid_Ke(X,lambda,mu);             % element stiffness matrix
 if e<NE                                     % if not last element 
  Ie=[2*e-1 2*e 2*e+1];                      % sets nodal degrees of freedom
  K(Ie,Ie)=K(Ie,Ie)+Ke;                      % the whole Ke goes into K 
 else
  Ie=[2*e-1 2*e];                            % sets nodal degrees of freedom
  Le0=[1 2];
  K(Ie,Ie)=K(Ie,Ie)+Ke(Le0,Le0);             % else only the first coefficients
  F(Ie)=-Ke(Le0,3)*U(NN);                    % and Ke also contributes to rhs
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
% post phase: evaluates exact solution on a finer grid (ONLY DISPLACEMENTS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coor2=[R1:(R2-R1)/100:R2];                   % finer grid
a=(-(1+nu)*p*R1^3+2*d*E*R2^2)/...            
         ((1+nu)*R1^3+2*(1-2*nu)*R2^3);
b=R1^3*R2^2*(d*E+(1-2*nu)*p*R2)/...
         ((1+nu)*R1^3+2*(1-2*nu)*R2^3);
pa=a*(1-2*nu)/E;
pb=b*(1+nu)/E;
U_ex=pa.*coor2+pb./coor2.^2;                 % exact displacements 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% post phase: plots comparison (ONLY DISPLACEMENTS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bkgcol=.93*[1 1 1];
figure('Color',bkgcol)                       % plots numerical and exact radial U.
plot(coor,U,'ko','MarkerSize',8)
hold on
plot(coor2,U_ex,'k')
xlabel('Radial coordinate','FontSize',12)
ylabel('Displacement u_r','FontSize',12)
legend('FEM','exact') 
set(gca,'Color','w')                         % sets color of current axes

