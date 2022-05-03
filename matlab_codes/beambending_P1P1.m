
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Beam bending (displacement approach - exact and reduced integration)
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

addpath(genpath('.'))

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% user input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lb=10;                                    % slenderness ratio

fd=1;                                     % value of load
L =10;                                    % Beam length
NE=10;                                    % number of elements
b=L/Lb;                                   % slenderness ratio
a=1;
E=200000.;                                % Young modulus
nu=.3;                                    % Poisson coefficient
mu=E/(2*(1+nu));                          % definition of Lame constants
I = a*(b^3)/12;                           % area moment of inertia 
S = a*b;                                  % section area
EI = E*I;
kmuS=mu*S;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pre-processor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=L/NE;                                   % size of one element
NN=NE+1;                                  % number of nodes
coor=[0:h:L]';                            % nodal coordinates: uniform discretization 
neq=2*(NN-1);                             % number of unknown nodal values

U=zeros(2*NN,1);                          % initialisation of displacements
U(1:2)=0;                                 % puts imposed displacement on first node

UR=zeros(2*NN,1);                         % initialisation of displacements (reduced version)
UR(1:2)=0;                                % puts imposed displacement on first node (reduced version)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assemblage of stiffness matrix and rhs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K=zeros(neq,neq);                         % allocates stiffness matrix (exact integration)
KR=zeros(neq,neq);                        % allocates stiffness matrix (reduced integration)
F=zeros(neq,1);                           % allocates rhs side
for e=1:NE                                % assemblage of stiffness matrix
 X=[coor(e) coor(e+1)];                   % creates segment
 Ke =B2_P1P1_beam_Ke(X,EI,kmuS);          % element stiffness matrix
 KRe=B2_P1P1R_beam_Ke(X,EI,kmuS);         % element stiffness matrix (reduced integration)
 Fe=B2_P1P1_beam_Fe(X,fd);                % rhs
 Ie=[2*e-3:1:2*e];                        % sets nodal degrees of freedom
 if e>1                                   % if not first element 
  K(Ie,Ie)=K(Ie,Ie)+Ke;                   % the whole Ke  goes into K 
  KR(Ie,Ie)=KR(Ie,Ie)+KRe;                % the whole KRe goes into KR 
  F(Ie)=F(Ie)+Fe;                         % the whole Fe goes into F 
 else
  K(1:2,1:2)=K(1:2,1:2)+Ke(3:4,3:4);      % else only the last coefficients
  KR(1:2,1:2)=KR(1:2,1:2)+KRe(3:4,3:4);
  F(1:2)=F(1:2)+Fe(3:4);                  
 end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solution phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U(3:2*NN)=K\F;                                % solution of linear system
UR(3:2*NN)=KR\F;                              % solution of linear system

clear K KR F                                  % clears memory

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% post phase: evaluates exact solution on a finer grid (Euler Bernoulli beam theory)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coor2=[0:L/100:L];                   % finer grid
fEI=fd/EI;
U_ex=fEI.*((1/24).*coor2.^4 -(L/6).*coor2.^3 + (L^2/4).*coor2.^2);  % exact solution

% Tip deflection
TIPUex=fEI*(1/8)*L^4;
TIPU  =U(2*NN-1);
TIPUR =UR(2*NN-1);
%---------------------------------------------------------------------------
disp(['Slenderness ratio          (L/b):', num2str(Lb)])
disp(['Tip deflection (Euler-Bernoulli):', num2str(TIPUex)])
disp(['Tip deflection      (Ke^{exact}):', num2str(TIPU)])
disp(['Tip deflection        (Ke^{red}):', num2str(TIPUR)])
%---------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% post phase: plots comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

UEF= U(1:2:2*NN);
UREF= UR(1:2:2*NN);
bkgcol=1.*[1 1 1];
figure('Color',bkgcol)                       % plots numerical and exact radial displ.
plot(coor2,U_ex,'k-')
hold on
plot(coor,UEF,'ko','MarkerSize',8)
hold on
plot(coor,UREF,'kx','MarkerSize',8)
xlabel('Beam axis','FontSize',12)
ylabel('Transverse displacement (u)','FontSize',12)
legend('exact','FEM_{exact}','FEM_{Reduced}') 
set(gca,'Color','w')                         % sets color of current axes

