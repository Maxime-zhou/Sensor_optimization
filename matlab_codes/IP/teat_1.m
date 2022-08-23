clear
addpath(genpath('../'))   
P = [130e3 10e3 0.3 5e3]; % The units of length, force and modulus is mm, N and Mpa
Pscale = [260e3 20e3 0.6 10e3];
P0 = [0.51 0.51 0.51 0.51];


% P0 = [131e3 10.5e3 0.31 5.1e3];

% normalization 
P = P./Pscale;

MyGeometry_rect_dx01;
coor_c = msh.POS;
element = msh.QUADS(:,1:4);

% [U, K, dKdp, dUdp, coor, element] = Plate_shear_c(P);
[U2, K2, dKdp2, dUdp2, coor2, element2] = Plate_shear(P);

% IS = find(coor2(:,2)==0);
% Pos_ini = IS(randperm(length(IS),1)); % randomly choose a initial point at bottom boundary  
% Pos_ini = 418;
% Pos_ini = 385;

Pos_inic = 29;


e = 1e-6;
P_ind_temp = find(abs(coor_c(Pos_inic,1)-coor2(:,1))<e);
Ind = find(abs(coor_c(Pos_inic,2)-coor2(P_ind_temp,2))<e);
Pos_ini = P_ind_temp(Ind);
%% greedy algorithm
L0 = 29;
Lmax = 1;
L = greedy_tr(Lmax,element,coor_c,coor2,L0,8,dUdp2);
% L = greedy_det(Lmax,element,coor_c,coor2,L0,8,dUdp2);
%% genetic algorithm
L0 =3*ones(1,1);  % initial angles
% L2 = [1 7 4];
L2 = [1     8     2     4     5     5     7     7     7     7];
% best result in coarse mesh
% L2 =  [2 7 4 6 3];
% L2 = [3,6,4,7,2,3,3,3,1,1];
% L2 = [7     7     2     7     2     8     2     1     3     5];

dx = 0.1;
ddx = 0.0125;
% fun = @(L) -Q_fiber_2m(Pos_ini,L,dUdp2,coor2,dx,ddx);
fun = @(L) -Q_fiber_wc(Pos_ini,L,dUdp2,coor2,dx,ddx);
lb = ones(1,length(L0));
ub = 8*ones(1,length(L0));
intcon = 1:length(L0);  % constrian design variables are integer
options = optimoptions('ga','PopulationSize',100,'CrossoverFraction',0.9,'Display','diagnose', ...
      'MaxStallGenerations',1000,'ConstraintTolerance',1e-6, 'FunctionTolerance',1e-9,'PlotFcn', @gaplotbestf);

options = optimoptions(options,'MaxGenerations',5000,'EliteCount',2);
options = optimoptions(options,'InitialPopulationMatrix',[1]);
% options = optimoptions(options,'PlotFcn','gaplotselection');
[Lsol,fval] = ga(fun,length(L0),[],[],[],[],lb,ub,[],intcon,options)
 
% [Q_f1_t,dDdp1] = Q_fiber_dev(Pos_ini,L0,dUdp,coor2);

% Q_f = -Q_fiber_wc(Pos_ini,L2,dUdp2,coor2,dx,ddx);

% [Q_f,dDdp2] = Q_fiber_2m(Pos_ini,L0,dUdp,coor2,0.05,0.00625);
% Q_c = Q_fiber_dev(Pos_inic,L2,dUdp,coor);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% using numerical method to compute the derivative of strain
[Dn, coor, element2] = Plate_shear_Dn(P);
dDdp = zeros(size(Dn,1),size(Dn,2),length(P));
 
for i = 1:length(P)

    dpi = 1e-3;     % add a small perturbation at one component of P each time
    dp = zeros(1,length(P));
    dp(i) = dpi;
    dDdp(:,:,i) = (Plate_shear_Dn(P+dp) - Dn)/dpi; 

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

