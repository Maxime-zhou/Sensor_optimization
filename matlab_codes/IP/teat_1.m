clear
P = [130e3 10e3 0.3 5e3]; % The units of length, force and modulus is mm, N and Mpa
Pscale = [260e3 20e3 0.6 10e3];
P0 = [0.51 0.51 0.51 0.51];


% P0 = [131e3 10.5e3 0.31 5.1e3];

% normalization 
P = P./Pscale;


coarse_mesh;
coor_c = msh.POS;

[U, K, dKdp, dUdp, coor, element] = Plate_shear_c(P);
[U2, K2, dKdp2, dUdp2, coor2, element2] = Plate_shear(P);

IS = find(coor2(:,1)==0);
% Pos_ini = IS(randperm(length(IS),1)); % randomly choose a initial point at right boundary  
% Pos_ini = 418;
Pos_ini = 385;

e = 1e-6;
P_ind_temp =  find(abs(coor_c(:,1)-coor2(Pos_ini,1))<e);
Ind =  find(abs(coor_c(P_ind_temp,2)-coor2(Pos_ini,2))<e);
P_indc = P_ind_temp(Ind);

Pos_inic = 33;

L0 =3*ones(1,5);  % initial angles
L2 = [7,7,7,7,7];
% L2 =  [3     7     3     3     3];

% [Q_f1_t,dDdp1] = Q_fiber_dev(Pos_ini,L0,dUdp,coor2);
Q_f = -Q_fiber_2m(Pos_ini,L2,dUdp2,coor2,0.05,0.00625);
% [Q_f,dDdp2] = Q_fiber_2m(Pos_ini,L0,dUdp,coor2,0.05,0.00625);
Q_c = Q_fiber_dev(Pos_inic,L2,dUdp,coor);

dx = 0.05;
ddx = 0.00625;
fun = @(L) -Q_fiber_2m(Pos_ini,L,dUdp2,coor2,dx,ddx);
lb = ones(1,length(L0));
ub = 8*ones(1,length(L0));
intcon = 1:length(L0);  % constrian design variables are integer
options = optimoptions('ga','Display','iter','ConstraintTolerance',1e-9, 'FunctionTolerance',1e-9,'PlotFcn', @gaplotbestf);


[Lsol,fval] = ga(fun,length(L0),[],[],[],[],lb,ub,[],intcon,options)
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