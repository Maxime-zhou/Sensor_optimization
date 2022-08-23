clear 
% isotrpic case
% P = [10 0.3];           % Young and Poisson
% Pscale = [20 0.6];
% P0 = [0.51 0.51];
% P0 = [10.1 0.31];       % prior guess of parameters

% orthotropic case
P = [50e3 50e3 0.3 5e3]; % The units of length, force and modulus is mm, N and Mpa
Pscale = [260e3 20e3 0.6 10e3];
P0 = [0.51 0.51 0.51 0.51];
% P0 = [131e3 10.5e3 0.31 5.1e3];

% normalization 
P = P./Pscale;

alpha = 1e-5;            % dimensionless regularization paramete

% Plate_shear(P);

[U, K, dKdp, dUdp, coor, element] = Plate_shear_c(P); % observed data without noise, Dn is the strain

% dUdp1_x = dUdp1(1:2:end,:);
% dUdp1_y = dUdp1(2:2:end,:);

%-------------------------------------------------------------
% add a noise
%
% mu = 0.1*mean(abs(U)); % mean value of noise
% sigma = 0.2*mu; 
% R = diag(sigma.^2); % covariance matrix of noise
% U_noise = repmat(mu,size(U,1),1) + randn(size(U,1),2)*R; 
% U = U + U_noise; % with noise
mu = 0.01*abs(U);
sigma = 0.2*mu;
U_noise = mu + randn*sigma;
U_n = U + U_noise;


%----------------------------------------------------------------------
% calculate the derivates of U with respect to P, numberical method
%
% dUdp = zeros(size(U,1),length(P));
% for i = 1:length(P)
%     dpi = 1e-3;     % add a small perturbation at one component of P each time
% 
%     dp = zeros(1,length(P));
%     dp(i) = dpi;
%     dUdp(:,i) = (Plate_shear(P+dp) - Plate_shear(P))/dpi;
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot dUdp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% for j = 1:length(P)
%     subplot(2,2,j)
%     for i=1:size(element,1)
%       x = coor(element(i,:),1);
%       y = coor(element(i,:),2);
%       dUdpx = dUdp1_x(element(i,:),j);
%       patch(x,y,abs(dUdpx));
%       hold on
%       subtitle(['Sensitivity contour ' num2str(j)]);
%       xlabel('X'); ylabel('Y'); 
%       colorbar
%       set (gca,'DataAspectRatio', [1 1 1])
%     end   
% end
% 
% 
% 
% figure    
% for j = 1:length(P)
%     subplot(2,2,j)
%     for i=1:size(element,1)
%       x = coor(element(i,:),1);
%       y = coor(element(i,:),2);
%       dUdpy = dUdp1_y(element(i,:),j);
%       patch(x,y,abs(dUdpy));
%       hold on
%       title(['Sensitivity contour ' num2str(j)]);
%       xlabel('X'); ylabel('Y'); 
%       colorbar
%       set (gca,'DataAspectRatio', [1 1 1])
%     end
% end

% a plot of dUdp1_x
% [xq,yq] = meshgrid(0:.01:0.5, 0:.02:1);
% vq = griddata(coor(:,1), coor(:,2),dUdp1_x(:,1), xq,yq);
% contour(xq,yq,vq,Fill="on",ShowText="on")
% colorbar
% hold on
% plot3(coor(:,1),coor(:,2),dUdp1_x(:,1),'o')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% single point sensor optimization

L = [33];
% L1 = 1:100;
L1 = [20 31 58 76]; % L indicate the measured dofs.
Q0 = dUdp'*dUdp;    % Fisher information matrix based on all dofs of the finite element mesh
Q = dUdp(L,:)'*dUdp(L,:); % Fisher information matrix based on measured dofs.
Q1 = dUdp(L1,:)'*dUdp(L1,:);
b = trace(Q);   
b1 = trace(Q1);
b_det = det(Q);
%---------------------------------------------------------------------------
% FSSP algo.
%
N0 = 1;
% L_E = FSSP_FIM(L0,N0,U,dUdp(:,1));
% L_nu = FSSP_FIM(L0,N0,U,dUdp(:,2));
% L_global = FSSP_FIM(L0,N0,U,dUdp);


L_E1_tr = FSSP_FIM_Tr(N0,dUdp(:,1));
L_E2_tr = FSSP_FIM_Tr(N0,dUdp(:,2));
L_nu12_tr = FSSP_FIM_Tr(N0,dUdp(:,3));
L_G12_tr = FSSP_FIM_Tr(N0,dUdp(:,4));

% normalize the dUdp
dUdp_n = normalize(dUdp);
L_global_n = FSSP_FIM_Tr(N0,dUdp_n);


% without normalization
L_global = FSSP_FIM_Tr(N0,dUdp);

%--------------------------------------------------------------------
% using det(Q)
L_E1_det = FSSP_FIM_det(N0,dUdp(:,1));
L_E2_det = FSSP_FIM_det(N0,dUdp(:,2));
L_nu12_det = FSSP_FIM_det(N0,dUdp(:,3));
L_G12_det = FSSP_FIM_det(N0,dUdp(:,4));
L_global_det = FSSP_FIM_det(N0,dUdp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% using numerical method to compute the derivative of strain
%
% [Dn, coor, element2] = Plate_shear_Dn(P);
% dDdp = zeros(size(Dn,1),size(Dn,2),length(P));
%  
% for i = 1:length(P)
% 
%     dpi = 1e-3;     % add a small perturbation at one component of P each time
%     dp = zeros(1,length(P));
%     dp(i) = dpi;
%  
%     dDdp(:,:,i) = (Plate_shear_Dn(P+dp) - Dn)/dpi; 
% 
% end
% 
% F = scatteredInterpolant(coor(:,1),coor(:,2),dDdp(:,1,1));
% F2 = scatteredInterpolant(coor(:,1),coor(:,2),dDdp(:,2,1));
% F3 = scatteredInterpolant(coor(:,1),coor(:,2),dDdp(:,3,1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% continuous fiber sensor optimization
%
% using greedy algorithm to find a optimal fiber placement

L_fiber = 0;
Lmax = 2;
L = [];

IS = find(coor(:,2)==0); % index of start node, at the bottom boundary 
L0 = IS(randperm(length(IS),1));
% L0 = 33;
L = [L,L0];
while (L_fiber<Lmax)
    [row,~]=find(element == L(end));

    nodes_nb = unique(element(row,:)); % neighborhood nodes 
    Q_tr = zeros(length(nodes_nb),1);
    for i = 1:length(nodes_nb)
        if (ismember(nodes_nb(i),L))
            Q_tr(i)=0;
        else
            dx = coor(nodes_nb(i),1)-coor(L(end),1);
            dy = coor(nodes_nb(i),2)-coor(L(end),2);
            dL = sqrt(dx^2+dy^2);
            dDdp = ( (dUdp(2*nodes_nb(i)-1,:) - dUdp(2*L(end)-1,:)) * (dx/dL) + ...
                      (dUdp(2*nodes_nb(i),:) - dUdp(2*L(end),:)) * (dy/dL) ) / dL;
            Q_tr(i) = trace(dDdp'*dDdp);
        end
    end

    [~,ind]=max(Q_tr);
    L = [L,nodes_nb(ind)];
    L_fiber = L_fiber+dL;
end
%--------------------------------------------------------------------




%--------------------------------------------------------------------
% First try of Genetic algorithm
% here, the only design valiables are the angles of fiber

% [U2, K2, dKdp2, dUdp2, coor2, element2] = Plate_shear(P);

% IS = find(coor2(:,2)==0);
% Pos_ini = IS(randperm(length(IS),1)); % randomly choose a initial point at bottom boundary  
Pos_ini = 39;

fun = @(L) -Q_fiber (Pos_ini,L,dUdp,coor);
fun_dev = @(L) -Q_fiber_dev(Pos_ini,L,dUdp,coor);



L0 =1*ones(1,10);  % initial angles
L2 = [6 6 6 7 7 7 7 7 7 7];
% [Q_f,dDdp2] = Q_fiber_2m(Pos_ini,L0,dUdp,coor2,0.05,0.00625);
 
% Q_f1 = Q_fiber_dev(Pos_ini,L0,dUdp,coor);
Q_f1_t = Q_fiber_dev(Pos_ini,L2,dUdp,coor);

lb = ones(1,length(L0));
ub = 8*ones(1,length(L0));
intcon = 1:length(L0);  % constrian design variables are integer
options = optimoptions('ga','Display','iter','ConstraintTolerance',1e-6, 'FunctionTolerance',1e-9,'PlotFcn', @gaplotbestf);


[Lsol,fval] = ga(fun_dev,length(L0),[],[],[],[],lb,ub,[],intcon,options)
% [Lsol,fval] = ga(fun,length(L0),[],[],[],[],lb,ub,[],intcon,optimset('Display','iter'))

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a = 1;

Usim = Plate_shear(P);
% Usim = zeros(size(U,1)*2,1);
% Usim(1:2:end) = Usim_temp(:,1);
% Usim(2:2:end) = Usim_temp(:,2);

% Uobs = zeros(size(U,1)*2,1);
% Uobs(1:2:end) = U(:,1);
% Uobs(2:2:end) = U(:,2);
Uobs = U_n;
dJdU = abs(Usim-Uobs);   % dJdU is a function of P. Here compute it at P true value.
% dJdU = zeros(size(U,1)*2,1);
% dJdU(1:2:end) = dJdU_temp(:,1);
% dJdU(2:2:end) = dJdU_temp(:,2);
Ie=find(Usim==0);
% Uobs(Ie) = [];
Usim(Ie) = [];
dJdU(Ie) = [];



% try adjoint method to compute the gradient of cost function. 
phi = K'\(-dJdU);
dJdp = zeros(1,length(P));
for i=1:length(P)
    dJdp(i) = phi'*dKdp(:,:,i)*Usim + alpha*abs(P(i)-P0(i));
end
L = 1:160;
Uobs = U_n(L,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare with adjoint method 
% C = nchoosek(L0',4)
% L0 = randperm(size(U,1),10);
dFdE = df_misfit(P,P0,alpha,Uobs,L,dUdp1);  % derivative with respect to Young's modulus.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% griddata
% surf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameter identification 
% Pini = [100e3 8e3 0.4 3e3];
% Pini = [0.4 0.4]; % isotropic case
Pini = [0.4 0.4 0.4 0.4]; % orthotropic case

L = L_global_det;
Uobs = U_n(L,:);
[Psol,fval]=fminsearch(@(P) misfit(P,P0,alpha,Uobs,L),Pini,optimset('Display','iter' ...
                ,'TolFun',1e-6,'GradObj','off'))

Psol = Psol.*Pscale
