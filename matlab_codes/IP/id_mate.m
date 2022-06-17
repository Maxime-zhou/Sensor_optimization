clear 
% isotrpic case
% P = [10 0.3];           % Young and Poisson
% Pscale = [20 0.6];
% P0 = [0.51 0.51];
% P0 = [10.1 0.31];       % prior guess of parameters

% orthotropic case
P = [130e3 10e3 0.3 5e3]; % The units of length, force and modulus is mm, N and Mpa
Pscale = [260e3 20e3 0.6 10e3];
P0 = [0.51 0.51 0.51 0.51];
% P0 = [131e3 10.5e3 0.31 5.1e3];

% normalization 
P = P./Pscale;

alpha = 1e-5;            % dimensionless regularization paramete

% Plate_shear(P);

[U, K, dKdp, dUdp1, coor, element] = Plate_shear(P); % observed data without noise, Dn is the strain
dUdp1_x = dUdp1(1:2:end,:);
dUdp1_y = dUdp1(2:2:end,:);

% mu = 0.1*mean(abs(U)); % mean value of noise
% sigma = 0.2*mu; 
% R = diag(sigma.^2); % covariance matrix of noise
% U_noise = repmat(mu,size(U,1),1) + randn(size(U,1),2)*R; 
% U = U + U_noise; % with noise
mu = 0.01*abs(U);
sigma = 0.2*mu;
U_noise = mu + randn*sigma;
U_n = U + U_noise;

% calculate the derivates of U with respect to P
% P = [0.8 0.2 0.6 0.6];
dUdp = zeros(size(U,1),length(P));
for i = 1:length(P)
    dpi = 1e-3;     % add a small perturbation at one component of P each time

    dp = zeros(1,length(P));
    dp(i) = dpi;
    dUdp(:,i) = (Plate_shear(P+dp) - Plate_shear(P))/dpi;
end

% dUdp_x = dUdp(1:size(U,1),:);
% L = randperm(size(U,1),100); % choose 100 measurement points randomly


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot dUdp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for j = 1:length(P)
%     figure
%     for i=1:size(element,1)
%     x = coor(element(i,:),1);
%     y = coor(element(i,:),2);
%     dUdpx = dUdp1_x(element(i,:),j);
%     patch(x,y,abs(dUdpx));
%     end
% end
% 
% for j = 1:length(P)
%     figure
%     for i=1:size(element,1)
%     x = coor(element(i,:),1);
%     y = coor(element(i,:),2);
%     dUdpy = dUdp1_y(element(i,:),j);
%     patch(x,y,abs(dUdpx));
%     end
% end


%%%%%%%%%%%%%%%%%%%%%%%%%
L = [33];
% L1 = 1:100;
L1 = [20 31 58 76]; % L indicate the measured dofs.
Q0 = dUdp'*dUdp;    % Fisher information matrix based on all dofs of the finite element mesh
Q = dUdp(L,:)'*dUdp(L,:); % Fisher information matrix based on measured dofs.
Q1 = dUdp(L1,:)'*dUdp(L1,:);
b = trace(Q);   
b1 = trace(Q1);
b_det = det(Q);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N0 = 10;
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

%----------------------------------------
% using det(Q)
L_E1_det = FSSP_FIM_det(N0,dUdp(:,1));
L_E2_det = FSSP_FIM_det(N0,dUdp(:,2));
L_nu12_det = FSSP_FIM_det(N0,dUdp(:,3));
L_G12_det = FSSP_FIM_det(N0,dUdp(:,4));
L_global_det = FSSP_FIM_det(N0,dUdp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% using greedy algorithm to find a optimal fiber placement
IS = find(coor(:,1)==0); % index of start node, at the right boundary 
L_fiber = 0;
Lmax = 2;
L = [];
% L0 = IS(randperm(length(IS),1));
L0 = 4;
L = [L,L0];
while (L_fiber<Lmax)
    [row,~]=find(element == L(end));

    nodes_nb = unique(element(row,:)); % neighborhood nodes
%     dDdp = zeros(length(nodes_nb),1);
    Q_tr = zeros(length(nodes_nb),1);
    for i = 1:length(nodes_nb)
        if (ismember(nodes_nb(i),L))
%             dDdp(i)=0;
            Q_tr(i)=0;
        else
            dx = coor(nodes_nb(i),1)-coor(L(end),1);
            dy = coor(nodes_nb(i),2)-coor(L(end),2);
            dL = sqrt(dx^2+dy^2);
            dDdp = ( (dUdp(2*nodes_nb(i)-1,:) - dUdp(2*L(end)-1,:)) * (dx/dL) + ...
                      (dUdp(2*nodes_nb(i),:) - dUdp(2*L(end),:)) * (dy/dL) ) / dL;
            Q_tr(i) =trace(dDdp'*dDdp);
        end
    end
%     [~,ind] = max(dDdp);
    [~,ind]=max(Q_tr);
    L = [L,nodes_nb(ind)];
    L_fiber = L_fiber+dL;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First try of Genetic algorithm
% here, the only design valiables are the angles of fiber
intcon = [1,2];  % constrian design variables are integer
IS = find(coor(:,1)==0);
Pos_ini = IS(randperm(length(IS),1)); % randomly choose a initial point at right boundary  
fun = @(L) Q_fiber(Pos_ini,L,dUdp,coor);

L0 = 2*ones(1,10);  % initial angles
Q_f = Q_fiber(Pos_ini,L0,dUdp,coor);
ga(fun,10,[],[],[],[],1,8,[],intcon)


 
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
% conpare with adjoint method 
% C = nchoosek(L0',4)
% L0 = randperm(size(U,1),10);
dFdE = df_misfit(P,P0,alpha,Uobs,L,dUdp1);  % derivative with respect to Young's modulus.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
