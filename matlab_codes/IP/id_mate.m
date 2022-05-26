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

alpha = 1e-2;            % dimensionless regularization paramete

[U, K, dKdp, dUdp1] = Plate_shear(P); % observed data without noise
dUdp1_x = dUdp1(1:2:end,:);

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
%  P = [0.8 0.2 0.6 0.6];
dUdp = zeros(size(U,1),length(P));
for i = 1:length(P)
    dpi = 1e-3;     % add a small perturbation at one component of P each time

    dp = zeros(1,length(P));
    dp(i) = dpi;
%     dp = [1 0.1 1e-4 1e-2];
%     dUdp(:,i) = reshape((Plate_shear(P+dp) - Plate_shear(P)),size(U,1)*2,1)/dpi;
    
    dUdp(:,i) = (Plate_shear(P+dp) - Plate_shear(P))/dpi;
end
% dUdp_x = dUdp(1:size(U,1),:);
% L = randperm(size(U,1),100); % choose 100 measurement points randomly
L = [33];
% L1 = 1:100;
L1 = [20 31 58 76];
Q0 = dUdp'*dUdp;
% Q = dUdp(L,:)'*dUdp(L,:) + dUdp(L+size(U,1),:)'*dUdp(L+size(U,1),:); % Fisher information matrix
% Q1 = dUdp(L1,:)'*dUdp(L1,:) + dUdp(L1+size(U,1),:)'*dUdp(L1+size(U,1),:); 
Q = dUdp(L,:)'*dUdp(L,:);
Q1 = dUdp(L1,:)'*dUdp(L1,:);
b = trace(Q);   
b1 = trace(Q1);
b_det = det(Q);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N0 = 10;
% L_E = FSSP_FIM(L0,N0,U,dUdp(:,1));
% L_nu = FSSP_FIM(L0,N0,U,dUdp(:,2));
% L_global = FSSP_FIM(L0,N0,U,dUdp);


% L_E1 = FSSP_FIM_Tr(N0,dUdp(:,1));
% L_E2 = FSSP_FIM_Tr(N0,dUdp(:,2));
% L_nu12 = FSSP_FIM_Tr(N0,dUdp(:,3));
% L_G12 = FSSP_FIM_Tr(N0,dUdp(:,4));
L_global = FSSP_FIM_Tr(N0,dUdp);


%----------------------------------------
% using det(Q)
L_E1_det = FSSP_FIM_det(N0,dUdp(:,1));
L_E2 = FSSP_FIM_det(N0,dUdp(:,2));
L_nu12 = FSSP_FIM_det(N0,dUdp(:,3));
L_G12 = FSSP_FIM_det(N0,dUdp(:,4));
L_global_det = FSSP_FIM_det(N0,dUdp);
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


phi = K'\(-dJdU);
dJdp = zeros(1,length(P));
for i=1:length(P)
    dJdp(i) = phi'*dKdp(:,:,i)*Usim + alpha*abs(P(i)-P0(i));
end
L = 1:160;
Uobs = U_n(L,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adjoint method to compute the gradient of cost function% L0 = 1:size(U);
% C = nchoosek(L0',4)
% L0 = randperm(size(U,1),10);
dFdE = df_misfit(P,P0,alpha,Uobs,L,dUdp1);  % derivative with respect to Young's modulus.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameter identification 
% Pini = [100e3 8e3 0.4 3e3];
Pini = [0.4 0.4];
% Pini = [0.4 0.4 0.4 0.4];

[Psol,fval]=fminsearch(@(P) misfit(P,P0,alpha,Uobs,L),Pini,optimset('Display','iter' ...
                ,'TolFun',1e-6,'GradObj','off'))

Psol = Psol.*Pscale
