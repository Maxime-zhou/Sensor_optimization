clear 
P = [10 0.3];           % Young and Poisson
% P0 = [10.1 0.31];       % prior guess of parameters
% P = [130e3 10e3 0.3 5e3]; % The units of length, force and modulus is mm, N and Mpa

% normalization 
% Pscale = [260e3 20e3 0.6 10e3];
Pscale = [20 0.6];
P = P./Pscale;
% P0 = [0.51 0.51 0.51 0.51];
P0 = [0.51 0.51];
% P0 = [131e3 10.5e3 0.31 5.1e3];
alpha = 1e-2;            % dimensionless regularization paramete

[U, coor, ndof,dUdp1] = Plate_shear(P); % observed data without noise

dUdp1_x = dUdp1(1:2:end);

mu = 0.1*mean(abs(U)); % mean value of noise
sigma = 0.2*mu; 
R = diag(sigma.^2); % covariance matrix of noise
U_noise = repmat(mu,size(U,1),1) + randn(size(U,1),2)*R; 
U = U + U_noise; % with noise


% calculate the derivates of U respect to P
dUdp = zeros(size(U,1)*ndof,length(P));
for i = 1:length(P)
    dpi = 1e-3;
    dp = zeros(1,length(P));
    dp(i) = dpi;
    % dp = [1 0.1 1e-4 1e-2];
    dUdp(:,i) = reshape((Plate_shear(P+dp) - Plate_shear(P)),size(U,1)*2,1)./dpi;
    
    % dUdp = reshape((Plate_shear(P+dp) - Plate_shear(P)),760*2,1)/dp(2);
end

% L = randperm(size(U,1),100); % choose 100 measurement points randomly
L = [3,4];
% L1 = 1:100;
L1 = [20 31 58 76];
Q0 = dUdp'*dUdp;
Q = dUdp(L,:)'*dUdp(L,:) + dUdp(L+size(U,1),:)'*dUdp(L+size(U,1),:); % Fisher information matrix
Q1 = dUdp(L1,:)'*dUdp(L1,:) + dUdp(L1+size(U,1),:)'*dUdp(L1+size(U,1),:); 
% Q = dUdp(L,2)'*dUdp(L,2) + dUdp(L+size(U,1),2)'*dUdp(L+size(U,1),2); % Fisher information matrix
% Q1 = dUdp(L1,2)'*dUdp(L1,2) + dUdp(L1+size(U,1),2)'*dUdp(L1+size(U,1),2); 
b = trace(Q);   
b1 = trace(Q1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N0 = 10;
% N = 0;
% L = [];
% for i = 1:N0
%     N = N+1;
%     for j = 1:size(U,1)
%         
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adjoint method to compute the gradient of cost function
% L0 = 1:size(U);
% C = nchoosek(L0',4)
% L0 = randperm(size(U,1),10);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Uobs = U(L,:);
dFdE = df_misfit(P,P0,alpha,Uobs,L,dUdp1);
% parameter identification 
% Pini = [100e3 8e3 0.4 3e3];
Pini = [0.4 0.4];
% Pini = [0.4 0.4 0.4 0.4];

[Psol,fval]=fminsearch(@(P) misfit(P,P0,alpha,Uobs,L),Pini,optimset('Display','iter' ...
                ,'TolFun',1e-6,'GradObj','off'))

Psol = Psol.*Pscale