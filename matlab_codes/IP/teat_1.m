P = [130e3 10e3 0.3 5e3]; % The units of length, force and modulus is mm, N and Mpa
Pscale = [260e3 20e3 0.6 10e3];
P0 = [0.51 0.51 0.51 0.51];
% P0 = [131e3 10.5e3 0.31 5.1e3];

% normalization 
P = P./Pscale;


[U, K, dKdp, dUdp,coor2, element] = Plate_shear(P);

% IS = find(coor2(:,1)==0);
% Pos_ini = IS(randperm(length(IS),1)); % randomly choose a initial point at right boundary  
Pos_ini = 418;

L0 =7*ones(5);  % initial angles
L2 = [3,1,3,2,3,2];
Q_f = -Q_fiber_2m(Pos_ini,L0,dUdp,coor2,0.05,0.00625);