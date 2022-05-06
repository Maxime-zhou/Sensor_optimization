function dF = df_misfit(P,P0,alpha,Uobs,L,dUdp)
U = Plate_shear(P);
U = U(L,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U1 = zeros(size(L,2)*2,1);
U1(1:2:end) = U(1);
U1(2:2:end) = U(2);
Uobs1 = zeros(size(L,2)*2,1);
Uobs1(1:2:end) = Uobs(1);
Uobs1(2:2:end) = Uobs(2);

dUdp1 = zeros(size(L,2)*2,1);
for i = 1:size(L,2)
    dUdp1(2*i-1:2*i) = dUdp(2*L(i)-1:2*L(i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dF = abs(U1-Uobs1)'*dUdp1 + alpha*abs(P(1)-P0(1));

