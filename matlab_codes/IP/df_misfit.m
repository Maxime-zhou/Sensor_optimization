function dF = df_misfit(P,P0,alpha,Uobs,L,dUdp)
U = Plate_shear(P);
U = U(L,:);
dUdp1 = dUdp(L,1);

dF = abs(U-Uobs)'*dUdp1 + alpha*abs(P(1)-P0(1));
% dF = (U-Uobs)'*dUdp1 + alpha*abs(P(1)-P0(1));
