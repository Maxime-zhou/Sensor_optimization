function J = misfit(P,P0,alpha,Uobs,L)
U = Plate_shear(P);
U = U(L,:);
J = sum(1/2*(U-Uobs).^2,'all') + sum(1/2*alpha*(P-P0).^2,'all');
