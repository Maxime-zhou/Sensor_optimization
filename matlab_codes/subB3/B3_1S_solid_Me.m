
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% element mass matrix for spherical symmetry
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Me=B3_1S_solid_Me(X,rho,lumped);

a1=sqrt((3-2*sqrt(6/5))/7); a2=sqrt((3+2*sqrt(6/5))/7);
a_gauss=[-a2 -a1 a1 a2];
w1=1/36*(18-sqrt(30)); w2=1/36*(18+sqrt(30));
w_gauss=[w1 w2 w2 w1];
Me=zeros(3,3);
for g=1:length(a_gauss),                     % loop over Gauss points
 a=a_gauss(g);                               % param. coordinates for gauss point
 N=[.5*a*(a-1) (1-a^2) .5*a*(1+a)]';
 r=N'*X;
 D=[a-.5 -2*a a+.5]';
 J=X'*D;                                     % jacobian matrix
 Me=Me+rho*r^2*N*N'*J*w_gauss(g);  
end
if lumped==1
 Me=diag(sum(Me));                           % diagonal "vector"
end