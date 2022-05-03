
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% element capacity matrix for linear 2 nodes element and spherical symmetry
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Me=B2_1S_therm_Me(X,c,lumped)

r1=X(1); r2=X(2);                            % radial coordinates of nodes
h=r2-r1;                                     % element length
Me=c*h/60*[12*r1^2+6*r1*r2+2*r2^2, ...       % elementary mass matrix
           3*r1^2+4*r1*r2+3*r2^2;
           3*r1^2+4*r1*r2+3*r2^2, ...
           2*r1^2+6*r1*r2+12*r2^2];
if lumped==1
 Me=diag(sum(Me));                           % diagonal matrix
end

