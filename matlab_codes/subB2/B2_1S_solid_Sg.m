
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Radial stresses at element center for 2node line element and spherical symmetry
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Se=B2_1S_solid_Sg(X,lambda,mu,Ue);

r1=X(1); r2=X(2);                            % radial coordinates of nodes 
h=r2-r1;                                     % length of element
rm=.5*(r1+r2);                               % average radius
B=1/h*[-1 1];
N=.5*[1 1];
Se=((lambda+2*mu)*B+2*lambda/rm*N)*Ue;

