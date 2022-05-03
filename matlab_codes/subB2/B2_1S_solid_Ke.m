
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% element stiffness matrix for linear 2 nodes element and spherical symmetry
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function Ke=B2_1S_solid_Ke(X,lambda,mu);

r1=X(1); r2=X(2);                            % radial coordinates of nodes
h=r2-r1;                                     % length of element
Ke=(lambda+2*mu)*...                         % element stiffness matrix 
   (r2^2+r1*r2+r1^2)/(3*h)*[1 -1;-1 1]+...   
   2*lambda/6*...
   [-2*(h+3*r1) -h; -h 2*(2*h+3*r1)]+...
   4*(lambda+mu)*h/6*[2 1;1 2];


