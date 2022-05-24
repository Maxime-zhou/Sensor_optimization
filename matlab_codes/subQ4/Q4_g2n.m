%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Extrapolation from Gauss points to nodes: Q4 
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function qN=Q4_g2n(qG)

aN=[-1 1 1 -1; ...                 % value of a at nodes
    -1 -1 1 1];

for n=1:4                          % for the 4 nodes
 
 a1=aN(1,n);
 a2=aN(2,n);
 
 f1=1/2*(1-sqrt(3)*a1);            % 1D functions requied to build
 f2=1/2*(1+sqrt(3)*a1);            % 2D functions by cartesian product
 g1=1/2*(1-sqrt(3)*a2);
 g2=1/2*(1+sqrt(3)*a2);
 
 qN(n,1:3)=qG(1,:)*f1*g1+qG(2,:)*f1*g2+ ...
           qG(3,:)*f2*g1+qG(4,:)*f2*g2;
 
end


