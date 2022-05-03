
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Extrapolation from Gauss points to nodes: T4 
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function qN=Q8_g2n(qG)

b=sqrt(3/5); 
aN=[-1 1 1 -1 0 1 0 -1; ...                 % value of a at nodes
     -1 -1 1 1 -1 0 1 0];
for n=1:8
 a1=aN(1,n);
 a2=aN(2,n);
 
 f1=1/(2*b^2)*a1*(a1-b);
 f2=1/b^2*(b-a1)*(b+a1);
 f3=1/(2*b^2)*a1*(a1+b);
 
 g1=1/(2*b^2)*a2*(a2-b);
 g2=1/b^2*(b-a2)*(b+a2);
 g3=1/(2*b^2)*a2*(a2+b);
 
 qN(n,1:3)=qG(1,:)*f1*g1+qG(2,:)*f1*g2+qG(3,:)*f1*g3+ ...
           qG(4,:)*f2*g1+qG(5,:)*f2*g2+qG(6,:)*f2*g3+ ...
           qG(7,:)*f3*g1+qG(8,:)*f3*g2+qG(9,:)*f3*g3;
 
end
