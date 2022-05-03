
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Extrapolation from Gauss points to nodes: T6 
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function qN=T6_g2n(qG)

%N=1/9*[2 -1 -1 4 1 4;
%       -1 2 -1 4 4 1;
%       -1 -1 2 1 4 4];
%xg=[ones(3,1) N*T];                     % physical coordinates of gauss points
%M=[ones(6,1) T]; 
%qN=M*(xg\qG);                           % linear extrapolation to nodes

aN=[1 0 0 1/2 0 1/2; ...                 % value of a at nodes
    0 1 0 1/2 1/2 0; 0 0 1 0 1/2 1/2];
N=1/3*[5 -1 -1; -1 5 -1; -1 -1 5];       % coefficients of Mg
fN=N*aN;                                 % values of Mg at nodes
qN(1,:)=qG(1,:)*fN(1,1)+qG(2,:)*fN(2,1)+qG(3,:)*fN(3,1);
qN(2,:)=qG(1,:)*fN(1,2)+qG(2,:)*fN(2,2)+qG(3,:)*fN(3,2);
qN(3,:)=qG(1,:)*fN(1,3)+qG(2,:)*fN(2,3)+qG(3,:)*fN(3,3);
qN(4,:)=qG(1,:)*fN(1,4)+qG(2,:)*fN(2,4)+qG(3,:)*fN(3,4);
qN(5,:)=qG(1,:)*fN(1,5)+qG(2,:)*fN(2,5)+qG(3,:)*fN(3,5);
qN(6,:)=qG(1,:)*fN(1,6)+qG(2,:)*fN(2,6)+qG(3,:)*fN(3,6);


