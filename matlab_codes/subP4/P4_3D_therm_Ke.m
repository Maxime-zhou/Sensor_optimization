
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Linear elastic stiffness matrix for scalar problems: P4
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Ke=P4_3D_therm_Ke(X,mate)

k=mate(1);

n1=cross(X(4,:)-X(2,:),X(3,:)-X(2,:));
n2=cross(X(4,:)-X(3,:),X(1,:)-X(3,:));
n3=cross(X(4,:)-X(1,:),X(2,:)-X(1,:));
n4=cross(X(2,:)-X(1,:),X(3,:)-X(1,:));
vol6=(X(4,:)-X(1,:))*n4';
G=[n1; n2; n3; n4]/vol6;
Ke=k*G*G'*vol6/6;