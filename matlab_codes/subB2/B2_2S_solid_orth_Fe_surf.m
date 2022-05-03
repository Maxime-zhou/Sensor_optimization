
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Nodal forces due to surface tractions: B2
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Fe=B2_2S_solid_Fe_surf(X,val,idir)

x11=X(1,1); x12=X(1,2);                      % first node
x21=X(2,1); x22=X(2,2);                      % second node
L=sqrt((x21-x11)^2+(x22-x12)^2);             % element length
NL=L/2*[1 0 1 0;
        0 1 0 1];
TD=zeros(2,1);                               % traction vector
if idir>0,                                   % if cartesian direction
 TD(idir)=val;
else                                         % if along normal
 n=[(x22-x12); -(x21-x11)]/H;                % unit normal vector
 TD=val*n;
end
Fe=NL'*TD;                                   % elementary nodal forces due to tractions

