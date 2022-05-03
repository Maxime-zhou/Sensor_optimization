
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Element contributution to Gtheta: T3
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function GBe=T3_2A_solid_gtheta(X,mate,Ue,Ge)

sigma=T3_2A_solid_Sg(X,mate,Ue);             % computes stresses in element
x11=X(1,1); x12=X(1,2);                     % coordinates of first node
x21=X(2,1); x22=X(2,2);                      % coordinates of second node
x31=X(3,1); x32=X(3,2);                      % coordinates of third node
S=((x21-x11)*(x32-x12)-...                   % area of element
   (x31-x11)*(x22-x12))/2;
G=[x22-x32,0,x32-x12,0,x12-x22,0;            % gradient matrix
   x31-x21,0,x11-x31,0,x21-x11,0;
   0,x22-x32,0,x32-x12,0,x12-x22;
   0,x31-x21,0,x11-x31,0,x21-x11]/(2*S);
G_d=(G*Ue);                                  % displacement gradient
G_t=(G*Ge);                                  % theta gradient
prod_G=[G_d(1)*G_t(1)+G_d(2)*G_t(3) ...      % contraction of the two gradients
        G_d(1)*G_t(2)+G_d(2)*G_t(4) ...
        G_d(3)*G_t(1)+G_d(4)*G_t(3) ...
        G_d(3)*G_t(2)+G_d(4)*G_t(4)];
div_t=G_t(1)+G_t(4);
eps=[G_d(1) G_d(4) G_d(2)+G_d(3)];           % strains in engineering notation
symprod_G=[prod_G(1) prod_G(4) ...           % symmetric part of prod_G 
           prod_G(2)+prod_G(3)];
GBe=(symprod_G*sigma'-.5*eps*sigma'*div_t)*S;


