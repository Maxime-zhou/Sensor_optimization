%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% heat-capacity matrix for T6 AXI
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function Me=T6_AX_therm_Me(X,mate)

rhoc=mate(2);

a1=.3333333333333;                           % seven points Gauss rule
a2=.4701420641051;
a3=.1012865073235;
a_gauss=[a1 a1 a1; 
         a2 a2 1-2*a2; a2 1-2*a2 a2; 1-2*a2 a2 a2;
         a3 a3 1-2*a3; a3 1-2*a3 a3; 1-2*a3 a3 a3];       
w1=.1125;                   
w2=0.066197076394253;
w3=0.062969590272414;
w_gauss=[w1 w2 w2 w2 w3 w3 w3];              % Gauss weights 
                  
%a_gauss=1/6*[4 1 1; 1 4 1; 1 1 4];           % three points Gauss abscissae 
%w_gauss=[1/6 1/6 1/6];                       % three points Gauss weights 

Me=zeros(6,6);
for g=1:length(w_gauss),                     % loop over Gauss points
  a=a_gauss(g,:);                            % param. coordinates for gauss point
  D=[4*a(1)-1 0 -4*a(3)+1 4*a(2)...          % derivative of shape functions...
      -4*a(2) 4*(a(3)-a(1));                 % w.r.t. a_1,a_2
      0 4*a(2)-1 -4*a(3)+1 4*a(1) ...
      4*(a(3)-a(2)) -4*a(1)]';
  J=X'*D;                                    % jacobian matrix
  detJ=J(1,1)*J(2,2)-J(1,2)*J(2,1);          % jacobian
  N=[a(1)*(2*a(1)-1) a(2)*(2*a(2)-1) ...     % shape functiona
    a(3)*(2*a(3)-1) 4*a(1)*a(2) ...
    4*a(2)*a(3) 4*a(1)*a(3)];
  r=N*X(:,1);                                % r coord
  Me=Me+rhoc*N'*N*detJ*w_gauss(g)*r;         % mass matrix
end
