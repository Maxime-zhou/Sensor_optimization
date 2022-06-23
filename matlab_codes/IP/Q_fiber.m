function Q = Q_fiber(Pos_ini,L,dUdp,coor)

% theta = 0:pi/4:2*pi; % theta contain 8 angles corresponding to 8 orientation.
Lb = min(coor(:,1)); % left boundary coor. of the domain
Rb = max(coor(:,1)); % right boundary of the domain
Bb = min(coor(:,2)); % bottom boundary of the domain 
Tb = max(coor(:,2)); % top boundary of the domain

dx = coor(5,1) - coor(1,1); % mesh size in x direction 
dy = coor(44,2) - coor(1,2); % mesh size in y direction
a = atan(dx/dy);

theta = [0, a, pi/2, pi-a, pi, pi+a, 3*pi/2, 2*pi-a];
% theta contain 8 angles corresponding to 8 orientation.

P_coor = zeros(length(L)+1,2); % initialize the coordinates of the connect points
P_ind = zeros(length(L)+1,1);  % initialize the indices of the coor.
P_coor(1,:) = coor(Pos_ini,:);
P_ind(1) = Pos_ini;

for i = 1:length(L)
    if (cos(theta(L(i)))==1||cos(theta(L(i)))==-1)
       P_coor(i+1,1) = P_coor(i,1)+dx*cos(theta(L(i)));
       P_coor(i+1,2) = P_coor(i,2);
    elseif(sin(theta(L(i)))==1||sin(theta(L(i)))==-1)
       P_coor(i+1,1) = P_coor(i,1);
       P_coor(i+1,2) = P_coor(i,2)+dy*sin(theta(L(i)));
    else
       P_coor(i+1,1) = P_coor(i,1)+dx*sign(cos(theta(L(i))));
       P_coor(i+1,2) = P_coor(i,2)+dy*sign(sin(theta(L(i))));
    end

    if (P_coor(i+1,1)<Lb || P_coor(i+1,1)>Rb || P_coor(i+1,2)<Bb || P_coor(i+1,2)>Tb)
%         Q = 0;
        break
    end
%     P_ind(i+1) = find(coor(:,1)==P_coor(i+1,1) && coor(:,2)==P_coor(i+1,2));
%     P_ind = coor(:,1)==P_coor(i+1,1) && coor(:,2)==P_coor(i+1,2);
    e = 1e-6;
    P_ind_temp =  find(abs(coor(:,1)-P_coor(i+1,1))<e);
    Ind =  find(abs(coor(P_ind_temp,2)-P_coor(i+1,2))<e);
    P_ind(i+1) = P_ind_temp(Ind);
end


dDdp = zeros(length(L),size(dUdp,2));

% IS = find(~P_ind);
% if (~IS==[])
P_ind = unique(P_ind);
if (length(P_ind)<length(L)+1)
    Q = 0;
elseif(~all(P_ind))  % if any element of P_ind equal to zero, Q is zero
    Q = 0;
else
    for i=1:length(L)
    dL = sqrt((P_coor(i+1,1)-P_coor(i,1))^2 + (P_coor(i+1,2)-P_coor(i,2))^2);
    
    dDdp(i,:) = ( (dUdp(2*P_ind(i+1)-1,:) - dUdp(2*P_ind(i)-1,:)) * cos(theta(L(i))) + ...
                      (dUdp(2*P_ind(i+1),:) - dUdp(2*P_ind(i),:)) * sin(theta(L(i))) ) / dL;
    end
    Q = det(dDdp'*dDdp);
%     Q = trace(dDdp'*dDdp);
end

