% replace if else to swith case.
function Q = Q_fiber_2m(Pos_ini,L,dUdp,coor2,dx,ddx)

Lb = min(coor2(:,1)); % left boundary coor. of the domain
Rb = max(coor2(:,1)); % right boundary of the domain
Bb = min(coor2(:,2)); % bottom boundary of the domain 
Tb = max(coor2(:,2)); % top boundary of the domain

% dx = coor1(5,1) - coor1(1,1); % mesh size in x direction 
% dy = coor1(44,2) - coor1(1,2); % mesh size in y direction
dy = dx;
ddy = ddx;
e = 1e-6;

P_coor = zeros(length(L)+1,2); % initialize the coordinates of the connect points
P_ind = zeros(dx/ddx+1,1);  % initialize the indices of the coor.
% theta = zeros(length(L),1);

P_coor(1,:) = coor2(Pos_ini,:);
P_ind(1) = Pos_ini;      % need to use athe index of fine mesh.

% if L indicate a overlapping path, set Q = 0
% for i=1:length(L)-1
%     if L(i+1) - L(i) == 4
%         Q = 0;
%     end
% end
dDdp = zeros(length(L)*dx/ddx,size(dUdp,2));

for i = 1:length(L)
    switch L(i)
        case 1 
            P_coor(i+1,1) = P_coor(i,1)+dx;
            P_coor(i+1,2) = P_coor(i,2);
            if (P_coor(i+1,1)<Lb || P_coor(i+1,1)>Rb || P_coor(i+1,2)<Bb || P_coor(i+1,2)>Tb) % outside the domain
                break 
            else 
                P_coor2 = zeros(dx/ddx,2);
                P_coor2(1,:) = P_coor(i,:);
                for j = 1:dx/ddx-1
                    P_coor2(j+1,1) = P_coor2(j,1)+ddx;
                    P_coor2(j+1,2) = P_coor2(j,2);
                    
                    P_ind_temp =  find(abs(coor2(:,1)-P_coor2(j+1,1))<e);
                    Ind =  find(abs(coor2(P_ind_temp,2)-P_coor2(j+1,2))<e);
                    P_ind(j+1) = P_ind_temp(Ind);
                       
                    dL = ddx;
                    dDdp((i-1)*dx/ddx+j,:) =   ( dUdp(2*P_ind(j+1)-1,:) - dUdp(2*P_ind(j)-1,:) ) / dL;  
                end
            end
        case 2 
            P_coor(i+1,1) = P_coor(i,1)+dx;
            P_coor(i+1,2) = P_coor(i,2)+dy;
            if (P_coor(i+1,1)<Lb || P_coor(i+1,1)>Rb || P_coor(i+1,2)<Bb || P_coor(i+1,2)>Tb)
                break 
            else
                P_coor2 = zeros(dx/ddx,2);
                P_coor2(1,:) = P_coor(i,:);
                for j = 1:dx/ddx-1
                    P_coor2(j+1,1) = P_coor2(j,1)+ddx;
                    P_coor2(j+1,2) = P_coor2(j,2)+ddy;
                     
                    P_ind_temp =  find(abs(coor2(:,1)-P_coor2(j+1,1))<e);
                    Ind =  find(abs(coor2(P_ind_temp,2)-P_coor2(j+1,2))<e);
                    P_ind(j+1) = P_ind_temp(Ind);

                    dL = sqrt(ddx^2+ddy^2);
                    dDdp((i-1)*dx/ddx+j,:) =   ( (dUdp(2*P_ind(j+1)-1,:) - dUdp(2*P_ind(j)-1,:)) * ddx/dL ...
                                                   + (dUdp(2*P_ind(j+1),:) - dUdp(2*P_ind(j),:)) * ddy/dL) / dL;    
                end
            end
 
        case 3
            P_coor(i+1,1) = P_coor(i,1);
            P_coor(i+1,2) = P_coor(i,2)+dy;

            if (P_coor(i+1,1)<Lb || P_coor(i+1,1)>Rb || P_coor(i+1,2)<Bb || P_coor(i+1,2)>Tb)
                break 
            else 
                P_coor2 = zeros(dx/ddx,2);
                P_coor2(1,:) = P_coor(i,:);
                for j = 1:dx/ddx-1
                    P_coor2(j+1,1) = P_coor2(j,1);
                    P_coor2(j+1,2) = P_coor2(j,2)+ddy;
                     
                    P_ind_temp =  find(abs(coor2(:,1)-P_coor2(j+1,1))<e);
                    Ind =  find(abs(coor2(P_ind_temp,2)-P_coor2(j+1,2))<e);
                    P_ind(j+1) = P_ind_temp(Ind);

                    dL = ddy;
                    dDdp((i-1)*dx/ddx+j,:) = ( dUdp(2*P_ind(j+1),:) - dUdp(2*P_ind(j),:) )  / dL;    
                end
            end
 
        case 4 
            P_coor(i+1,1) = P_coor(i,1)-dx;
            P_coor(i+1,2) = P_coor(i,2)+dy;
            if (P_coor(i+1,1)<Lb || P_coor(i+1,1)>Rb || P_coor(i+1,2)<Bb || P_coor(i+1,2)>Tb)
                break 
            else
                P_coor2 = zeros(dx/ddx,2);
                P_coor2(1,:) = P_coor(i,:);
                for j = 1:dx/ddx-1
                    P_coor2(j+1,1) = P_coor2(j,1)-ddx;
                    P_coor2(j+1,2) = P_coor2(j,2)+ddy;
                   
                    P_ind_temp =  find(abs(coor2(:,1)-P_coor2(j+1,1))<e);
                    Ind =  find(abs(coor2(P_ind_temp,2)-P_coor2(j+1,2))<e);
                    P_ind(j+1) = P_ind_temp(Ind);

                    dL = sqrt(ddx^2+ddy^2);
                    dDdp((i-1)*dx/ddx+j,:) =   -( (dUdp(2*P_ind(j+1)-1,:) - dUdp(2*P_ind(j)-1,:)) * ddx/dL ...
                                                   + (dUdp(2*P_ind(j+1),:) - dUdp(2*P_ind(j),:)) * ddy/dL) / dL;    
                end
            end
        case 5
            P_coor(i+1,1) = P_coor(i,1)-dx;
            P_coor(i+1,2) = P_coor(i,2);
            
            if (P_coor(i+1,1)<Lb || P_coor(i+1,1)>Rb || P_coor(i+1,2)<Bb || P_coor(i+1,2)>Tb)
                break 
            else 
                P_coor2 = zeros(dx/ddx,2);
                P_coor2(1,:) = P_coor(i,:);
                for j = 1:dx/ddx-1
                    P_coor2(j+1,1) = P_coor2(j,1)-ddx;
                    P_coor2(j+1,2) = P_coor2(j,2);
                    
                    P_ind_temp =  find(abs(coor2(:,1)-P_coor2(j+1,1))<e);
                    Ind =  find(abs(coor2(P_ind_temp,2)-P_coor2(j+1,2))<e);
                    P_ind(j+1) = P_ind_temp(Ind);

                    dL = ddx;
                    dDdp((i-1)*dx/ddx+j,:) = -( dUdp(2*P_ind(j+1)-1,:) - dUdp(2*P_ind(j)-1,:) )  / dL;    
                end
            end
 
        case 6 
            P_coor(i+1,1) = P_coor(i,1)-dx;
            P_coor(i+1,2) = P_coor(i,2)-dy;
            if (P_coor(i+1,1)<Lb || P_coor(i+1,1)>Rb || P_coor(i+1,2)<Bb || P_coor(i+1,2)>Tb)
                break 
            else
                P_coor2 = zeros(dx/ddx,2);
                P_coor2(1,:) = P_coor(i,:);
                for j = 1:dx/ddx-1
                    P_coor2(j+1,1) = P_coor2(j,1)-ddx;
                    P_coor2(j+1,2) = P_coor2(j,2)-ddy;
                    
                    P_ind_temp =  find(abs(coor2(:,1)-P_coor2(j+1,1))<e);
                    Ind =  find(abs(coor2(P_ind_temp,2)-P_coor2(j+1,2))<e);
                    P_ind(j+1) = P_ind_temp(Ind);

                    dL = sqrt(ddx^2+ddy^2);
                    dDdp((i-1)*dx/ddx+j,:) =   -( (dUdp(2*P_ind(j+1)-1,:) - dUdp(2*P_ind(j)-1,:)) * ddx/dL ...
                                                   - (dUdp(2*P_ind(j+1),:) - dUdp(2*P_ind(j),:)) * ddy/dL) / dL;    
                end
            end

        case 7
            P_coor(i+1,1) = P_coor(i,1);
            P_coor(i+1,2) = P_coor(i,2)-dy;
            if (P_coor(i+1,1)<Lb || P_coor(i+1,1)>Rb || P_coor(i+1,2)<Bb || P_coor(i+1,2)>Tb) % outside the domain
                break 
            else 
                P_coor2 = zeros(dx/ddx,2);
                P_coor2(1,:) = P_coor(i,:);
                for j = 1:dx/ddx-1
                    P_coor2(j+1,1) = P_coor2(j,1);
                    P_coor2(j+1,2) = P_coor2(j,2)-ddy;
                   
                    P_ind_temp =  find(abs(coor2(:,1)-P_coor2(j+1,1))<e);
                    Ind =  find(abs(coor2(P_ind_temp,2)-P_coor2(j+1,2))<e);
                    P_ind(j+1) = P_ind_temp(Ind);
                       
                    dL = ddy;
                    dDdp((i-1)*dx/ddx+j,:) =   -( dUdp(2*P_ind(j+1)-1,:) - dUdp(2*P_ind(j)-1,:) ) / dL;  
                end
            end

%             theta(i) = 3*pi/2;
        case 8
            P_coor(i+1,1) = P_coor(i,1)+dx;
            P_coor(i+1,2) = P_coor(i,2)-dy;
            if (P_coor(i+1,1)<Lb || P_coor(i+1,1)>Rb || P_coor(i+1,2)<Bb || P_coor(i+1,2)>Tb)
                break 
            else
                P_coor2 = zeros(dx/ddx,2);
                P_coor2(1,:) = P_coor(i,:);

                for j = 1:dx/ddx-1
                    P_coor2(j+1,1) = P_coor2(j,1)+ddx;
                    P_coor2(j+1,2) = P_coor2(j,2)-ddy;
                  
                    P_ind_temp =  find(abs(coor2(:,1)-P_coor2(j+1,1))<e);
                    Ind =  find(abs(coor2(P_ind_temp,2)-P_coor2(j+1,2))<e);
                    P_ind(j+1) = P_ind_temp(Ind);

                    dL = sqrt(ddx^2+ddy^2);
                    dDdp((i-1)*dx/ddx+j,:) =   ( (dUdp(2*P_ind(j+1)-1,:) - dUdp(2*P_ind(j)-1,:)) * ddx/dL ...
                                                   - (dUdp(2*P_ind(j+1),:) - dUdp(2*P_ind(j),:)) * ddy/dL) / dL;    
                end
            end
 
    end

 
end
Q = det(dDdp'*dDdp);


% dDdp = zeros(length(L),size(dUdp,2));

% IS = find(~P_ind);
% if (~IS==[])
% P_ind = unique(P_ind);
% if (length(P_ind)<length(L)+1)
%     Q = 0;
% elseif(~all(P_ind))  % if any element of P_ind equal to zero, Q is zero
%     Q = 0;
% else
%     for i=1:length(L)
%         dL = sqrt((P_coor(i+1,1)-P_coor(i,1))^2 + (P_coor(i+1,2)-P_coor(i,2))^2);
%     
%         dDdp(i,:) = ( ( dUdp(2*P_ind(i+1)-1,:) - dUdp(2*P_ind(i)-1,:) ) * cos(theta(i) ) + ...
%                       ( dUdp(2*P_ind(i+1),:) - dUdp(2*P_ind(i),:) ) * sin(theta(i) ) ) / dL;
%     end
% %     Q = det(dDdp'*dDdp);
%     Q = trace(dDdp'*dDdp);
% end

