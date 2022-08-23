function Q = Q_fiber_2m_cc(Pos_ini,L,dUdp,coor2,dx,ddx)

Lb = min(coor2(:,1)); % left boundary coor. of the domain
Rb = max(coor2(:,1)); % right boundary of the domain
Bb = min(coor2(:,2)); % bottom boundary of the domain 
Tb = max(coor2(:,2)); % top boundary of the domain

% dx = coor1(5,1) - coor1(1,1); % mesh size in x direction 
% dy = coor1(44,2) - coor1(1,2); % mesh size in y direction
dy = dx;
ddy = ddx;
e = 1e-6;
n = round(dx/ddx);

P_coor = zeros(length(L)+1,2); % initialize the coordinates of the connect points.
P_ind =  zeros(length(L)+1,1); % the indice of fine mesh nodes of endpoints of segments.
P_inds = zeros(n+1,1);  % the indice of fine mesh nodes in each segments, contain the internal nodes.
% theta = zeros(length(L),1);

P_coor(1,:) = coor2(Pos_ini,:);
P_ind(1) = Pos_ini;
P_inds(1) = Pos_ini;     

dDdp = zeros(length(L)*n,size(dUdp,2));

for i = 1:length(L)
    
    % constrain of curvation 
    if i==1
       
    elseif abs(L(i)-L(i-1))==3 || abs(L(i)-L(i-1))==5
        break
    end

    switch L(i)
        case 1 
            P_coor(i+1,1) = P_coor(i,1)+dx;
            P_coor(i+1,2) = P_coor(i,2);
            if (P_coor(i+1,1)<Lb-e || P_coor(i+1,1)>Rb+e || P_coor(i+1,2)<Bb-e || P_coor(i+1,2)>Tb+e) % outside the domain
                break 
            else
                P_coor2 = zeros(n+1,2); % The coor. of fine mesh nodes of each segments
                P_coor2(1,:) = P_coor(i,:);
                for j = 1:n
                    P_coor2(j+1,1) = P_coor2(j,1)+ddx;
                    P_coor2(j+1,2) = P_coor2(j,2);
                    
                    P_ind_temp =  find(abs(coor2(:,1)-P_coor2(j+1,1))<e);
                    Ind =  find(abs(coor2(P_ind_temp,2)-P_coor2(j+1,2))<e);
                    P_inds(j+1) = P_ind_temp(Ind);
                       
                    dL = ddx;
                    dDdp((i-1)*n+j,:) = ( dUdp(2*P_inds(j+1)-1,:) - dUdp(2*P_inds(j)-1,:) ) / dL;  
                end
                P_ind(i+1) = P_inds(end);
                P_inds(1) = P_inds(end);  %the last point of current segment is the start point of next segment.
            end

        case 2 
            P_coor(i+1,1) = P_coor(i,1)+dx;
            P_coor(i+1,2) = P_coor(i,2)+dy;
            if (P_coor(i+1,1)<Lb-e || P_coor(i+1,1)>Rb+e || P_coor(i+1,2)<Bb-e || P_coor(i+1,2)>Tb+e)
                break 
            else
                P_coor2 = zeros(n+1,2);
                P_coor2(1,:) = P_coor(i,:);
                for j = 1:n
                    P_coor2(j+1,1) = P_coor2(j,1)+ddx;
                    P_coor2(j+1,2) = P_coor2(j,2)+ddy;
                     
                    P_ind_temp =  find(abs(coor2(:,1)-P_coor2(j+1,1))<e);
                    Ind =  find(abs(coor2(P_ind_temp,2)-P_coor2(j+1,2))<e);
                    P_inds(j+1) = P_ind_temp(Ind);

                    dL = sqrt(ddx^2+ddy^2);
                    dDdp((i-1)*n+j,:) =   ( (dUdp(2*P_inds(j+1)-1,:) - dUdp(2*P_inds(j)-1,:)) * ddx/dL ...
                                                   + (dUdp(2*P_inds(j+1),:) - dUdp(2*P_inds(j),:)) * ddy/dL) / dL;    
                end
                P_ind(i+1) = P_inds(end);
                P_inds(1) = P_inds(end);  %the last point of current segment is the start point of next segment.
            end
 
        case 3
            P_coor(i+1,1) = P_coor(i,1);
            P_coor(i+1,2) = P_coor(i,2)+dy;

            if (P_coor(i+1,1)<Lb-e || P_coor(i+1,1)>Rb+e || P_coor(i+1,2)<Bb-e || P_coor(i+1,2)>Tb+e)
                break 
            else 
                P_coor2 = zeros(n+1,2);
                P_coor2(1,:) = P_coor(i,:);
                for j = 1:n
                    P_coor2(j+1,1) = P_coor2(j,1);
                    P_coor2(j+1,2) = P_coor2(j,2)+ddy;
                     
                    P_ind_temp =  find(abs(coor2(:,1)-P_coor2(j+1,1))<e);
                    Ind =  find(abs(coor2(P_ind_temp,2)-P_coor2(j+1,2))<e);
                    P_inds(j+1) = P_ind_temp(Ind);

                    dL = ddy;
                    dDdp((i-1)*n+j,:) = ( dUdp(2*P_inds(j+1),:) - dUdp(2*P_inds(j),:) )  / dL;    
                end
                P_ind(i+1) = P_inds(end);
                P_inds(1) = P_inds(end);  %the last point of current segment is the start point of next segment.
            end
 
        case 4 
            P_coor(i+1,1) = P_coor(i,1)-dx;
            P_coor(i+1,2) = P_coor(i,2)+dy;
            if (P_coor(i+1,1)<Lb-e || P_coor(i+1,1)>Rb+e || P_coor(i+1,2)<Bb-e || P_coor(i+1,2)>Tb+e)
                break
            else
                P_coor2 = zeros(n+1,2);
                P_coor2(1,:) = P_coor(i,:);
                for j = 1:n
                    P_coor2(j+1,1) = P_coor2(j,1)-ddx;
                    P_coor2(j+1,2) = P_coor2(j,2)+ddy;
                   
                    P_ind_temp =  find(abs(coor2(:,1)-P_coor2(j+1,1))<e);
                    Ind =  find(abs(coor2(P_ind_temp,2)-P_coor2(j+1,2))<e);
                    P_inds(j+1) = P_ind_temp(Ind);

                    dL = sqrt(ddx^2+ddy^2);
                    dDdp((i-1)*n+j,:) =   ( (dUdp(2*P_inds(j+1)-1,:) - dUdp(2*P_inds(j)-1,:)) * (-ddx/dL) ...
                                                   + (dUdp(2*P_inds(j+1),:) - dUdp(2*P_inds(j),:)) * ddy/dL) / dL;    
                end
                P_ind(i+1) = P_inds(end);
                P_inds(1) = P_inds(end);  %the last point of current segment is the start point of next segment.
            end

        case 5
            P_coor(i+1,1) = P_coor(i,1)-dx;
            P_coor(i+1,2) = P_coor(i,2);
            
            if (P_coor(i+1,1)<Lb-e || P_coor(i+1,1)>Rb+e || P_coor(i+1,2)<Bb-e || P_coor(i+1,2)>Tb+e)
                break
            else 
                P_coor2 = zeros(n+1,2);
                P_coor2(1,:) = P_coor(i,:);
                for j = 1:n
                    P_coor2(j+1,1) = P_coor2(j,1)-ddx;
                    P_coor2(j+1,2) = P_coor2(j,2);
                    
                    P_ind_temp =  find(abs(coor2(:,1)-P_coor2(j+1,1))<e);
                    Ind =  find(abs(coor2(P_ind_temp,2)-P_coor2(j+1,2))<e);
                    P_inds(j+1) = P_ind_temp(Ind);

                    dL = ddx;
                    dDdp((i-1)*n+j,:) = -( dUdp(2*P_inds(j+1)-1,:) - dUdp(2*P_inds(j)-1,:) )  / dL;    
                end
                P_ind(i+1) = P_inds(end);
                P_inds(1) = P_inds(end);  %the last point of current segment is the start point of next segment.
            end
 
        case 6 
            P_coor(i+1,1) = P_coor(i,1)-dx;
            P_coor(i+1,2) = P_coor(i,2)-dy;
            if (P_coor(i+1,1)<Lb-e || P_coor(i+1,1)>Rb+e || P_coor(i+1,2)<Bb-e || P_coor(i+1,2)>Tb+e)
                break
            else
                P_coor2 = zeros(n+1,2);
                P_coor2(1,:) = P_coor(i,:);
                for j = 1:n
                    P_coor2(j+1,1) = P_coor2(j,1)-ddx;
                    P_coor2(j+1,2) = P_coor2(j,2)-ddy;
                    
                    P_ind_temp =  find(abs(coor2(:,1)-P_coor2(j+1,1))<e);
                    Ind =  find(abs(coor2(P_ind_temp,2)-P_coor2(j+1,2))<e);
                    P_inds(j+1) = P_ind_temp(Ind);

                    dL = sqrt(ddx^2+ddy^2);
                    dDdp((i-1)*n+j,:) =   ( (dUdp(2*P_inds(j+1)-1,:) - dUdp(2*P_inds(j)-1,:)) * (-ddx/dL) ...
                                                   + (dUdp(2*P_inds(j+1),:) - dUdp(2*P_inds(j),:)) * (-ddy/dL) ) / dL;    
                end
                P_ind(i+1) = P_inds(end);
                P_inds(1) = P_inds(end);  %the last point of current segment is the start point of next segment.
            end

        case 7
            P_coor(i+1,1) = P_coor(i,1);
            P_coor(i+1,2) = P_coor(i,2)-dy;
            if (P_coor(i+1,1)<Lb-e || P_coor(i+1,1)>Rb+e || P_coor(i+1,2)<Bb-e || P_coor(i+1,2)>Tb+e) % outside the domain
                break
            else 
                P_coor2 = zeros(n+1,2);
                P_coor2(1,:) = P_coor(i,:);
                for j = 1:n
                    P_coor2(j+1,1) = P_coor2(j,1);
                    P_coor2(j+1,2) = P_coor2(j,2)-ddy;
                   
                    P_ind_temp =  find(abs(coor2(:,1)-P_coor2(j+1,1))<e);
                    Ind =  find(abs(coor2(P_ind_temp,2)-P_coor2(j+1,2))<e);
                    P_inds(j+1) = P_ind_temp(Ind);
                       
                    dL = ddy;
                    dDdp((i-1)*round(n)+j,:) = -( dUdp(2*P_inds(j+1),:) - dUdp(2*P_inds(j),:) ) / dL;  
                end
                P_ind(i+1) = P_inds(end);
                P_inds(1) = P_inds(end);  %the last point of current segment is the start point of next segment.
            end

        case 8
            P_coor(i+1,1) = P_coor(i,1)+dx;
            P_coor(i+1,2) = P_coor(i,2)-dy;
            if (P_coor(i+1,1)<Lb-e || P_coor(i+1,1)>Rb+e || P_coor(i+1,2)<Bb-e || P_coor(i+1,2)>Tb+e)
                break
            else
                P_coor2 = zeros(n+1,2);
                P_coor2(1,:) = P_coor(i,:);

                for j = 1:n
                    P_coor2(j+1,1) = P_coor2(j,1)+ddx;
                    P_coor2(j+1,2) = P_coor2(j,2)-ddy;
                  
                    P_ind_temp =  find(abs(coor2(:,1)-P_coor2(j+1,1))<e);
                    Ind =  find(abs(coor2(P_ind_temp,2)-P_coor2(j+1,2))<e);
                    P_inds(j+1) = P_ind_temp(Ind);

                    dL = sqrt(ddx^2+ddy^2);
                    dDdp((i-1)*n+j,:) =   ( (dUdp(2*P_inds(j+1)-1,:) - dUdp(2*P_inds(j)-1,:)) * ddx/dL ...
                                                   + (dUdp(2*P_inds(j+1),:) - dUdp(2*P_inds(j),:)) * (-ddy/dL) ) / dL;    
                end
                P_ind(i+1) = P_inds(end);
                P_inds(1) = P_inds(end);  %the last point of current segment is the start point of next segment.
            end
 
    end

 
end


% check the existence of repeat points, Q=0 if in that case. 
P_ind = unique(P_ind,'stable'); 

if (length(P_ind)<length(L)+1)
    Q = 0;
elseif (abs(P_coor(end,1)-Lb)<e || abs(P_coor(end,1)-Rb)<e || abs(P_coor(end,2)-Bb)<e || abs(P_coor(end,2)-Tb)<e)
    Q = det(dDdp'*dDdp);
else
    Q = 0;  % if the end point of fiber not at  boundary, Q = 0
%     Q = det(dDdp'*dDdp);
%     Q = trace(dDdp'*dDdp);
end

