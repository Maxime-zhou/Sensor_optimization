function L = greedy_det(Lmax,element,coor,coor2,L0,n,dUdp)
L = [];
L = [L,L0];
L_fiber = 0;
e = 1e-6;
while (L_fiber<Lmax)
    [row,~]=find(element == L(end));   
    nodes_nb = unique(element(row,:)); % neighborhood nodes 
    
    Q_det = zeros(length(nodes_nb),1);
    dDdp = zeros(n,size(dUdp,2));
    P_ind = zeros(n+1,1);
    for i = 1:length(nodes_nb)
        if (ismember(nodes_nb(i),L))
            Q_det(i)=0;
        else
            ddx = ( coor(nodes_nb(i),1)-coor(L(end),1) ) / n;
            ddy = ( coor(nodes_nb(i),2)-coor(L(end),2) ) / n;
            ddL = sqrt(ddx^2+ddy^2);
            
            coor_tem = coor(L(end),:);

            P_ind_temp = find(abs(coor2(:,1)-coor_tem(1))<e);
            Ind =  find(abs(coor2(P_ind_temp,2)-coor_tem(2))<e);
            P_ind(1) = P_ind_temp(Ind);
            for j = 1:n
                coor_tem(1) = coor_tem(1) + ddx; 
                coor_tem(2) = coor_tem(2) + ddy;

                P_ind_temp = find(abs(coor2(:,1)-coor_tem(1))<e);
                Ind =  find(abs(coor2(P_ind_temp,2)-coor_tem(2))<e);
                P_ind(j+1) = P_ind_temp(Ind);
                       
                dDdp(j,:) = ( (dUdp(2*P_ind(j+1)-1,:) - dUdp(2*P_ind(j)-1,:)) * (ddx/ddL) + ...
                      (dUdp(2*P_ind(j+1),:) - dUdp(2*P_ind(j),:)) * (ddy/ddL) ) / ddL; 
            end

            Q_det(i) = det(dDdp'*dDdp);

        end

    end

    [~,ind]=max(Q_det);
    L = [L,nodes_nb(ind)];
    dx = coor(L(end),1) - coor(L(end-1),1);
    dy = coor(L(end),2) - coor(L(end-1),2);
    dL = sqrt(dx^2+dy^2);
    L_fiber = L_fiber+dL;
end