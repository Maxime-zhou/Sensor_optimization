% Forward sequential sensor placement (FSSP) uisng det(FIM)
% as optimization criterion.
function L = FSSP_FIM_det(N0,dUdp)

Q_det = zeros(size(dUdp,1),1);
L = randperm(size(dUdp,1),size(dUdp,2)-1);
for n = 1:N0-size(dUdp,2)+1
%     Q_trace = zeros(size(U,1),1);
    for i = 1:size(dUdp,1)
        if (~ismember(i,L))
            L_temp = [L,i];
%             Q = dUdp(L_temp,:)'*dUdp(L_temp,:)+...
%             dUdp(L_temp+size(U,1),:)'*dUdp(L_temp+size(U,1),:);
            Q = dUdp(L_temp,:)'*dUdp(L_temp,:);
            Q_det(i) = det(Q);
        end
    end  
    [Q_det_max,ind] = max(Q_det);
    L = [L, ind];
end