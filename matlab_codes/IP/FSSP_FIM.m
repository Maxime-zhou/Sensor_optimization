% Forward sequential sensor placement (FSSP)
function L = FSSP_FIM(L,N0,dUdp)

Q_trace = zeros(size(dUdp,1),1);
for n = 1:N0
%     Q_trace = zeros(size(U,1),1);
    for i = 1:size(dUdp,1)
        if (~ismember(i,L))
            L_temp = [L,i];
%             Q = dUdp(L_temp,:)'*dUdp(L_temp,:)+...
%             dUdp(L_temp+size(U,1),:)'*dUdp(L_temp+size(U,1),:);
            Q = dUdp(L_temp,:)'*dUdp(L_temp,:);
            Q_trace(i) = trace(Q);
%             Q_trace(i) = det(Q);
        end
    end  
    [Q_trace_max,ind] = max(Q_trace);
    L = [L, ind];
end