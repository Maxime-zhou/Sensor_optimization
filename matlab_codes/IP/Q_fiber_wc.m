function Q = Q_fiber_wc(Pos_ini,L,dUdp,coor2,dx,ddx)
% compute the Q(sensor preformance matric) with constrains
if (Check_L(L))
    Q = fiber_path(Pos_ini,L,dUdp,coor2,dx,ddx);
else
    Q = 0;
end
