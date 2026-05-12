% Postprocessing

% 每根管的热负荷
heatload_tube = ((h_R_in(1,:) - h_R_out(end,:))').*mdot_R*1000;   % W

% 每根管的压力损失 为dp_tube

% 总热负荷
heatload = sum(heatload_tube)

% 总压降
dp_total = (1 - min(p_R_out(end,:)))*1e6