function hxPerf = exportHxPerf(TCinf, GeoCondition, BDCondition, ...
    h_R_in, h_R_out, p_R_in, p_R_out, ...
    T_MA_in, T_MA_out, mdot_R, mdot_MA, dp_tube, ...
    N, heatPaths, pdropPaths, residual_max, dp_loop_max, time, Prop_handle)
% 导出标准化换热器性能数据结构 hxPerf
% 供 fmincon / 遗传算法等优化器直接调用

hxPerf = struct();

% ---- 几何信息 ----
hxPerf.Geo = struct();
hxPerf.Geo.row        = GeoCondition.row;
hxPerf.Geo.col        = GeoCondition.col;
hxPerf.Geo.Tube_num   = GeoCondition.Tube_num;
hxPerf.Geo.L           = GeoCondition.L;
hxPerf.Geo.D_inner     = GeoCondition.D_inner;
hxPerf.Geo.D_outer     = GeoCondition.D_outer;
hxPerf.Geo.A_R         = GeoCondition.A_R;
hxPerf.Geo.A_MA        = GeoCondition.A_MA;
hxPerf.Geo.Fin_pitch   = GeoCondition.Fin_pitch;

% ---- 拓扑信息 ----
hxPerf.Topo = struct();
hxPerf.Topo.TC_matrix       = TCinf.TC_matrix;
hxPerf.Topo.FlowDirection   = TCinf.FlowDirection;
hxPerf.Topo.inlet_num       = TCinf.inlet_num;
hxPerf.Topo.outlet_num      = TCinf.outlet_num;
hxPerf.Topo.heatPaths       = heatPaths;
hxPerf.Topo.pdropPaths      = pdropPaths;
hxPerf.Topo.n_loops         = size(N, 2);

% ---- 边界条件 ----
hxPerf.Boundary = struct();
hxPerf.Boundary.h_R_inlet       = BDCondition.BD_R.h_R_inlet;
hxPerf.Boundary.p_R_inlet       = BDCondition.BD_R.p_R_inlet;
hxPerf.Boundary.mdot_R_inlet    = BDCondition.BD_R.mdot_R_inlet;
hxPerf.Boundary.T_MA_inlet      = BDCondition.BD_MA.T_MA_inlet;
hxPerf.Boundary.p_MA_inlet      = BDCondition.BD_MA.p_MA_inlet;
hxPerf.Boundary.mdot_MA_inlet   = BDCondition.BD_MA.mdot_MA_inlet;

% ---- 流动结果 ----
hxPerf.Flow = struct();
hxPerf.Flow.mdot_R     = mdot_R;
hxPerf.Flow.mdot_MA    = mdot_MA;
hxPerf.Flow.dp_tube    = dp_tube;

% ---- 热力场 ----
hxPerf.Thermal = struct();
hxPerf.Thermal.h_R_in  = h_R_in;
hxPerf.Thermal.h_R_out = h_R_out;
hxPerf.Thermal.p_R_in  = p_R_in;
hxPerf.Thermal.p_R_out = p_R_out;
hxPerf.Thermal.T_MA_in  = T_MA_in;
hxPerf.Thermal.T_MA_out = T_MA_out;

% ---- 性能汇总 ----
Tube_num = GeoCondition.Tube_num;
heatload_tube = zeros(1, Tube_num);
for t = 1:Tube_num
    heatload_tube(t) = abs((h_R_in(1,t) - h_R_out(end,t)) * mdot_R(t) * 1000);
end

T_R_out = zeros(1, Tube_num);
x_R_out = zeros(1, Tube_num);
for t = 1:Tube_num
    [~, ~, ~, T_R_out(t), x_R_out(t)] = Prop1(p_R_out(end,t)/1e6, h_R_out(end,t), Prop_handle);
end

hxPerf.Performance = struct();
hxPerf.Performance.Q_total       = sum(heatload_tube);
hxPerf.Performance.Q_tube        = heatload_tube;
hxPerf.Performance.dp_total      = (max(p_R_in(1,:)) - min(p_R_out(end,:))) * 1e6;
hxPerf.Performance.T_R_out       = T_R_out;
hxPerf.Performance.x_R_out       = x_R_out;
hxPerf.Performance.residual_max   = residual_max;
hxPerf.Performance.dp_loop_max    = dp_loop_max;
hxPerf.Performance.time           = time;

fprintf('hxPerf 结构体已生成，包含字段: Geo, Topo, Boundary, Flow, Thermal, Performance\n');
end
