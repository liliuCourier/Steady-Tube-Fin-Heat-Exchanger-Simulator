function summaryTable(TCinf, GeoCondition, BDCondition, h_R_in, h_R_out, ...
    p_R_in, p_R_out, T_MA_in, T_MA_out, mdot_R, dp_tube, heatPaths, time, Prop_handle)
% 换热器整体性能汇总表：输出关键性能指标到命令行

Tube_num = GeoCondition.Tube_num;
row = GeoCondition.row;
col = GeoCondition.col;

% ---- 基本参数 ----
inlet_num = TCinf.inlet_num;
outlet_num = TCinf.outlet_num;

% ---- 换热量计算 ----
heatload_tube = zeros(1, Tube_num);
for t = 1:Tube_num
    heatload_tube(t) = abs((h_R_in(1,t) - h_R_out(end,t)) * mdot_R(t) * 1000);  % W
end
heatload_total = sum(heatload_tube);

% ---- 压降计算 ----
dp_total = (max(p_R_in(1,:)) - min(p_R_out(end,:))) * 1e6;  % Pa
dp_tube_total = sum(dp_tube) * 1e6;  % Pa (各管压降之和)

% ---- 流量分配 ----
mdot_total = sum(mdot_R);

% ---- 温差 ----
[~, ~, ~, T_R_in] = Prop1(BDCondition.BD_R.p_R_inlet, BDCondition.BD_R.h_R_inlet, Prop_handle);
T_MA_inlet = BDCondition.BD_MA.T_MA_inlet;
T_MA_out_avg = mean(T_MA_out(end, :));

% 各管出口温度
T_R_out = zeros(1, Tube_num);
for t = 1:Tube_num
    [~, ~, ~, T_R_out(t)] = Prop1(p_R_out(end,t)/1e6, h_R_out(end,t), Prop_handle);
end

% 最小传热温差 (pinch point)
dT_tube = zeros(1, Tube_num);
for t = 1:Tube_num
    dT_tube(t) = min(abs(mean(T_MA_in(:,t)) - T_R_out(t)));
end

% ---- 打印汇总表 ----
fprintf('\n');
fprintf('╔══════════════════════════════════════════════╗\n');
fprintf('║     换热器整体性能汇总表 (Demo1.01)          ║\n');
fprintf('╠══════════════════════════════════════════════╣\n');
fprintf('║  基本参数                                    ║\n');
fprintf('║    管排布局: %d排 × %d列 = %d管              ║\n', row, col, Tube_num);
fprintf('║    进口管号: %d   出口管号: %d                ║\n', inlet_num, outlet_num);
fprintf('║    计算耗时: %.2f s                           ║\n', time);
fprintf('╠══════════════════════════════════════════════╣\n');
fprintf('║  热力性能                                    ║\n');
fprintf('║    总换热量: %8.1f W                       ║\n', heatload_total);
fprintf('║    工质进口温度: %6.1f K                     ║\n', T_R_in);
fprintf('║    工质出口温度: %6.1f K (平均)              ║\n', mean(T_R_out));
fprintf('║    空气进口温度: %6.1f K                     ║\n', T_MA_inlet);
fprintf('║    空气出口温度: %6.1f K (平均)              ║\n', T_MA_out_avg);
fprintf('║    最小传热温差: %6.2f K                     ║\n', min(dT_tube));
fprintf('╠══════════════════════════════════════════════╣\n');
fprintf('║  流动性能                                    ║\n');
fprintf('║    工质总流量: %8.4f kg/s                   ║\n', mdot_total);
fprintf('║    工质总压降: %8.1f Pa                     ║\n', dp_total);
fprintf('║    各管压降之和: %6.1f Pa                    ║\n', dp_tube_total);
fprintf('║    各管流量: 均匀分配 (%.4f kg/s/管)         ║\n', mean(mdot_R));
fprintf('╠══════════════════════════════════════════════╣\n');
fprintf('║  各管换热量占比                              ║\n');
for t = 1:Tube_num
    pct = heatload_tube(t) / heatload_total * 100;
    bar_len = round(pct / 2);
    fprintf('║  管%2d: %6.1f W (%5.1f%%) %s║\n', ...
        t, heatload_tube(t), pct, repmat('█', 1, bar_len));
end
fprintf('╚══════════════════════════════════════════════╝\n');
fprintf('\n');

% 管出口状态矩阵
fprintf('\n各管出口状态:\n');
fprintf('  管号  出口焓(kJ/kg)  出口压力(MPa)  出口温度(K)   压降(Pa)   换热量(W)\n');
fprintf('  ----  -------------  --------------  -----------  ---------  ----------\n');
for t = 1:Tube_num
    fprintf('  %4d  %13.2f  %14.6f  %11.2f  %9.2f  %10.1f\n', ...
        t, h_R_out(end,t), p_R_out(end,t)/1e6, T_R_out(t), ...
        dp_tube(t)*1e6, heatload_tube(t));
end
fprintf('\n');
end
