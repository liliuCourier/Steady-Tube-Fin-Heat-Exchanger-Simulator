% PostProcessing — 后处理集成入口 (Demo1.02)
% 仿真完成后运行此脚本，一键生成全部可视化与性能报告

% 添加子模块路径
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'PostProcessing'));
addpath(fullfile(root, 'Lib'));
addpath(fullfile(root, 'Solver'));

%% 1. 沿线物性分布曲线
plotAlongPath(heatPaths, h_R_in, h_R_out, p_R_in, p_R_out, ...
    T_MA_in, T_MA_out, dp_tube, mdot_R, Prop_handle, row, TCinf, predecessors_in);

%% 2. 环路性能对比（有环路时生效）
plotLoopBalance(N, R_flow, mdot0, mdot_R, dp_tube, u0, options);

%% 3. 整体性能汇总表
summaryTable(TCinf, GeoCondition, BDCondition, h_R_in, h_R_out, ...
    p_R_in, p_R_out, T_MA_in, T_MA_out, mdot_R, dp_tube, heatPaths, time, Prop_handle, N);

%% 3.5 各管流量分布
Tube_num = GeoCondition.Tube_num;
row = GeoCondition.row;
col = GeoCondition.col;

figure('Name', '各管流量分布', 'NumberTitle', 'off');

% 流量柱状图
bar(1:Tube_num, mdot_R * 1000, 'FaceColor', [0.3 0.6 0.9], 'EdgeColor', 'k', 'LineWidth', 0.8);
hold on;

% 平均流量线
mdot_avg = mean(mdot_R) * 1000;
yline(mdot_avg, 'r--', sprintf('均值 %.2f g/s', mdot_avg), 'LineWidth', 1.2);

xlabel('管号');
ylabel('质量流量 (g/s)');
title(sprintf('各管流量分配 (%d排×%d列, 共%d管)', row, col, Tube_num));
grid on;

% 每根管上标注数值
for t = 1:Tube_num
    text(t, mdot_R(t)*1000 + max(mdot_R)*15, ...
        sprintf('%.2f', mdot_R(t)*1000), ...
        'HorizontalAlignment', 'center', 'FontSize', 7, 'Color', [0.2 0.2 0.2]);
end
hold off;

%% 4. 收敛历史
% 修剪未使用的预分配
residual_history(i1+1:end) = [];
dp_loop_history(i1+1:end)   = [];

figure('Name', '收敛历史', 'NumberTitle', 'off');

subplot(1, 2, 1);
semilogy(1:i1, residual_history, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('迭代次数'); ylabel('能量残差');
title(sprintf('能量残差收敛历史 (最终: %.2e)', residual_history(end)));
grid on;

subplot(1, 2, 2);
if isempty(N)
    text(0.5, 0.5, '无环路，压降平衡不适用', ...
        'HorizontalAlignment', 'center', 'FontSize', 12);
    title('环路压降平衡 (N/A)');
else
    semilogy(1:i1, dp_loop_history*1e6, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 8);
    xlabel('迭代次数'); ylabel('环路压降残差 (Pa)');
    title(sprintf('环路压降收敛历史 (最终: %.2e Pa)', dp_loop_history(end)));
    grid on;
end

sgtitle(sprintf('收敛历史 (共 %d 次迭代, %.2f s)', i1, time));

%% 4.5 收敛轨迹：每次扫描/每次迭代与最终解的相对偏差
figure('Name', '收敛轨迹', 'NumberTitle', 'off');

% --- 子图1：管扫描级误差 (tube_cal 为横轴) ---
subplot(1, 2, 1);
snapshot_h_out = snapshot_h_out(1:tube_cal);
snapshot_p_out = snapshot_p_out(1:tube_cal);
snapshot_mdot  = snapshot_mdot(1:tube_cal);
snapshot_dp    = snapshot_dp(1:tube_cal);
snapshot_flag  = snapshot_flag(1:tube_cal);

err_h  = zeros(tube_cal, 1);
err_p  = zeros(tube_cal, 1);
err_m  = zeros(tube_cal, 1);
err_dp = zeros(tube_cal, 1);

for j = 1:tube_cal
    err_h(j)  = norm(snapshot_h_out{j} - h_R_out, 'fro') / norm(h_R_out, 'fro');
    err_p(j)  = norm(snapshot_p_out{j} - p_R_out, 'fro') / norm(p_R_out, 'fro');
    err_m(j)  = norm(snapshot_mdot{j}  - mdot_R) / norm(mdot_R);
    err_dp(j) = norm(snapshot_dp{j}    - dp_tube) / norm(dp_tube);
end

err_h  = max(err_h,  eps);
err_p  = max(err_p,  eps);
err_m  = max(err_m,  eps);
err_dp = max(err_dp, eps);

semilogy(1:tube_cal, err_h,  'b-o', 'LineWidth', 1.2, 'MarkerSize', 6); hold on;
semilogy(1:tube_cal, err_p,  'r-s', 'LineWidth', 1.2, 'MarkerSize', 6);
semilogy(1:tube_cal, err_m,  'g-^', 'LineWidth', 1.2, 'MarkerSize', 6);
semilogy(1:tube_cal, err_dp, 'm-d', 'LineWidth', 1.2, 'MarkerSize', 6);

% 标注扫描类型（换热/压力）
heat_idx = find(snapshot_flag == 1);
pres_idx = find(snapshot_flag == 2);
for k = 1:length(heat_idx)
    xline(heat_idx(k), '--k', 'Alpha', 0.3);
end
for k = 1:length(pres_idx)
    xline(pres_idx(k), '-.', 'Color', [0.7 0.7 0.7], 'Alpha', 0.3);
end

xlabel('管扫描累计次数 (tube\_cal)');
ylabel('与最终解的相对偏差');
title(sprintf('管扫描级收敛轨迹 (%d 次扫描)', tube_cal));
legend({'\epsilon_h', '\epsilon_p', '\epsilon_m', '\epsilon_{\Delta p}'}, ...
    'Interpreter', 'tex', 'Location', 'northeast');
grid on; hold off;

% --- 子图2：迭代级误差 (i1 为横轴) ---
subplot(1, 2, 2);
iter_h_out = iter_h_out(1:i1);
iter_p_out = iter_p_out(1:i1);
iter_mdot  = iter_mdot(1:i1);
iter_dp    = iter_dp(1:i1);

iter_err_h  = zeros(i1, 1);
iter_err_p  = zeros(i1, 1);
iter_err_m  = zeros(i1, 1);
iter_err_dp = zeros(i1, 1);

for k = 1:i1
    iter_err_h(k)  = norm(iter_h_out{k} - h_R_out, 'fro') / norm(h_R_out, 'fro');
    iter_err_p(k)  = norm(iter_p_out{k} - p_R_out, 'fro') / norm(p_R_out, 'fro');
    iter_err_m(k)  = norm(iter_mdot{k}  - mdot_R) / norm(mdot_R);
    iter_err_dp(k) = norm(iter_dp{k}    - dp_tube) / norm(dp_tube);
end

iter_err_h  = max(iter_err_h,  eps);
iter_err_p  = max(iter_err_p,  eps);
iter_err_m  = max(iter_err_m,  eps);
iter_err_dp = max(iter_err_dp, eps);

semilogy(1:i1, iter_err_h,  'b-o', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
semilogy(1:i1, iter_err_p,  'r-s', 'LineWidth', 1.5, 'MarkerSize', 8);
semilogy(1:i1, iter_err_m,  'g-^', 'LineWidth', 1.5, 'MarkerSize', 8);
semilogy(1:i1, iter_err_dp, 'm-d', 'LineWidth', 1.5, 'MarkerSize', 8);

xlabel('外层迭代次数');
ylabel('与最终解的相对偏差');
title(sprintf('迭代级收敛轨迹 (%d 次迭代)', i1));
legend({'\epsilon_h', '\epsilon_p', '\epsilon_m', '\epsilon_{\Delta p}'}, ...
    'Interpreter', 'tex', 'Location', 'northeast');
grid on; hold off;

sgtitle(sprintf('收敛轨迹 (求解耗时 %.2f s)', time));

%% 4.6 综合收敛指标 — 管扫描级，迭代分段着色
snapshot_iter = snapshot_iter(1:tube_cal);

err_total = sqrt(0.5*err_h.^2 + 0.2*err_p.^2 + 0.1*err_m.^2 + 0.2*err_dp.^2);

figure('Name', '综合收敛指标', 'NumberTitle', 'off');

% 迭代背景色带
iter_edges = [1, find(diff(snapshot_iter) ~= 0)' + 1, tube_cal + 1];
band_colors = [0.94 0.94 1.0; 1.0 0.94 0.94; 0.94 1.0 0.94; 1.0 1.0 0.85];
for k = 1:length(iter_edges)-1
    xs = iter_edges(k) - 0.5;
    xe = iter_edges(k+1) - 0.5;
    patch([xs xe xe xs], [1e-12 1e-12 1 1], ...
        band_colors(mod(k-1,4)+1,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    hold on;
end

% 综合指标曲线
semilogy(1:tube_cal, err_total, 'k.-', 'LineWidth', 2.5, 'MarkerSize', 18);

% 标注扫描类型
heat_idx = find(snapshot_flag == 1);
pres_idx = find(snapshot_flag == 2);
h_heat = scatter(heat_idx, err_total(heat_idx), 80, 'b', 'filled', 'o');
h_pres = scatter(pres_idx, err_total(pres_idx), 80, 'r', 'filled', 's');

% 迭代间分隔线
for k = 2:length(iter_edges)-1
    xline(iter_edges(k) - 0.5, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
end

xlabel('管扫描累计次数 (tube\_cal)');
ylabel('综合相对偏差 \epsilon_{total}', 'Interpreter', 'tex');
title(sprintf('综合收敛轨迹 (%d 次扫描, %d 次迭代)', tube_cal, i1));
h_total = semilogy(nan, nan, 'k.-', 'LineWidth', 2.5, 'MarkerSize', 18);
legend([h_heat, h_pres, h_total], ...
    {'换热扫描', '压力扫描', '\epsilon_{total}'}, ...
    'Interpreter', 'tex', 'Location', 'northeast');
xlim([0.5, tube_cal + 0.5]);
grid on; hold off;

%% 5. 导出标准化性能数据
hxPerf = exportHxPerf(TCinf, GeoCondition, BDCondition, ...
    h_R_in, h_R_out, p_R_in, p_R_out, ...
    T_MA_in, T_MA_out, mdot_R, mdot_MA, dp_tube, ...
    N, heatPaths, pdropPaths, residual_max, dp_loop_max, time, Prop_handle);

fprintf('\n后处理完成。工作区变量: hxPerf (标准化性能数据)\n');
