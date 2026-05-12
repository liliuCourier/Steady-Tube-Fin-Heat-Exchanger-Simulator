function plotAlongPath(heatPaths, h_R_in, h_R_out, p_R_in, p_R_out, ...
    T_MA_in, T_MA_out, dp_tube, mdot_R, Prop_handle, row)
% 沿广度优先路径绘制关键物性分布曲线

Tube_num = length(heatPaths);
path_order = heatPaths;

% 沿管序提取出口值
T_R_out = zeros(1, Tube_num);
p_R_mid = zeros(1, Tube_num);
q_tube = zeros(1, Tube_num);
T_MA_mid = zeros(1, Tube_num);

for idx = 1:Tube_num
    t = path_order(idx);
    % 工质出口温度（用出口焓和出口压力反推）
    h_out = h_R_out(end, t);
    p_out = p_R_out(end, t);  
    [~, ~, ~, T_R_out(idx)] = Prop1(p_out, h_out, Prop_handle);
    % 平均压力
    p_R_mid(idx) = (p_R_in(1, t) + p_R_out(end, t)) / 2;
    % 单管换热量
    q_tube(idx) = abs((h_R_in(1, t) - h_R_out(end, t)) * mdot_R(t) * 1000);  % W
    % 空气平均温度
    T_MA_mid(idx) = (T_MA_in(1, t) + T_MA_out(end, t)) / 2;
end

figure('Name', '沿线物性分布', 'NumberTitle', 'off');

% ---- 子图1: 温度分布 ----
subplot(2, 2, 1);
plot(1:Tube_num, T_R_out, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 8);
hold on;
plot(1:Tube_num, T_MA_mid, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('管序 (BFS)'); ylabel('温度 (K)');
title('温度沿线分布');
legend({'工质出口温度', '空气平均温度'}, 'Location', 'best');
xticks(1:Tube_num);
xticklabels(cellstr(num2str(path_order')));
grid on;

% ---- 子图2: 压力/压降分布 ----
subplot(2, 2, 2);
yyaxis left;
plot(1:Tube_num, p_R_mid, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 8);
ylabel('平均压力 (MPa)');
yyaxis right;
bar(1:Tube_num, dp_tube(path_order) * 1e6, 0.4, 'FaceColor', [0.9 0.5 0.2]);
ylabel('管压降 (Pa)');
xlabel('管序 (BFS)');
title('压力与压降沿线分布');
xticks(1:Tube_num);
xticklabels(cellstr(num2str(path_order')));
grid on;

% ---- 子图3: 换热量分布 ----
subplot(2, 2, 3);
bar(1:Tube_num, q_tube, 0.6, 'FaceColor', [0.2 0.6 0.8]);
hold on;
q_total = sum(q_tube);
for idx = 1:Tube_num
    pct = q_tube(idx) / q_total * 100;
    text(idx, q_tube(idx) + max(q_tube)*0.02, sprintf('%.1f%%', pct), ...
        'HorizontalAlignment', 'center', 'FontSize', 8);
end
xlabel('管序 (BFS)'); ylabel('换热量 (W)');
title(sprintf('各管换热量 (总计: %.1f W)', q_total));
xticks(1:Tube_num);
xticklabels(cellstr(num2str(path_order')));
grid on;

% ---- 子图4: 压降累积曲线 ----
subplot(2, 2, 4);
dp_cum = cumsum(dp_tube(path_order)) * 1e6;  % Pa (dp_tube 为 MPa)
plot(1:Tube_num, dp_cum, 'k-o', 'LineWidth', 2, 'MarkerSize', 8, ...
    'MarkerFaceColor', 'k');
xlabel('管序 (BFS)'); ylabel('累积压降 (Pa)');
title(sprintf('压降累积曲线 (总计: %.1f Pa)', dp_cum(end)));
xticks(1:Tube_num);
xticklabels(cellstr(num2str(path_order')));
grid on;

sgtitle(sprintf('沿线物性分布 (%d管, %d排)', Tube_num, row));
end
