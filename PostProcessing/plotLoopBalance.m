function plotLoopBalance(N, R_flow, mdot0, mdot_R, dp_tube, u0, options)
% 环路性能对比：压降平衡验证
% 仅在存在环路 (N 非空) 时有输出
%
% 输入：
%   N        - 零空间基向量 (Tube_num × n_loops)
%   R_flow   - 各管流阻 (Tube_num×1)
%   mdot0    - 最小范数流量解 (Tube_num×1)
%   mdot_R   - 当前流量分配 (Tube_num×1)
%   dp_tube  - 各管压降 (Tube_num×1)
%   u0       - 自由变量初值
%   options  - fsolve 选项

if isempty(N)
    disp('无环路，跳过环路性能对比。');
    return;
end

n_loops = size(N, 2);

figure('Name', '环路性能对比', 'NumberTitle', 'off');

% ---- 子图1: 各环路压降平衡 ----
subplot(1, 2, 1);
dp_loop_pos = zeros(1, n_loops);
dp_loop_neg = zeros(1, n_loops);
loop_names = cell(1, n_loops);

for i = 1:n_loops
    tubes_pos = find(N(:, i) == 1)';
    tubes_neg = find(N(:, i) == -1)';
    dp_loop_pos(i) = sum(dp_tube(tubes_pos));
    dp_loop_neg(i) = sum(dp_tube(tubes_neg));
    loop_names{i} = sprintf('L%d', i);
end

bar_data = [dp_loop_pos; dp_loop_neg]';
b = bar(1:n_loops, bar_data * 1e6, 'grouped');
b(1).FaceColor = [0.2 0.6 0.8];
b(2).FaceColor = [0.8 0.3 0.2];
set(gca, 'XTickLabel', loop_names);
xlabel('环路'); ylabel('支路压降 (Pa)');
title('各环路两支路压降对比');
legend({'正向支路', '反向支路'}, 'Location', 'best');
grid on;

% 标注不平衡量
for i = 1:n_loops
    imbalance = abs(dp_loop_pos(i) - dp_loop_neg(i)) * 1e6;
    y_max = max(bar_data(i, :));
    text(i, y_max + range(ylim)*0.05, ...
        sprintf('Δ=%.2fPa', imbalance), ...
        'HorizontalAlignment', 'center', 'FontSize', 8);
end

% ---- 子图2: 自由变量收敛轨迹 ----
subplot(1, 2, 2);
if ~isempty(u0)
    u0_norm = u0 / (norm(u0) + eps);
    u_scan = linspace(-2, 2, 50) * norm(u0);
    F_scan = zeros(n_loops, length(u_scan));
    for j = 1:length(u_scan)
        F_scan(:, j) = uF_demo(u_scan(j) * u0_norm, R_flow, N, mdot0);
    end
    for i = 1:n_loops
        plot(u_scan, F_scan(i, :), 'LineWidth', 1.5);
        hold on;
    end
    xline(0, 'k--', 'LineWidth', 1);
    yline(0, 'k--', 'LineWidth', 1);
    xlabel('自由变量 u'); ylabel('环路压降残差 F(u)');
    title('环路残差函数');
    legend(cellstr(num2str((1:n_loops)', '环路 %d')), 'Location', 'best');
    grid on;
end

sgtitle(sprintf('环路性能分析 (%d个环路)', n_loops));
end

function F = uF_demo(u_vec, R_flow, N, mdot0)
    dp_t = (R_flow) .* (mdot0 + N * u_vec).^2;
    F = (dp_t') * N;
end
