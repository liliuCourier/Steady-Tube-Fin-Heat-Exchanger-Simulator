function plotLoopBalance(N, R_flow, mdot0, mdot_R, dp_tube, u0, options)
% 环路性能对比：压降平衡验证
% 仅在存在环路 (N 非空) 时有输出

if isempty(N)
    disp('无环路，跳过环路性能对比。');
    return;
end

n_loops = size(N, 2);

figure('Name', '环路性能对比', 'NumberTitle', 'off');

% ---- 子图1: 各环路两支路压降对比 ----
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
    y_max = max(bar_data(i, :)) * 1e6;
    dy = diff(ylim);
    text(i, y_max + dy * 0.05, ...
        sprintf('Delta=%.2fPa', imbalance), ...
        'HorizontalAlignment', 'center', 'FontSize', 8);
end

% ---- 子图2: 各环路压降不平衡量 ----
subplot(1, 2, 2);
imbalance = abs(dp_loop_pos - dp_loop_neg) * 1e6;  % Pa
b2 = bar(1:n_loops, imbalance, 0.5, 'FaceColor', [0.95 0.6 0.2]);
xlabel('环路'); ylabel('压降不平衡量 (Pa)');
title('各环路压降不平衡量 (越小越平衡)');
set(gca, 'XTickLabel', loop_names);
grid on;
for i = 1:n_loops
    text(i, imbalance(i) + max(imbalance)*0.05, sprintf('%.2f', imbalance(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
end

sgtitle(sprintf('环路性能分析 (%d个环路)', n_loops));
end
