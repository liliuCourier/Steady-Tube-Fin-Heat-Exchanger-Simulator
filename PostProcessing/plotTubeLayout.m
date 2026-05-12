function plotTubeLayout(TCinf, GeoCondition)
% 流路拓扑可视化：管排网格 + 流向箭头 + 进出口标记
%
% 输入：
%   TCinf        - 管路连接信息结构体 (含 TC_matrix, FlowDirection, row, col, inlet_num, outlet_num)
%   GeoCondition - 几何信息结构体 (含 row, col)

row = GeoCondition.row;
col = GeoCondition.col;
Tube_num = GeoCondition.Tube_num;
TC_matrix = TCinf.TC_matrix;
FlowDirection = TCinf.FlowDirection;
inlet_num = TCinf.inlet_num;
outlet_num = TCinf.outlet_num;

% 管排网格坐标：每根管的 (x, y) 位置
% 列沿 x 方向，排（行）沿 y 方向
x = zeros(1, Tube_num);
y = zeros(1, Tube_num);
for t = 1:Tube_num
    c = ceil(t / row);          % 第几列
    r = mod(t - 1, row) + 1;    % 第几排 (从上到下)
    x(t) = c;
    y(t) = row - r + 1;         % 翻转 y 使第一排在上
end

figure('Name', '流路拓扑可视化', 'NumberTitle', 'off');
hold on;

% 绘制管路连接线
for i = 1:size(TC_matrix, 1)
    node_tubes = find(TC_matrix(i, :) ~= 0);
    for j = 1:length(node_tubes)
        for k = j+1:length(node_tubes)
            t1 = node_tubes(j);
            t2 = node_tubes(k);
            plot([x(t1) x(t2)], [y(t1) y(t2)], 'k-', 'LineWidth', 1.5);
        end
    end
end

% 绘制管路圆点 + 流向箭头
for t = 1:Tube_num
    % 管路圆点
    plot(x(t), y(t), 'o', 'MarkerSize', 18, ...
        'MarkerFaceColor', [0.8 0.9 1], 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

    % 管号标注
    text(x(t), y(t), num2str(t), 'HorizontalAlignment', 'center', ...
        'FontSize', 10, 'FontWeight', 'bold');

    % 流向箭头
    arrow_len = 0.2;
    if FlowDirection(t) == 1
        quiver(x(t) - arrow_len, y(t), arrow_len * 2, 0, 0, ...
            'MaxHeadSize', 0.5, 'Color', [0 0.4 0.8], 'LineWidth', 1.5);
    else
        quiver(x(t) + arrow_len, y(t), -arrow_len * 2, 0, 0, ...
            'MaxHeadSize', 0.5, 'Color', [0.8 0.2 0], 'LineWidth', 1.5);
    end
end

% 进出口高亮
plot(x(inlet_num), y(inlet_num), 's', 'MarkerSize', 22, ...
    'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
text(x(inlet_num), y(inlet_num) + 0.3, 'IN', 'HorizontalAlignment', 'center', ...
    'FontSize', 10, 'FontWeight', 'bold', 'Color', 'g');

for o = 1:length(outlet_num)
    ot = outlet_num(o);
    plot(x(ot), y(ot), '^', 'MarkerSize', 22, ...
        'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
    text(x(ot), y(ot) - 0.3, 'OUT', 'HorizontalAlignment', 'center', ...
        'FontSize', 10, 'FontWeight', 'bold', 'Color', 'r');
end

% 图例
plot(nan, nan, 'o', 'MarkerSize', 12, 'MarkerFaceColor', [0.8 0.9 1], 'MarkerEdgeColor', 'k');  % dummy for legend
plot(nan, nan, 's', 'MarkerSize', 12, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k');
plot(nan, nan, '^', 'MarkerSize', 12, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
legend({'管路', '进口', '出口'}, 'Location', 'bestoutside');

axis equal;
xlim([min(x)-0.8, max(x)+0.8]);
ylim([min(y)-0.8, max(y)+0.8]);
title(sprintf('流路拓扑 (%d排 × %d列, %d管)', row, col, Tube_num));
xlabel('列'); ylabel('排');
set(gca, 'XTick', 1:col, 'YTick', 1:row);
grid on;
hold off;

fprintf('流路拓扑图已生成 (%d管, %d排×%d列)\n', Tube_num, row, col);
end
