% PostProcessing — 后处理集成入口 (Demo1.02)
% 仿真完成后运行此脚本，一键生成全部可视化与性能报告

% 添加后处理函数路径
addpath(fullfile(fileparts(mfilename('fullpath')), 'PostProcessing'));

%% 1. 流路拓扑可视化
%plotTubeLayout(TCinf, GeoCondition);

%% 2. 沿线物性分布曲线
plotAlongPath(heatPaths, h_R_in, h_R_out, p_R_in, p_R_out, ...
    T_MA_in, T_MA_out, dp_tube, mdot_R, Prop_handle, row);

%% 3. 环路性能对比（有环路时生效）
plotLoopBalance(N, R_flow, mdot0, mdot_R, dp_tube, u0, options);

%% 4. 整体性能汇总表
summaryTable(TCinf, GeoCondition, BDCondition, h_R_in, h_R_out, ...
    p_R_in, p_R_out, T_MA_in, T_MA_out, mdot_R, dp_tube, heatPaths, time, Prop_handle, N);

%% 5. 收敛历史
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
    semilogy(1:i1, dp_loop_history, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 8);
    xlabel('迭代次数'); ylabel('环路压降残差 (Pa)');
    title(sprintf('环路压降收敛历史 (最终: %.2e Pa)', dp_loop_history(end)));
    grid on;
end

sgtitle(sprintf('收敛历史 (共 %d 次迭代, %.2f s)', i1, time));

%% 6. 导出标准化性能数据
hxPerf = exportHxPerf(TCinf, GeoCondition, BDCondition, ...
    h_R_in, h_R_out, p_R_in, p_R_out, ...
    T_MA_in, T_MA_out, mdot_R, mdot_MA, dp_tube, ...
    N, heatPaths, pdropPaths, residual_max, dp_loop_max, time, Prop_handle);

fprintf('\n后处理完成。工作区变量: hxPerf (标准化性能数据)\n');
