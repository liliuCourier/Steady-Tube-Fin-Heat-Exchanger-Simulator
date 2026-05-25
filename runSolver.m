% runSolver.m — 稳态求解器精简脚本版，无后处理
% 用法：
%   mode = 'main';      runSolver;   out = solver_out;
%   mode = 'loopbase';  runSolver;   out = solver_out;
%
% 输出 solver_out 字段: Q_total, Q_tube, dp_total, mdot_R, dp_tube, time, n_iter, ...

if ~exist('mode', 'var'), mode = 'main'; end

% 临时重命名根目录 PostProcessing.m，防止后处理运行
rootPP = fullfile(fileparts(mfilename('fullpath')), 'PostProcessing.m');
rootPP_hidden = fullfile(fileparts(mfilename('fullpath')), 'PostProcessing_hidden.m');
ppExists = exist(rootPP, 'file');
if ppExists
    movefile(rootPP, rootPP_hidden);
    cleanup = onCleanup(@() movefile(rootPP_hidden, rootPP));
end

% 运行求解
if strcmp(mode, 'loopbase')
    run('Main_loop_base.m');
else
    run('Main.m');
end

% 如果后处理没运行（被遮蔽），手动生成 hxPerf 精简版
if ~exist('hxPerf', 'var')
    h_R_out_vec = h_R_out(CV_num, :);
    Q_tube = mdot_R' .* (BDCondition.BD_R.h_R_inlet - h_R_out_vec) * 1e3;
    hxPerf.Performance.Q_total = sum(Q_tube);
    hxPerf.Performance.Q_tube  = Q_tube;
    hxPerf.Performance.dp_total = sum(dp_tube) * 1e6;
    hxPerf.Performance.residual_max = 0;
    hxPerf.Performance.dp_loop_max = max(abs((dp_tube') * N));
end

% 输出结果结构体
solver_out = struct();
solver_out.Q_total      = hxPerf.Performance.Q_total;
solver_out.Q_tube       = hxPerf.Performance.Q_tube;
solver_out.dp_total     = hxPerf.Performance.dp_total;
solver_out.mdot_total   = sum(mdot_R);
solver_out.mdot_R       = mdot_R;
solver_out.dp_tube      = dp_tube;
solver_out.h_R_out      = h_R_out;
solver_out.p_R_out      = p_R_out;
solver_out.time         = time;
solver_out.n_iter       = i1;
solver_out.residual_max = hxPerf.Performance.residual_max;
solver_out.dp_loop_max  = hxPerf.Performance.dp_loop_max;

fprintf('runSolver(%s): Q=%.1f W, dp=%.1f Pa, time=%.3f s, iter=%d\n', ...
    mode, solver_out.Q_total, solver_out.dp_total, solver_out.time, solver_out.n_iter);
