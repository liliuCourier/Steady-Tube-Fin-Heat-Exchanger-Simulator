% Main_loop_base — 环路优先求解器 (Demo1.12 同步)
% Phase1: 仅热力扫描 + 提取压阻 + 流量更新（不更新压力场），迭代至流量稳定
% Phase2: 完整热力+压力+流量重分配，与 Main 一致
% 运行前请先执行 PreProcessing

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'Lib'));
addpath(fullfile(root, 'Solver'));
addpath(fullfile(root, 'PostProcessing'));
addpath(fullfile(root, 'PreProc'));

if ~exist('TCinf','var'),       error('TCinf 未找到，请先运行 PreProcessing'); end
if ~exist('GeoCondition','var'), error('GeoCondition 未找到，请先运行 PreProcessing'); end
if ~exist('BDCondition','var'),  error('BDCondition 未找到，请先运行 PreProcessing'); end

% options = optimoptions('fsolve','Display','none',...
%     'Algorithm','levenberg-marquardt',...
%     'FunctionTolerance',1e-12,...
%     'MaxFunctionEvaluations',5e4,...
%     'StepTolerance',1e-8,...
%     'UseParallel',false,...
%     'ScaleProblem','jacobian');

%% 初始化
[N,u0,mdot0,mdot_R_init,exitflag] = mdot_Initial(BDCondition,GeoCondition,TCinf);
[heatPaths,pdropPaths,predecessors_in,predecessors_out] = buildPath(TCinf.TC_matrix,N);

CV_num   = BDCondition.CV_num;
row      = GeoCondition.row;
Tube_num = GeoCondition.Tube_num;

h_R_inlet = BDCondition.BD_R.h_R_inlet;
p_R_inlet = BDCondition.BD_R.p_R_inlet;
mdot_MA_inlet = BDCondition.BD_MA.mdot_MA_inlet;
T_MA_inlet    = BDCondition.BD_MA.T_MA_inlet;
p_MA_inlet    = BDCondition.BD_MA.p_MA_inlet;

p_R_in  = p_R_inlet*ones(CV_num,Tube_num);
p_R_out = p_R_inlet*ones(CV_num,Tube_num);
h_R_in  = h_R_inlet*ones(CV_num,Tube_num);
h_R_out = h_R_inlet*ones(CV_num,Tube_num);
mdot_R  = mdot_R_init;

T_MA_in  = T_MA_inlet*ones(CV_num,Tube_num);
T_MA_out = T_MA_inlet*ones(CV_num,Tube_num);
p_MA_in  = p_MA_inlet*ones(CV_num,Tube_num);
p_MA_out = p_MA_inlet*ones(CV_num,Tube_num);
mdot_MA  = mdot_MA_inlet/row/CV_num*ones(CV_num,Tube_num);

dp_tube = zeros(Tube_num,1);
R_flow  = [];
R_coef  = 1.81;

loopmax        = 10;
residual_limit = 1e-3;

residual_history = zeros(loopmax, 1);
dp_loop_history  = zeros(loopmax, 1);
tube_cal = 0;

max_scans = loopmax * 3 + 1;
snapshot_h_out = cell(max_scans, 1);
snapshot_p_out = cell(max_scans, 1);
snapshot_mdot  = cell(max_scans, 1);
snapshot_dp    = cell(max_scans, 1);
snapshot_flag  = zeros(max_scans, 1);
snapshot_iter  = zeros(max_scans, 1);

iter_h_out = cell(loopmax, 1);
iter_p_out = cell(loopmax, 1);
iter_mdot  = cell(loopmax, 1);
iter_dp    = cell(loopmax, 1);

% 流量分配演进记录
u0_history = cell(loopmax+1, 1);
u0_history{1} = u0;

has_loops = ~isempty(N) && size(N,2) > 0;

%% Phase 1: 仅热力扫描 + 提取压阻，迭代至流量稳定
% 压力场全程不变（p_R_inlet），11 项饱和物性所有 CV/迭代共用
sat_global = cell(1,11);
sat_global{1} = Prop_handle.v_liq(0, p_R_inlet);   sat_global{2} = Prop_handle.v_vap(1, p_R_inlet);
sat_global{3} = Prop_handle.Pr_liq(0, p_R_inlet);  sat_global{4} = Prop_handle.Pr_vap(1, p_R_inlet);
sat_global{5} = Prop_handle.Nu_liq(0, p_R_inlet);  sat_global{6} = Prop_handle.Nu_vap(1, p_R_inlet);
sat_global{7} = Prop_handle.k_liq(0, p_R_inlet);   sat_global{8} = Prop_handle.k_vap(1, p_R_inlet);
sat_global{9} = Prop_handle.T_liq(0, p_R_inlet);
sat_global{10} = Prop_handle.h_sat_liq(p_R_inlet); sat_global{11} = Prop_handle.h_sat_vap(p_R_inlet);
tic
if has_loops
    for i1 = 1:loopmax
        iter_h_out{i1} = h_R_out;
        iter_p_out{i1} = p_R_out;
        iter_mdot{i1}  = mdot_R;
        iter_dp{i1}    = dp_tube;

        % 热力扫描（提取各 CV 压阻，压力场不变）
        residual_max = 0;
        [h_R_in, h_R_out, T_MA_in, T_MA_out, ...
         p_R_in, p_R_out, p_MA_in, p_MA_out, ...
         R_flow, dp_tube, residual_max] = scanTubes_phase1(...
            heatPaths, predecessors_in, ...
            h_R_in, h_R_out, T_MA_in, T_MA_out, ...
            p_R_in, p_R_out, p_MA_in, p_MA_out, ...
            mdot_R, mdot_MA, TCinf, GeoCondition, ...
            CV_num, row, Prop_handle, ...
            h_R_inlet, p_R_inlet, T_MA_inlet, p_MA_inlet, ...
            residual_max, R_coef, sat_global);
        tube_cal = tube_cal + 1;
        snapshot_h_out{tube_cal} = h_R_out;
        snapshot_p_out{tube_cal} = p_R_out;
        snapshot_mdot{tube_cal}  = mdot_R;
        snapshot_dp{tube_cal}    = dp_tube;
        snapshot_flag(tube_cal)  = 1;
        snapshot_iter(tube_cal)  = i1;

        % 流量更新 — Newton 迭代
        mdot_R_prev = mdot_R;
        dp_Pa = dp_tube*1e6;
        R_flow_coef = dp_Pa ./ (mdot_R.^R_coef);
        u = u0;
        for k = 1:10
            m_k   = mdot0 + N * u;
            F     = N' * (R_flow_coef .* m_k.^R_coef);
            S     = R_coef * R_flow_coef .* m_k.^(R_coef - 1);
            J     = N' * (S .* N);
            du    = -J \ F;
            u     = u + du;
            if norm(du) < 1e-8, break; end
        end
        mdot_R = mdot0 + N * u;
        u0 = u;
        u0_history{i1+1} = u0;

        residual_history(i1) = residual_max;
        dp_loop_history(i1) = max(abs((dp_tube') * N));

        % Phase1 步出判据：前后两次流量相对变化 < 1e-3
        delta_mdot = max(abs(mdot_R - mdot_R_prev) ./ mdot_R);
        if i1 > 1 && delta_mdot < 1e-2
            fprintf('Phase1 流量稳定 (delta_mdot=%.2e)，进入 Phase2\n', delta_mdot);
            break
        end
    end
    i1_phase1 = i1;
else
    i1_phase1 = 0;
end

%% Phase 2: 完整扫描 — 与 Main 一致（热力+压力+流量更新）
if has_loops
    for i2 = i1_phase1+1:loopmax
        iter_h_out{i2} = h_R_out;
        iter_p_out{i2} = p_R_out;
        iter_mdot{i2}  = mdot_R;
        iter_dp{i2}    = dp_tube;

        % 热力扫描
        residual_max = 0;
        [h_R_in, h_R_out, T_MA_in, T_MA_out, ...
         p_R_in, p_R_out, p_MA_in, p_MA_out, ...
         dp_tube, residual_max] = scanTubes(...
            heatPaths, predecessors_in, ...
            h_R_in, h_R_out, T_MA_in, T_MA_out, ...
            p_R_in, p_R_out, p_MA_in, p_MA_out, ...
            mdot_R, mdot_MA, TCinf, GeoCondition, ...
            CV_num, row, 2, Prop_handle, ...
            h_R_inlet, p_R_inlet, T_MA_inlet, p_MA_inlet, ...
            residual_max, dp_tube);
        tube_cal = tube_cal + 1;
        snapshot_h_out{tube_cal} = h_R_out;
        snapshot_p_out{tube_cal} = p_R_out;
        snapshot_mdot{tube_cal}  = mdot_R;
        snapshot_dp{tube_cal}    = dp_tube;
        snapshot_flag(tube_cal)  = 2;
        snapshot_iter(tube_cal)  = i2;

        % 压力扫描
        [h_R_in, h_R_out, T_MA_in, T_MA_out, ...
         p_R_in, p_R_out, p_MA_in, p_MA_out, ...
         dp_tube, residual_max] = scanTubes(...
            heatPaths, predecessors_in, ...
            h_R_in, h_R_out, T_MA_in, T_MA_out, ...
            p_R_in, p_R_out, p_MA_in, p_MA_out, ...
            mdot_R, mdot_MA, TCinf, GeoCondition, ...
            CV_num, row, 1, Prop_handle, ...
            h_R_inlet, p_R_inlet, T_MA_inlet, p_MA_inlet, ...
            residual_max, dp_tube);
        tube_cal = tube_cal + 1;
        snapshot_h_out{tube_cal} = h_R_out;
        snapshot_p_out{tube_cal} = p_R_out;
        snapshot_mdot{tube_cal}  = mdot_R;
        snapshot_dp{tube_cal}    = dp_tube;
        snapshot_flag(tube_cal)  = 1;
        snapshot_iter(tube_cal)  = i2;

        % 流量更新 — Newton 迭代
        dp_Pa = dp_tube * 1e6;
        R_flow_coef = dp_Pa ./ (mdot_R.^R_coef);
        u = u0;
        for k = 1:10
            m_k   = mdot0 + N * u;
            F     = N' * (R_flow_coef .* m_k.^R_coef);
            S     = R_coef * R_flow_coef .* m_k.^(R_coef - 1);
            J     = N' * (S .* N);
            du    = -J \ F;
            u     = u + du;
            if norm(du) < 1e-8, break; end
        end
        mdot_R = mdot0 + N * u;
        u0 = u;
        u0_history{i2+1} = u0;

        residual_history(i2) = residual_max;
        dp_loop_history(i2)  = max(abs((dp_tube')*N));

        if max(abs((dp_tube')*N)) < 1e-6 && (residual_max < residual_limit)
            disp("压降收敛")
            i1 = i2;
            break
        end
    end
else
    % 无环路：热力+压力，残差收敛即步出
    for i1 = 1:loopmax
        iter_h_out{i1} = h_R_out;
        iter_p_out{i1} = p_R_out;
        iter_mdot{i1}  = mdot_R;
        iter_dp{i1}    = dp_tube;

        residual_max = 0;
        [h_R_in, h_R_out, T_MA_in, T_MA_out, ...
         p_R_in, p_R_out, p_MA_in, p_MA_out, ...
         dp_tube, residual_max] = scanTubes(...
            heatPaths, predecessors_in, ...
            h_R_in, h_R_out, T_MA_in, T_MA_out, ...
            p_R_in, p_R_out, p_MA_in, p_MA_out, ...
            mdot_R, mdot_MA, TCinf, GeoCondition, ...
            CV_num, row, 1, Prop_handle, ...
            h_R_inlet, p_R_inlet, T_MA_inlet, p_MA_inlet, ...
            residual_max, dp_tube);
        tube_cal = tube_cal + 1;
        snapshot_h_out{tube_cal} = h_R_out;
        snapshot_p_out{tube_cal} = p_R_out;
        snapshot_mdot{tube_cal}  = mdot_R;
        snapshot_dp{tube_cal}    = dp_tube;
        snapshot_flag(tube_cal)  = 1;
        snapshot_iter(tube_cal)  = i1;

        [h_R_in, h_R_out, T_MA_in, T_MA_out, ...
         p_R_in, p_R_out, p_MA_in, p_MA_out, ...
         dp_tube, residual_max] = scanTubes(...
            heatPaths, predecessors_in, ...
            h_R_in, h_R_out, T_MA_in, T_MA_out, ...
            p_R_in, p_R_out, p_MA_in, p_MA_out, ...
            mdot_R, mdot_MA, TCinf, GeoCondition, ...
            CV_num, row, 2, Prop_handle, ...
            h_R_inlet, p_R_inlet, T_MA_inlet, p_MA_inlet, ...
            residual_max, dp_tube);
        tube_cal = tube_cal + 1;
        snapshot_h_out{tube_cal} = h_R_out;
        snapshot_p_out{tube_cal} = p_R_out;
        snapshot_mdot{tube_cal}  = mdot_R;
        snapshot_dp{tube_cal}    = dp_tube;
        snapshot_flag(tube_cal)  = 2;
        snapshot_iter(tube_cal)  = i1;

        dp_loop_history(i1) = 0;
        if residual_max < residual_limit
            disp("残差收敛")
            break
        end
    end
end
time = toc;

dp_loop_max = dp_loop_history(i1);
fprintf('求解完成，耗时 %.2f s，正在运行后处理...\n', time);
% 
% try
%     PostProcessing;
% catch ME
%     fprintf('后处理运行出错: %s\n', ME.message);
% end

%%
% function F = uF(u,R_flow,N,mdot0,R_coef)
%  dp_tube = (R_flow).*(mdot0 + N*u).^R_coef;
%  dp_loop = (dp_tube')*N;
%  F = dp_loop;
% end



function  [F, cache_R_out, cache_MA_out] = alg(x0,BD,InletBD,GeoCondition,CV,N,solver_flag,Prop_handle,cache_R_in,cache_MA_in)

loopmax = 100;
residual_Energy = 1e-3;
residual_dp = 1e-3;

if solver_flag == 3
    for i = 1:loopmax
        x0_R  = [BD(1);x0(1)];
        BD_R.h_R_inlet = BD(2);
        BD_R.mdot_R_inlet = InletBD(1);
        BD_R.p_R_inlet = x0(2);
        BD_R.Uband_L = InletBD(2);
        flag = 2;
        out = R_cal_10(x0_R,BD_R,GeoCondition,CV,Prop_handle,flag);
        dp_R = out;
        if max(abs((dp_R - (x0(2)-x0(1))*1e6)/dp_R))<residual_dp && i<loopmax
            x0(1) = x0(2) - dp_R/1e6;
            break
        elseif i==loopmax
            disp("Uband压降不动点迭代失败")
        else
            x0(1) = x0(2) - dp_R/1e6;
        end
    end
    F = [dp_R,x0(1)];
    cache_R_out = {}; cache_MA_out = {};
    return
end

h_R_inlet = InletBD(1);
p_R_inlet = InletBD(2);
T_MA_inlet = InletBD(3);
p_MA_inlet = InletBD(4);
mdot_R_inlet  = InletBD(5);
mdot_MA_inlet = InletBD(6);

BD_R.h_R_inlet    = h_R_inlet;
BD_R.mdot_R_inlet = mdot_R_inlet;
BD_R.p_R_inlet    = p_R_inlet;
BD_MA.T_MA_inlet    = T_MA_inlet;
BD_MA.mdot_MA_inlet = mdot_MA_inlet;
BD_MA.p_MA_inlet    = p_MA_inlet;

A_R = GeoCondition.A_R;
A_MA = GeoCondition.A_MA;
Tube_num = GeoCondition.Tube_num;
A_R_CV = A_R/CV/Tube_num;
A_MA_CV = A_MA/CV/Tube_num;

if solver_flag == 1

    x0_R  = [x0(1);BD(1)];
    x0_MA = [x0(2);BD(2)];

    % 热力扫描中压力不变，饱和物性仅取决于 p，提前计算一并用于入口 Prop1 和出口 sat_cache
    sat_in = cell(1,11);
    sat_in{1} = Prop_handle.v_liq(0, p_R_inlet); sat_in{2} = Prop_handle.v_vap(1, p_R_inlet);
    sat_in{3} = Prop_handle.Pr_liq(0, p_R_inlet); sat_in{4} = Prop_handle.Pr_vap(1, p_R_inlet);
    sat_in{5} = Prop_handle.Nu_liq(0, p_R_inlet); sat_in{6} = Prop_handle.Nu_vap(1, p_R_inlet);
    sat_in{7} = Prop_handle.k_liq(0, p_R_inlet); sat_in{8} = Prop_handle.k_vap(1, p_R_inlet);
    sat_in{9} = Prop_handle.T_liq(0, p_R_inlet);
    sat_in{10} = Prop_handle.h_sat_liq(p_R_inlet); sat_in{11} = Prop_handle.h_sat_vap(p_R_inlet);

    % 入口缓存：传入 sat_in 避免 Prop1 内重复计算饱和值
    if nargin >= 9 && ~isempty(cache_R_in)
        inlet_props = cache_R_in;
    else
        inlet_props = cell(1,14);
        [inlet_props{:}] = Prop1(p_R_inlet, h_R_inlet, Prop_handle, sat_in);
    end
    if nargin >= 10 && ~isempty(cache_MA_in)
        inlet_air_props = cache_MA_in;
    else
        Ra = 287.047;
        inlet_air_props = cell(1,6);
        inlet_air_props{1} = Prop_handle.hair(T_MA_inlet);
        inlet_air_props{2} = Prop_handle.Prair(T_MA_inlet);
        inlet_air_props{3} = Prop_handle.visair(T_MA_inlet);
        inlet_air_props{4} = Prop_handle.kair(T_MA_inlet);
        inlet_air_props{5} = Prop_handle.cpair(T_MA_inlet);
        inlet_air_props{6} = Ra * T_MA_inlet / (p_MA_inlet * 1e6);
    end

    % 出口饱和缓存：复用 sat_in，仅 Tsatliq 重取 p_out
    sat_cache = cell(1,11);
    sat_cache{1} = sat_in{1}; sat_cache{2} = sat_in{2};
    sat_cache{3} = sat_in{3}; sat_cache{4} = sat_in{4};
    sat_cache{5} = sat_in{5}; sat_cache{6} = sat_in{6};
    sat_cache{7} = sat_in{7}; sat_cache{8} = sat_in{8};
    sat_cache{9} = Prop_handle.T_liq(0, BD(1));
    sat_cache{10} = sat_in{10}; sat_cache{11} = sat_in{11};

    for i = 1:loopmax
        out = R_cal_10(x0_R,BD_R,GeoCondition,CV,Prop_handle,[],inlet_props,sat_cache);

        HTC_R = out{3};
        dEF_R = out{5};
        T_R   = out{6};

        out_MA = DryA_cal_10(x0_MA,BD_MA,GeoCondition,CV,N,Prop_handle,inlet_air_props);

        dEF_MA = out_MA{1};
        n_fin  = out_MA{2};
        HTC_MA = out_MA{3};
        cp_MA  = out_MA{4};

        dT = T_R -(T_MA_inlet + x0(2))/2;
        UA_R = 1./(1./(HTC_R*A_R_CV)+1./(n_fin.*HTC_MA*A_MA_CV));
        Q = dT.*UA_R;

        if max(abs([(dEF_R*1e3-Q)/Q,(dEF_MA+Q)/Q]))<residual_Energy && i<loopmax
            x0_R(1) = h_R_inlet - Q/mdot_R_inlet/1e3;
            x0_MA(1) = Q/cp_MA/(mdot_MA_inlet) + T_MA_inlet;
            break
        elseif i==loopmax
            disp("换热不动点迭代失败")
        else
            x0_R(1) = h_R_inlet - Q/mdot_R_inlet/1e3;
            x0_MA(1) = Q/cp_MA/(mdot_MA_inlet) + T_MA_inlet;
        end
    end

    cache_R_out = out{7};
    cache_MA_out = out_MA{7};
    F = [x0_R(1);x0_MA(1)];

elseif solver_flag == 2

    x0_R  = [BD(1);x0(1)];
    x0_MA = [BD(2);x0(2)];

    sat_in = cell(1,11);
    sat_in{1} = Prop_handle.v_liq(0, p_R_inlet); sat_in{2} = Prop_handle.v_vap(1, p_R_inlet);
    sat_in{3} = Prop_handle.Pr_liq(0, p_R_inlet); sat_in{4} = Prop_handle.Pr_vap(1, p_R_inlet);
    sat_in{5} = Prop_handle.Nu_liq(0, p_R_inlet); sat_in{6} = Prop_handle.Nu_vap(1, p_R_inlet);
    sat_in{7} = Prop_handle.k_liq(0, p_R_inlet); sat_in{8} = Prop_handle.k_vap(1, p_R_inlet);
    sat_in{9} = Prop_handle.T_liq(0, p_R_inlet);
    sat_in{10} = Prop_handle.h_sat_liq(p_R_inlet); sat_in{11} = Prop_handle.h_sat_vap(p_R_inlet);

    if nargin >= 9 && ~isempty(cache_R_in)
        inlet_props = cache_R_in;
    else
        inlet_props = cell(1,14);
        [inlet_props{:}] = Prop1(p_R_inlet, h_R_inlet, Prop_handle, sat_in);
    end
    if nargin >= 10 && ~isempty(cache_MA_in)
        inlet_air_props = cache_MA_in;
    else
        Ra = 287.047;
        inlet_air_props = cell(1,6);
        inlet_air_props{1} = Prop_handle.hair(T_MA_inlet);
        inlet_air_props{2} = Prop_handle.Prair(T_MA_inlet);
        inlet_air_props{3} = Prop_handle.visair(T_MA_inlet);
        inlet_air_props{4} = Prop_handle.kair(T_MA_inlet);
        inlet_air_props{5} = Prop_handle.cpair(T_MA_inlet);
        inlet_air_props{6} = Ra * T_MA_inlet / (p_MA_inlet * 1e6);
    end

    for i = 1:loopmax
        out = R_cal_10(x0_R,BD_R,GeoCondition,CV,Prop_handle,[],inlet_props);
        dp_R = out{4};
        out_MA = DryA_cal_10(x0_MA,BD_MA,GeoCondition,CV,N,Prop_handle,inlet_air_props);
        dp_MA = out_MA{5};
        if max(abs([(dp_R-(p_R_inlet-x0_R(2))*1e6)/dp_R,...
                    (dp_MA-(p_MA_inlet-x0_MA(2))*1e6)/dp_MA]))<residual_dp && i<loopmax
            x0_R(2) = p_R_inlet - dp_R/1e6;
            x0_MA(2) = p_MA_inlet - dp_MA/1e6;
            break
        elseif i==loopmax
            disp("压降不动点迭代失败")
        else
            x0_R(2) = p_R_inlet - dp_R/1e6;
            x0_MA(2) = p_MA_inlet - dp_MA/1e6;
        end
    end
    cache_R_out = out{7};
    cache_MA_out = out_MA{7};
    F = [x0_R(2);x0_MA(2)];
end
end



function  [F, dp_CV, cache_R_out, cache_MA_out] = alg_phase1(x0,BD,InletBD,GeoCondition,CV,N,Prop_handle,cache_R_in,cache_MA_in,sat_global)
% 与 alg(solver_flag=1) 相同，额外返回 dp_CV (Pa)
% Phase1 压力场不变，sat_global 全局预计算，所有 CV 共用

loopmax = 100;
residual_Energy = 1e-3;

h_R_inlet = InletBD(1);
p_R_inlet = InletBD(2);
T_MA_inlet = InletBD(3);
p_MA_inlet = InletBD(4);
mdot_R_inlet  = InletBD(5);
mdot_MA_inlet = InletBD(6);

BD_R.h_R_inlet    = h_R_inlet;
BD_R.mdot_R_inlet = mdot_R_inlet;
BD_R.p_R_inlet    = p_R_inlet;
BD_MA.T_MA_inlet    = T_MA_inlet;
BD_MA.mdot_MA_inlet = mdot_MA_inlet;
BD_MA.p_MA_inlet    = p_MA_inlet;

A_R = GeoCondition.A_R;
A_MA = GeoCondition.A_MA;
Tube_num = GeoCondition.Tube_num;
A_R_CV = A_R/CV/Tube_num;
A_MA_CV = A_MA/CV/Tube_num;

x0_R  = [x0(1);BD(1)];
x0_MA = [x0(2);BD(2)];

if nargin >= 8 && ~isempty(cache_R_in)
    inlet_props = cache_R_in;
else
    inlet_props = cell(1,14);
    [inlet_props{:}] = Prop1(p_R_inlet, h_R_inlet, Prop_handle, sat_global);
end
if nargin >= 9 && ~isempty(cache_MA_in)
    inlet_air_props = cache_MA_in;
else
    Ra = 287.047;
    inlet_air_props = cell(1,6);
    inlet_air_props{1} = Prop_handle.hair(T_MA_inlet);
    inlet_air_props{2} = Prop_handle.Prair(T_MA_inlet);
    inlet_air_props{3} = Prop_handle.visair(T_MA_inlet);
    inlet_air_props{4} = Prop_handle.kair(T_MA_inlet);
    inlet_air_props{5} = Prop_handle.cpair(T_MA_inlet);
    inlet_air_props{6} = Ra * T_MA_inlet / (p_MA_inlet * 1e6);
end

sat_cache = sat_global;  % Phase1 压力不变，入口=出口饱和物性

for i = 1:loopmax
    out = R_cal_10(x0_R,BD_R,GeoCondition,CV,Prop_handle,[],inlet_props,sat_cache);

    HTC_R = out{3};
    dEF_R = out{5};
    T_R   = out{6};
    dp_CV = out{4};   % Pa，当前 CV 压降

    out_MA = DryA_cal_10(x0_MA,BD_MA,GeoCondition,CV,N,Prop_handle,inlet_air_props);

    dEF_MA = out_MA{1};
    n_fin  = out_MA{2};
    HTC_MA = out_MA{3};
    cp_MA  = out_MA{4};

    dT = T_R -(T_MA_inlet + x0(2))/2;
    UA_R = 1./(1./(HTC_R*A_R_CV)+1./(n_fin.*HTC_MA*A_MA_CV));
    Q = dT.*UA_R;

    if max(abs([(dEF_R*1e3-Q)/Q,(dEF_MA+Q)/Q]))<residual_Energy && i<loopmax
        x0_R(1) = h_R_inlet - Q/mdot_R_inlet/1e3;
        x0_MA(1) = Q/cp_MA/(mdot_MA_inlet) + T_MA_inlet;
        break
    elseif i==loopmax
        disp("换热不动点迭代失败")
    else
        x0_R(1) = h_R_inlet - Q/mdot_R_inlet/1e3;
        x0_MA(1) = Q/cp_MA/(mdot_MA_inlet) + T_MA_inlet;
    end
end

cache_R_out = out{7};
cache_MA_out = out_MA{7};
F = [x0_R(1);x0_MA(1)];
end



function [h_R_in, h_R_out, T_MA_in, T_MA_out, ...
          p_R_in, p_R_out, p_MA_in, p_MA_out, ...
          R_flow, dp_tube, residual_max] = scanTubes_phase1(...
    tubePaths, predecessors_in, ...
    h_R_in, h_R_out, T_MA_in, T_MA_out, ...
    p_R_in, p_R_out, p_MA_in, p_MA_out, ...
    mdot_R, mdot_MA, TCinf, GeoCondition, ...
    CV_num, row, Prop_handle, ...
    h_R_inlet, p_R_inlet, T_MA_inlet, p_MA_inlet, ...
    residual_max, R_coef, sat_global)

nTubes = length(tubePaths);
R_flow  = zeros(nTubes, 1);
dp_tube = zeros(nTubes, 1);   % 每管总压降 (Pa)

for i2 = 1:nTubes
    tube = tubePaths(i2);

    p_R_in_tube  = p_R_in(:, tube);
    p_R_out_tube = p_R_out(:, tube);
    h_R_in_tube  = h_R_in(:, tube);
    h_R_out_tube = h_R_out(:, tube);

    tube_upper = predecessors_in{tube};
    if size(tube_upper, 2) == 0
        h_R_in_tube(1) = h_R_inlet;
        p_R_in_tube(1) = p_R_inlet;
    elseif size(tube_upper, 2) == 1
        h_R_in_tube(1) = h_R_out(CV_num, tube_upper);
        p_R_in_tube(1) = p_R_out(CV_num, tube_upper);
    else
        p_R_in_tube(1) = p_R_out(CV_num, tube_upper(1));
        h_R_in_tube(1) = h_R_out(CV_num, tube_upper) ...
            * mdot_R(tube_upper) / sum(mdot_R(tube_upper));
    end

    mdot_R_tube  = mdot_R(tube);
    mdot_MA_tube = mdot_MA(:, tube);

    T_MA_in_tube  = T_MA_in(:, tube);
    T_MA_out_tube = T_MA_out(:, tube);
    p_MA_in_tube  = p_MA_in(:, tube);
    p_MA_out_tube = p_MA_out(:, tube);

    flowDirection = TCinf.FlowDirection(tube);

    dp_acc = 0;

    cache_R_in = {}; cache_MA_in = {};
    for i3 = 1:CV_num
        rev = CV_num + 1 - i3;

        if flowDirection == 1
            x0  = [h_R_out_tube(i3); T_MA_out_tube(i3)];
            BD  = [p_R_out_tube(i3); p_MA_out_tube(i3)];
            InBD = [h_R_in_tube(i3); p_R_in_tube(i3); ...
                    T_MA_in_tube(i3); p_MA_in_tube(i3); ...
                    mdot_R_tube; mdot_MA_tube(i3)];
        else
            x0  = [h_R_out_tube(i3); T_MA_out_tube(rev)];
            BD  = [p_R_out_tube(i3); p_MA_out_tube(rev)];
            InBD = [h_R_in_tube(i3); p_R_in_tube(i3); ...
                    T_MA_in_tube(rev); p_MA_in_tube(rev); ...
                    mdot_R_tube; mdot_MA_tube(rev)];
        end

        N1 = ceil(tube / row);
        [xout, dp_CV, cache_R_out, cache_MA_out] = alg_phase1(x0, BD, InBD, GeoCondition, CV_num, N1, Prop_handle, cache_R_in, cache_MA_in, sat_global);
        residual_max = max(max(abs(xout(1:2) - x0) ./ x0), residual_max);

        h_R_out_tube(i3) = xout(1);
        dp_acc = dp_acc + dp_CV;

        if flowDirection == 1
            T_MA_out(i3, tube) = xout(2);
        else
            T_MA_out(rev, tube) = xout(2);
        end
        if i3 < CV_num
            h_R_in_tube(i3 + 1) = xout(1);
            cache_R_in = cache_R_out;  % 工质沿管流动，出口→下游入口
            % 空气横掠管束，同管所有 CV 入口相同，不转发 cache_MA
        end
    end

    h_R_out(:, tube) = h_R_out_tube;
    h_R_in(:, tube)  = h_R_in_tube;
    T_MA_in = [T_MA_inlet * ones(CV_num, row), T_MA_out(:, 1:end-row)];
    p_MA_in = [p_MA_inlet * ones(CV_num, row), p_MA_out(:, 1:end-row)];

    dp_tube(tube) = dp_acc / 1e6;                      % Pa → MPa
    R_flow(tube)  = dp_acc ./ (mdot_R_tube .^ R_coef);  % Pa / (kg/s)^e
end
end



function [h_R_in, h_R_out, T_MA_in, T_MA_out, ...
          p_R_in, p_R_out, p_MA_in, p_MA_out, ...
          dp_tube, residual_max] = scanTubes(...
    tubePaths, predecessors_in, ...
    h_R_in, h_R_out, T_MA_in, T_MA_out, ...
    p_R_in, p_R_out, p_MA_in, p_MA_out, ...
    mdot_R, mdot_MA, TCinf, GeoCondition, ...
    CV_num, row, solver_flag, Prop_handle, ...
    h_R_inlet, p_R_inlet, T_MA_inlet, p_MA_inlet, ...
    residual_max, dp_tube)

nTubes = length(tubePaths);

for i2 = 1:nTubes
    tube = tubePaths(i2);

    p_R_in_tube  = p_R_in(:, tube);
    p_R_out_tube = p_R_out(:, tube);
    h_R_in_tube  = h_R_in(:, tube);
    h_R_out_tube = h_R_out(:, tube);

    tube_upper = predecessors_in{tube};
    if size(tube_upper, 2) == 0
        h_R_in_tube(1) = h_R_inlet;
        p_R_in_tube(1) = p_R_inlet;
    elseif size(tube_upper, 2) == 1
        h_R_in_tube(1) = h_R_out(CV_num, tube_upper);
        p_R_in_tube(1) = p_R_out(CV_num, tube_upper);
    else
        p_R_in_tube(1) = p_R_out(CV_num, tube_upper(1));
        h_R_in_tube(1) = h_R_out(CV_num, tube_upper) * mdot_R(tube_upper) / sum(mdot_R(tube_upper));
    end

    mdot_R_tube = mdot_R(tube);
    mdot_MA_tube = mdot_MA(:, tube);

    T_MA_in_tube = T_MA_in(:, tube);
    T_MA_out_tube = T_MA_out(:, tube);
    p_MA_in_tube = p_MA_in(:, tube);
    p_MA_out_tube = p_MA_out(:, tube);

    flowDirection = TCinf.FlowDirection(tube);

    for i3 = 1:CV_num
        rev = CV_num + 1 - i3;

        if solver_flag == 1
            if flowDirection == 1
                x0  = [h_R_out_tube(i3); T_MA_out_tube(i3)];
                BD  = [p_R_out_tube(i3); p_MA_out_tube(i3)];
                InBD = [h_R_in_tube(i3); p_R_in_tube(i3); T_MA_in_tube(i3); p_MA_in_tube(i3); ...
                        mdot_R_tube; mdot_MA_tube(i3)];
            else
                x0  = [h_R_out_tube(i3); T_MA_out_tube(rev)];
                BD  = [p_R_out_tube(i3); p_MA_out_tube(rev)];
                InBD = [h_R_in_tube(i3); p_R_in_tube(i3); T_MA_in_tube(rev); p_MA_in_tube(rev); ...
                        mdot_R_tube; mdot_MA_tube(rev)];
            end
        else
            if flowDirection == 1
                BD  = [h_R_out_tube(i3); T_MA_out_tube(i3)];
                x0  = [p_R_out_tube(i3); p_MA_out_tube(i3)];
                InBD = [h_R_in_tube(i3); p_R_in_tube(i3); T_MA_in_tube(i3); p_MA_in_tube(i3); ...
                        mdot_R_tube; mdot_MA_tube(i3)];
            else
                BD  = [h_R_out_tube(i3); T_MA_out_tube(rev)];
                x0  = [p_R_out_tube(i3); p_MA_out_tube(rev)];
                InBD = [h_R_in_tube(i3); p_R_in_tube(i3); T_MA_in_tube(rev); p_MA_in_tube(rev); ...
                        mdot_R_tube; mdot_MA_tube(rev)];
            end
        end

        N1 = ceil(tube / row);
        [xout, ~, ~] = alg(x0, BD, InBD, GeoCondition, CV_num, N1, solver_flag, Prop_handle, {}, {});
        residual_max = max(max(abs(xout - x0) ./ x0), residual_max);

        if solver_flag == 1
            h_R_out_tube(i3) = xout(1);
            if flowDirection == 1
                T_MA_out(i3, tube) = xout(2);
            else
                T_MA_out(rev, tube) = xout(2);
            end
            if i3 < CV_num
                h_R_in_tube(i3 + 1) = xout(1);
            end
        else
            p_R_out_tube(i3) = xout(1);
            p_R_out(i3, tube) = xout(1);
            if flowDirection == 1
                p_MA_out(i3, tube) = xout(2);
            else
                p_MA_out(rev, tube) = xout(2);
            end
            if i3 < CV_num
                p_R_in_tube(i3 + 1) = xout(1);
            end
        end
    end

    if solver_flag == 1
        h_R_out(:, tube) = h_R_out_tube;
        h_R_in(:, tube) = h_R_in_tube;
        T_MA_in = [T_MA_inlet * ones(CV_num, row), T_MA_out(:, 1:end-row)];
    else
        dp_tube(tube) = p_R_in_tube(1) - p_R_out_tube(CV_num);
        p_R_in(:, tube) = p_R_in_tube;
        p_MA_in = [p_MA_inlet * ones(CV_num, row), p_MA_out(:, 1:end-row)];
    end
end
end
