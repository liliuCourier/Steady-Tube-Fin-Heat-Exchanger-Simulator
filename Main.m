% 主求解器 — 运行前请先执行 PreProcessing 完成流路/几何/边界设置

% 添加子模块路径
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'Lib'));
addpath(fullfile(root, 'Solver'));
addpath(fullfile(root, 'PostProcessing'));
addpath(fullfile(root, 'PreProc'));

%%
% 检查预处理数据是否就绪
if ~exist('TCinf','var')
    error('TCinf 未找到，请先运行 PreProcessing');
end
if ~exist('GeoCondition','var')
    error('GeoCondition 未找到，请先运行 PreProcessing');
end
if ~exist('BDCondition','var')
    error('BDCondition 未找到，请先运行 PreProcessing');
end

%% 初始化求解

[N,u0,mdot0,mdot_R_init,exitflag] = mdot_Initial(BDCondition,GeoCondition,TCinf);
[heatPaths,pdropPaths,predecessors_in,predecessors_out] = buildPath(TCinf.TC_matrix,N);
%%
CV_num = BDCondition.CV_num;
row = GeoCondition.row;
Tube_num = GeoCondition.Tube_num;

% 边界条件
h_R_inlet       = BDCondition.BD_R.h_R_inlet;
p_R_inlet       = BDCondition.BD_R.p_R_inlet;
mdot_MA_inlet   = BDCondition.BD_MA.mdot_MA_inlet;
T_MA_inlet      = BDCondition.BD_MA.T_MA_inlet;
p_MA_inlet      = BDCondition.BD_MA.p_MA_inlet;

% 物性场初始化
p_R_in  = p_R_inlet*ones(CV_num,Tube_num);
p_R_out = p_R_inlet*ones(CV_num,Tube_num);
h_R_in  = h_R_inlet*ones(CV_num,Tube_num);
h_R_out = h_R_inlet*ones(CV_num,Tube_num);
mdot_R  = mdot_R_init;

% 空气侧
T_MA_in     = T_MA_inlet*ones(CV_num,Tube_num);
T_MA_out    = T_MA_inlet*ones(CV_num,Tube_num);
p_MA_in     = p_MA_inlet*ones(CV_num,Tube_num);
p_MA_out    = p_MA_inlet*ones(CV_num,Tube_num);
mdot_MA     = mdot_MA_inlet/row/CV_num*ones(CV_num,Tube_num);

% 管进出口物性初始化
h_R_tube_inlet  = h_R_inlet * ones(1, Tube_num);
h_R_tube_outlet = h_R_inlet * ones(1, Tube_num);
p_R_tube_inlet  = p_R_inlet * ones(1, Tube_num);
p_R_tube_outlet = p_R_inlet * ones(1, Tube_num);

%%
dp_tube = zeros(Tube_num,1);
R_flow = [];            % 初始化（无环路时保持为空）
K_bend  = 1.5;          % U型弯单相阻力系数 (180°回弯, R/D≈1.5)

% 求解设置
% 迭代上限设置
loopmax = 10;
residual_limit = 1e-3;

% 换热路径长度
length_heat = size(heatPaths,2);
%length_dp   = size(pdropPaths,2);

% 收敛历史记录
residual_history = zeros(loopmax, 1);
dp_loop_history   = zeros(loopmax, 1);
tube_cal = 0;
R_coef = 1.81;

% 每次 scanTubes 后的场量快照（用于收敛轨迹分析）
max_scans = loopmax * 3 + 1;
snapshot_h_out = cell(max_scans, 1);
snapshot_p_out = cell(max_scans, 1);
snapshot_mdot  = cell(max_scans, 1);
snapshot_dp    = cell(max_scans, 1);
snapshot_flag  = zeros(max_scans, 1);  % 1=换热, 2=压力
snapshot_iter  = zeros(max_scans, 1);  % 所属外层迭代编号

% 每次外层迭代后的场量快照（用于迭代级收敛轨迹）
iter_h_out = cell(loopmax, 1);
iter_p_out = cell(loopmax, 1);
iter_mdot  = cell(loopmax, 1);
iter_dp    = cell(loopmax, 1);

% 流量分配演进记录（u0 为环路流量解）
u0_history = cell(loopmax+1, 1);
u0_history{1} = u0;

%%

tic
for i1 = 1:loopmax
    % 记录迭代级场量快照（本次迭代的出发点）
    iter_h_out{i1} = h_R_out;
    iter_p_out{i1} = p_R_out;
    iter_mdot{i1}  = mdot_R;
    iter_dp{i1}    = dp_tube;

    residual_max = 0;
    % 换热路径扫描（广度优先）
    [h_R_in, h_R_out, T_MA_in, T_MA_out, ...
        p_R_in, p_R_out, p_MA_in, p_MA_out, ...
        dp_tube, residual_max, ...
        h_R_tube_outlet,h_R_tube_inlet, p_R_tube_outlet,p_R_tube_inlet] = scanTubes(...
        heatPaths, predecessors_in, predecessors_out, ...
        h_R_in, h_R_out, T_MA_in, T_MA_out, ...
        p_R_in, p_R_out, p_MA_in, p_MA_out, ...
        mdot_R, mdot_MA, TCinf, GeoCondition, ...
        CV_num, row, 1, Prop_handle, ...
        h_R_inlet, p_R_inlet, T_MA_inlet, p_MA_inlet, ...
        h_R_tube_inlet, p_R_tube_inlet, ...
        h_R_tube_outlet, p_R_tube_outlet, ...
        residual_max, dp_tube);
    tube_cal = tube_cal + 1;
    snapshot_h_out{tube_cal} = h_R_out;
    snapshot_p_out{tube_cal} = p_R_out;
    snapshot_mdot{tube_cal}  = mdot_R;
    snapshot_dp{tube_cal}    = dp_tube;
    snapshot_flag(tube_cal)  = 1;
    snapshot_iter(tube_cal)  = i1;


    % 检查是否有环路，如果没有环路，流量固定，不需要调整流量
    if isempty(N)
        % 无环路：直接过渡到第二次广度优先压力场更新
        [h_R_in, h_R_out, T_MA_in, T_MA_out, ...
            p_R_in, p_R_out, p_MA_in, p_MA_out, ...
            dp_tube, residual_max, ...
            h_R_tube_outlet,h_R_tube_inlet, p_R_tube_outlet,p_R_tube_inlet] = scanTubes(...
            heatPaths, predecessors_in, predecessors_out, ...
            h_R_in, h_R_out, T_MA_in, T_MA_out, ...
            p_R_in, p_R_out, p_MA_in, p_MA_out, ...
            mdot_R, mdot_MA, TCinf, GeoCondition, ...
            CV_num, row, 2, Prop_handle, ...
            h_R_inlet, p_R_inlet, T_MA_inlet, p_MA_inlet, ...
            h_R_tube_inlet, p_R_tube_inlet, ...
            h_R_tube_outlet, p_R_tube_outlet, ...
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

        % 有环路时需要更新流量
    else

        % 再更新一次压力场
        %for i2 = 1:loopmax
        [h_R_in, h_R_out, T_MA_in, T_MA_out, ...
            p_R_in, p_R_out, p_MA_in, p_MA_out, ...
            dp_tube, residual_max, ...
            h_R_tube_outlet,h_R_tube_inlet, p_R_tube_outlet,p_R_tube_inlet] = scanTubes(...
            heatPaths, predecessors_in, predecessors_out, ...
            h_R_in, h_R_out, T_MA_in, T_MA_out, ...
            p_R_in, p_R_out, p_MA_in, p_MA_out, ...
            mdot_R, mdot_MA, TCinf, GeoCondition, ...
            CV_num, row, 2, Prop_handle, ...
            h_R_inlet, p_R_inlet, T_MA_inlet, p_MA_inlet, ...
            h_R_tube_inlet, p_R_tube_inlet, ...
            h_R_tube_outlet, p_R_tube_outlet, ...
            residual_max, dp_tube);
        tube_cal = tube_cal + 1;
        snapshot_h_out{tube_cal} = h_R_out;
        snapshot_p_out{tube_cal} = p_R_out;
        snapshot_mdot{tube_cal}  = mdot_R;
        snapshot_dp{tube_cal}    = dp_tube;
        snapshot_flag(tube_cal)  = 2;
        snapshot_iter(tube_cal)  = i1;

        % 第二次更新流量场 — Newton 迭代求解 N'·(R·m^e) = 0
        % R_flow 由当前 dp_tube 和 mdot_R 确定，求解过程中保持不变
        dp_Pa = dp_tube * 1e6;
        R_flow = dp_Pa ./ (mdot_R.^R_coef);
        u = u0;
        for k = 1:10
            m_k   = mdot0 + N * u;
            F     = N' * (R_flow .* m_k.^R_coef);          % 环路残差
            S     = R_coef * R_flow .* m_k.^(R_coef - 1);  % 灵敏度 d(dp)/dm
            J     = N' * (S .* N);                          % Jacobian
            du    = -J \ F;
            u     = u + du;
            if norm(du) < 1e-8, break; end
        end
        mdot_R = mdot0 + N * u;
        u0 = u;
        u0_history{i1+1} = u0;

        if any(mdot_R <0)
            disp("流量算出来个负数，请检查压降的计算是不是存在问题，是不是有的管初始流量给少了，算出来的管压损为负数")
        end
        % 环路压降收敛判断
        %     if max(abs((dp_tube')*N)) < 1e-6
        %         disp("压降收敛")
        %         break
        %     end
        % end

        residual_history(i1) = residual_max;
        dp_loop_history(i1) = max(abs((dp_tube')*N));
        % 收敛准则为环路压降最大差异小于1Pa &&
        % 前后两次计算的值的相对误差<residual_limit,所有值统一阈值限，包含空气侧的温度、压力、工质侧的焓、压力
        % 内部如果报了压降不收敛和能量守恒计算不收敛时，认为后续计算均错误
        % 内部能量守恒计算和压降计算的守恒性是需要保证的，这个直接根据不动点迭代是否收敛进行判断，在此基础上才是外层迭代的准确性
        if max(abs((dp_tube')*N)) < 1e-6 && (residual_max < residual_limit)
            disp("压降收敛")
            break
        end

    end
end
time = toc;

% 导出最终收敛指标
dp_loop_max = dp_loop_history(i1);

% 每根管的最终合理的压阻系数：
% log(dp_tube*1e6./R_flow)./log(mdot_R)
fprintf('求解完成，耗时 %.2f s，正在运行后处理...\n', time);

% % 自动运行后处理
try
    PostProcessing;
catch ME
    fprintf('后处理运行出错: %s\n', ME.message);
end


%%
function  F = alg(x0,BD,InletBD,GeoCondition,CV,N,solver_flag,Prop_handle)
% 迭代上限
loopmax = 100;
% 残差阈值上限：
residual_Energy = 1e-3;
residual_dp = 1e-3;

% if solver_flag == 3
%     h_u = InletBD(1);  p_in_u = InletBD(2);  mdot_u = InletBD(3);
%     D = GeoCondition.D_inner;  S = pi*D^2/4;  G = abs(mdot_u)/S;
%     p_out_u = x0(1);
%     for i = 1:15
%         p_avg = (p_in_u + p_out_u)/2;
%         [~,~,~,~,x_u,~,v_L,v_V,~,~,~,~,~,~] = Prop1(real(p_avg), h_u, Prop_handle);
%         rho_L = 1/v_L;  dp_u = K_bend * G^2 / (2*rho_L);
%         if x_u > 0.05 && x_u < 0.95
%             rho_G = 1/v_V;
%             dp_u = (1 + x_u*(rho_L-rho_G)/rho_G) * dp_u;
%         end
%         p_new = p_in_u - dp_u/1e6;
%         if abs(p_new-p_out_u) < 1e-8, p_out_u = p_new; break; end
%         p_out_u = p_new;
%     end
%     F = [0; p_out_u];  return
% end

h_R_inlet = InletBD(1);
p_R_inlet = InletBD(2);

T_MA_inlet = InletBD(3);
p_MA_inlet = InletBD(4);

mdot_R_inlet  = InletBD(5);
mdot_MA_inlet = InletBD(6);

BD_R.h_R_inlet      = h_R_inlet;
BD_R.mdot_R_inlet   = mdot_R_inlet;
BD_R.p_R_inlet      = p_R_inlet;

BD_MA.T_MA_inlet    = T_MA_inlet;
BD_MA.mdot_MA_inlet = mdot_MA_inlet;
BD_MA.p_MA_inlet    = p_MA_inlet;


A_R = GeoCondition.A_R;
A_MA = GeoCondition.A_MA;
Tube_num = GeoCondition.Tube_num;

A_R_CV = A_R/CV/Tube_num;
A_MA_CV = A_MA/CV/Tube_num;


if solver_flag == 1

    x0_R =  [x0(1);BD(1)];
    x0_MA = [x0(2);BD(2)];

    %     % inlet property cache: inlet p/h unchanged during fixed-point iteration
    % 进口物性缓存：饱和物性提前从 Prop_handle 计算，同时用于入口
    % Prop1（跳过 11 次插值）和出口 sat_cache（p_in≈p_out, dp~Pa）
    sat_in = cell(1,11);
    sat_in{1} = Prop_handle.v_liq(0, p_R_inlet);    % vsatliq
    sat_in{2} = Prop_handle.v_vap(1, p_R_inlet);    % vsatvap
    sat_in{3} = Prop_handle.Pr_liq(0, p_R_inlet);   % Prsatliq
    sat_in{4} = Prop_handle.Pr_vap(1, p_R_inlet);   % Prsatvap
    sat_in{5} = Prop_handle.Nu_liq(0, p_R_inlet);   % Nusatliq
    sat_in{6} = Prop_handle.Nu_vap(1, p_R_inlet);   % Nusatvap
    sat_in{7} = Prop_handle.k_liq(0, p_R_inlet);    % ksatliq
    sat_in{8} = Prop_handle.k_vap(1, p_R_inlet);    % ksatvap
    sat_in{9} = Prop_handle.T_liq(0, p_R_inlet);    % Tsatliq
    sat_in{10} = Prop_handle.h_sat_liq(p_R_inlet);  % hsatliq
    sat_in{11} = Prop_handle.h_sat_vap(p_R_inlet);  % hsatvap

    inlet_props = cell(1,14);
    [inlet_props{:}] = Prop1(p_R_inlet, h_R_inlet, Prop_handle, sat_in);
    Ra = 287.047;
    inlet_air_props = cell(1,6);
    inlet_air_props{1} = Prop_handle.hair(T_MA_inlet);
    inlet_air_props{2} = Prop_handle.Prair(T_MA_inlet);
    inlet_air_props{3} = Prop_handle.visair(T_MA_inlet);
    inlet_air_props{4} = Prop_handle.kair(T_MA_inlet);
    inlet_air_props{5} = Prop_handle.cpair(T_MA_inlet);
    inlet_air_props{6} = Ra * T_MA_inlet / (p_MA_inlet * 1e6);

    % 出口饱和缓存：复用 sat_in，仅 Tsatliq 重取 p_out
    sat_cache = cell(1,11);
    sat_cache{1} = sat_in{1}; sat_cache{2} = sat_in{2};
    sat_cache{3} = sat_in{3}; sat_cache{4} = sat_in{4};
    sat_cache{5} = sat_in{5}; sat_cache{6} = sat_in{6};
    sat_cache{7} = sat_in{7}; sat_cache{8} = sat_in{8};
    sat_cache{9} = Prop_handle.T_liq(0, BD(1));  % Tsatliq at p_out
    sat_cache{10} = sat_in{10};  % hsatliq: p_in≈p_out
    sat_cache{11} = sat_in{11};  % hsatvap: p_in≈p_out

    for i = 1:loopmax
        out = R_cal_10(x0_R,BD_R,GeoCondition,CV,Prop_handle,inlet_props,sat_cache);

        % 获得计算结果
        %Tin_R        = out{1};
        %Tout_R       = out{2};
        HTC_R        = out{3};
        dEF_R        = out{5};
        T_R          = out{6};

        out_MA = DryA_cal_10(x0_MA,BD_MA,GeoCondition,CV,N,Prop_handle,inlet_air_props);

        dEF_MA =    out_MA{1};
        n_fin =     out_MA{2};
        HTC_MA =    out_MA{3};
        cp_MA =     out_MA{4};

        dT = T_R -(T_MA_inlet + x0(2))/2;
        UA_R = 1./(1./(HTC_R*A_R_CV)+1./(n_fin.*HTC_MA*A_MA_CV));

        Q = dT.*UA_R;
        if max(abs([(dEF_R*1e3 - Q)/Q,(dEF_MA + Q)/Q]))<residual_Energy && i<loopmax
            x0_R(1) = h_R_inlet - Q/mdot_R_inlet/1e3;
            x0_MA(1) = Q/cp_MA/(mdot_MA_inlet) + T_MA_inlet;
            break
        elseif i==loopmax
            disp("换热不动点迭代失败")
        else
            % 根据换热反算焓值
            x0_R(1) = h_R_inlet - Q/mdot_R_inlet/1e3;
            x0_MA(1) = Q/cp_MA/(mdot_MA_inlet) + T_MA_inlet;
        end
    end

    F = [x0_R(1);x0_MA(1)];

elseif solver_flag == 2


    x0_R  = [BD(1);x0(1)];
    x0_MA = [BD(2);x0(2)];

    % 压力扫描中做对应的进口物性缓存，但是饱和物性不能缓存了
    % inlet property cache: inlet p/h unchanged during fixed-point iteration
    inlet_props = cell(1,14);
    [inlet_props{:}] = Prop1(p_R_inlet, h_R_inlet, Prop_handle);
    Ra = 287.047;
    inlet_air_props = cell(1,6);
    inlet_air_props{1} = Prop_handle.hair(T_MA_inlet);
    inlet_air_props{2} = Prop_handle.Prair(T_MA_inlet);
    inlet_air_props{3} = Prop_handle.visair(T_MA_inlet);
    inlet_air_props{4} = Prop_handle.kair(T_MA_inlet);
    inlet_air_props{5} = Prop_handle.cpair(T_MA_inlet);
    inlet_air_props{6} = Ra * T_MA_inlet / (p_MA_inlet * 1e6);

    for i = 1:loopmax

        out = R_cal_10(x0_R,BD_R,GeoCondition,CV,Prop_handle,inlet_props);
        dp_R        = out{4};

        out_MA = DryA_cal_10(x0_MA,BD_MA,GeoCondition,CV,N,Prop_handle,inlet_air_props);

        dp_MA =    out_MA{5};
        if max(abs([(dp_R - (p_R_inlet - x0_R(2))*1e6)/dp_R,(dp_MA - (p_MA_inlet - x0_MA(2))*1e6)/dp_MA]))<residual_dp && i<loopmax
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
    F = [x0_R(2);x0_MA(2)];
end

end

%%

function [h_R_in, h_R_out, T_MA_in, T_MA_out, ...
    p_R_in, p_R_out, p_MA_in, p_MA_out, ...
    dp_tube, residual_max, ...
    h_R_tube_outlet,h_R_tube_inlet, p_R_tube_outlet,p_R_tube_inlet] = scanTubes(...
    tubePaths, predecessors_in, predecessors_out, ...
    h_R_in, h_R_out, T_MA_in, T_MA_out, ...
    p_R_in, p_R_out, p_MA_in, p_MA_out, ...
    mdot_R, mdot_MA, TCinf, GeoCondition, ...
    CV_num, row, solver_flag, Prop_handle, ...
    h_R_inlet, p_R_inlet, T_MA_inlet, p_MA_inlet, ...
    h_R_tube_inlet, p_R_tube_inlet, ...
    h_R_tube_outlet, p_R_tube_outlet, ...
    residual_max, dp_tube)

nTubes = length(tubePaths);

% 逐管计算
for i2 = 1:nTubes
    tube = tubePaths(i2);

    % 传入控制体参数，作为迭代初值使用
    mdot_R_tube = mdot_R(tube);
    mdot_MA_tube = mdot_MA(:, tube);

    p_R_in_tube  = p_R_in(:, tube);
    p_R_out_tube = p_R_out(:, tube);
    h_R_in_tube  = h_R_in(:, tube);
    h_R_out_tube = h_R_out(:, tube);

    p_MA_in_tube =  p_MA_in(:, tube);
    p_MA_out_tube = p_MA_out(:, tube);
    T_MA_in_tube =  T_MA_in(:, tube);
    T_MA_out_tube = T_MA_out(:, tube);


    % 根据上游信息调整管的进口条件
    tube_upper = predecessors_in{tube};
    if size(tube_upper, 2) == 0
        h_R_tube_inlet(tube) = h_R_inlet;
        p_R_tube_inlet(tube) = p_R_inlet;
    elseif size(tube_upper, 2) == 1
        h_R_tube_inlet(tube) = h_R_tube_outlet(tube_upper);
        p_R_tube_inlet(tube) = p_R_tube_outlet(tube_upper);
    else
        p_R_tube_inlet(tube) = p_R_tube_outlet(tube_upper(1));
        h_R_tube_inlet(tube) = h_R_tube_outlet(tube_upper) * mdot_R(tube_upper) / sum(mdot_R(tube_upper));
    end

    % 管的第一个控制体暂且先使用管进口参数
    h_R_in_tube(1) = h_R_tube_inlet(tube);
    p_R_in_tube(1) = p_R_tube_inlet(tube);

    flowDirection = TCinf.FlowDirection(tube);

    % === bend_in: 上游分流 → CV1入口含弯管压降 ===
    % b AI就喜欢写子函数来进行调用
    is_bi = ~isempty(predecessors_in{tube}) && ...
        length(predecessors_out{predecessors_in{tube}(1)}) > 1;
    if is_bi  && solver_flag == 2
        Lb = bend_len(tube, predecessors_in{tube}(1), GeoCondition);
        dp_u = bend_cal(h_R_in_tube(1), p_R_in_tube(1), ...
            mdot_R_tube, Lb, GeoCondition, Prop_handle);            % MPa
        p_R_in_tube(1) = p_R_tube_inlet(tube) - dp_u;
    end

    % 逐控制体计算
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
        xout = alg(x0, BD, InBD, GeoCondition, CV_num, N1, solver_flag, Prop_handle);
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

    % 这里是一定要把参数传递进来，不然变量就会一直变
    h_R_tube_outlet(tube) = h_R_out_tube(end);
    p_R_tube_outlet(tube) = p_R_out_tube(end);

    % === bend_out: 单一下游 → 管出口含弯管压降 ===
    is_bo = (length(predecessors_out{tube}) == 1);
    if is_bo && solver_flag == 2
        Lb = bend_len(tube, predecessors_out{tube}(1), GeoCondition);
        dp_u = bend_cal(h_R_tube_outlet(tube), ...
            p_R_out_tube(end), mdot_R_tube, Lb, GeoCondition, Prop_handle);
        p_R_tube_outlet(tube) = p_R_out_tube(end) - dp_u;
    end

    if solver_flag == 1
        h_R_out(:, tube) = h_R_out_tube;
        h_R_in(:, tube) = h_R_in_tube;
        T_MA_in = [T_MA_inlet * ones(CV_num, row), T_MA_out(:, 1:end-row)];
    else
        % 这里AI把dp_tube弄错了，显然是管的进口减出口
        dp_tube(tube) = p_R_tube_inlet(tube) - p_R_tube_outlet(tube);
        p_R_in(:, tube) = p_R_in_tube;
        p_MA_in = [p_MA_inlet * ones(CV_num, row), p_MA_out(:, 1:end-row)];
    end
    % if any(dp_tube<0)
    %     pause()
    % end
end
end

function dp_out = bend_cal(h_in, p_in, mdot, L_geom, Geo, Ph)
% 弯管压降: L_total=几何长+当量长, Xu-Fang 2013 两相修正, 返回 MPa
D_inner = Geo.D_inner;  S = pi*D_inner^2/4;  G = abs(mdot)/S;
p_out = p_in;
for i = 1:100
    p_avg = (p_in + p_out)/2;
    [~,~,~,~, x_CV, ~,vsatliq_CV,vsatvap_CV,...
         ~, ~, Nuliq, Nuvap,~, ~] = Prop1(p_avg, h_in, Ph);

    Re_go = G*D_inner/(Nuvap*1e-6/vsatvap_CV);
    Re_lo = G*D_inner/(Nuliq*1e-6/vsatliq_CV);

    Re_lam = 2000;  Re_tur = 3000;

    r = 0;
    % f_go = 8*((8/Re_go)^12 + 1/((2.457*log((7/Re_go)^0.9 + 0.27*(r/D_inner)))^16 + (37530/Re_go)^16)^(3/2))^(1/12);
    % f_lo = 8*((8/Re_lo)^12 + 1/((2.457*log((7/Re_lo)^0.9 + 0.27*(r/D_inner)))^16 + (37530/Re_lo)^16)^(3/2))^(1/12);
    if Re_go <= Re_lam
        f_go = 64 / Re_go;
    elseif Re_go > Re_tur
        f_go = 0.25*(log10(150.39/Re_go^0.98865-152.66/Re_go))^(-2);
    else
        f_go = (1.1525*Re_go + 895)*1e-5;
    end

    if Re_lo <= Re_lam
        f_lo = 64 / Re_lo;
    elseif Re_lo > Re_tur
        f_lo = 0.25*(log10(150.39/Re_lo^0.98865-152.66/Re_lo))^(-2);
    else
        f_lo = (1.1525*Re_lo + 895)*1e-5;
    end

    dpdL_lo = f_lo*(abs(mdot)/S)^2/(2*D_inner/vsatliq_CV);
    dpdL_go = f_go*(abs(mdot)/S)^2/(2*D_inner/vsatvap_CV);

    % Xu-Fang 2013 冷凝两相摩擦压降(NED 263, 87-96)
    % 适用范围: R134a,R22,R410A 等, Dh 0.1–10mm, G 20–800, q 2–55.3
    S_slip = (vsatvap_CV/vsatliq_CV)^(1/3);
    alpha = x_CV*vsatvap_CV/(x_CV*vsatvap_CV + (1-x_CV)*vsatliq_CV);
    x_dyn = 1/(1 + (1-alpha)/alpha * (vsatvap_CV/vsatliq_CV)/S_slip);
   
    Y = sqrt(dpdL_go / dpdL_lo);
    rho_tp = 1 / (x_dyn*vsatvap_CV + (1-x_dyn)*vsatliq_CV);
    g_acc = 9.81;
    Fr_tp = (abs(mdot)/S)^2 / (g_acc * D_inner * rho_tp^2);
    sigma = 0.002;
    We_tp = (abs(mdot)/S)^2 * D_inner / (rho_tp * sigma);
    phi2_lo = Y^2 * x_dyn^3 + (1 - x_dyn^2.59)^0.632 * ...
        (1 + 2*x_dyn^1.17*(Y^2 - 1) + 0.00775*x_dyn^(-0.475)*Fr_tp^0.535*We_tp^0.188);

    L = 30*D_inner;
    dp_b = phi2_lo * dpdL_lo * L;

    % if dp_b <0
    %     pause()
    % end

    p_new = p_in - dp_b/1e6;
    if abs(p_new-p_out) < 1e-6, p_out = p_new; break;
    elseif i == 50, disp("U型弯压力迭代失败");
    end
    p_out = p_new;
end
dp_out = (p_in - p_out);
end

function Lb = bend_len(k1, k2, Gc)
% 管间U型弯几何长度: 弧长 + 直线段
col = Gc.col;  P_row = Gc.P_row;  P_col = Gc.P_col;
r1 = floor((k1-1)/col)+1;  c1 = mod(k1-1, col)+1;
r2 = floor((k2-1)/col)+1;  c2 = mod(k2-1, col)+1;
dist = sqrt(((r2-r1)*P_row)^2 + ((c2-c1)*P_col)^2);
Rb = 1.5 * Gc.D_inner;
straight = max(0, dist - 2*Rb);
arc = pi * Rb;
Lb = straight + arc;
end