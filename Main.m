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

% 生成广度优先路径和环路优先路径
[heatPaths,pdropPaths,predecessors_in,predecessors_out] = buildPath(TCinf.TC_matrix,N);
%%

CV_num = BDCondition.CV_num;
row = GeoCondition.row;
Tube_num = GeoCondition.Tube_num;
%Uband_length = GeoCondition.Uband_length;

% 初始化其余热力场和压力场
h_R_inlet       = BDCondition.BD_R.h_R_inlet;
p_R_inlet       = BDCondition.BD_R.p_R_inlet;

% 非均匀的风场请在这设置
mdot_MA_inlet   = BDCondition.BD_MA.mdot_MA_inlet;
T_MA_inlet      = BDCondition.BD_MA.T_MA_inlet;
p_MA_inlet      = BDCondition.BD_MA.p_MA_inlet;

% 仍然使用流向判断，空气侧的排列向内，工质则是沿着流向，同排列为1，非同为0，生成的初始条件矩阵
% 工质侧
p_R_in  = p_R_inlet*ones(CV_num,Tube_num);
p_R_out = p_R_inlet*ones(CV_num,Tube_num);

h_R_in  = h_R_inlet*ones(CV_num,Tube_num);
h_R_out = h_R_inlet*ones(CV_num,Tube_num);

mdot_R  = mdot_R_init;

% 管进出口物性初始化
h_R_tube_inlet  = h_R_inlet * ones(1, Tube_num);
h_R_tube_outlet = h_R_inlet * ones(1, Tube_num);
p_R_tube_inlet  = p_R_inlet * ones(1, Tube_num);
p_R_tube_outlet = p_R_inlet * ones(1, Tube_num);

% 空气侧
T_MA_in     = T_MA_inlet*ones(CV_num,Tube_num);
T_MA_out    = T_MA_inlet*ones(CV_num,Tube_num);

p_MA_in     = p_MA_inlet*ones(CV_num,Tube_num);
p_MA_out    = p_MA_inlet*ones(CV_num,Tube_num);

% 非均匀的风场请在这设置
mdot_MA     = mdot_MA_inlet/row/CV_num*ones(CV_num,Tube_num);

dp_tube = zeros(Tube_num,1);
R_flow = [];  % 初始化（无环路时保持为空）
%dp_Uband = zeros(Tube_num,1);

% 求解
% 迭代上限设置
loopmax = 10;
% 残差设置
residual_limit = 1e-3;

% 换热路径长度
length_heat = size(heatPaths,2);
length_dp   = size(pdropPaths,2);

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
         h_R_tube_outlet, p_R_tube_outlet] = scanTubes(...
        heatPaths, predecessors_in, ...
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
            dp_tube, residual_max] = scanTubes(...
            heatPaths, predecessors_in, ...
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
            dp_tube, residual_max] = scanTubes(...
            heatPaths, predecessors_in, ...
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

        % 原 fsolve 非线性求解（保留备查）
        % u = fsolve(@(u)uF(u,R_flow,N,mdot0,R_coef),u0,options);
        % mdot_R = mdot0 + N*u;

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
% 
% % 自动运行后处理
try
    PostProcessing;
catch ME
    fprintf('后处理运行出错: %s\n', ME.message);
end

%%
% function F = uF(u,R_flow,N,mdot0,R_coef)
% dp_tube = (R_flow).*(mdot0 + N*u).^R_coef;
% dp_loop = (dp_tube')*N;
% F = dp_loop;
% end

function  F = alg(x0,BD,InletBD,GeoCondition,CV,N,solver_flag,Prop_handle)
% 迭代上限
loopmax = 100;
% 残差阈值上限：
residual_Energy = 1e-3;
residual_dp = 1e-3;

% if solver_flag == 3
% 
%     for i = 1:loopmax
%         x0_R  = [BD(1);x0(1)];
%         BD_R.h_R_inlet      = BD(2);
%         BD_R.mdot_R_inlet   = InletBD(1);
%         BD_R.p_R_inlet      = x0(2);
%         BD_R.Uband_L =      InletBD(2);
% 
%         flag = 2;
%         out = R_cal_10(x0_R,BD_R,GeoCondition,CV,Prop_handle,flag);
%         dp_R        = out;
%         if max(abs((dp_R - (x0(2) - x0(1))*1e6)/dp_R))<residual_dp && i<loopmax
%             x0(1) = x0(2) - dp_R/1e6;
%             break
%         elseif i==loopmax
%             disp("Uband压降不动点迭代失败")
%         else
%             x0(1) = x0(2) - dp_R/1e6;
%         end
%     end
% 
%     F = [dp_R,x0(1)];
%     return
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
        out = R_cal_10(x0_R,BD_R,GeoCondition,CV,Prop_handle,[],inlet_props,sat_cache);

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

        out = R_cal_10(x0_R,BD_R,GeoCondition,CV,Prop_handle,[],inlet_props);
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


function [h_R_in, h_R_out, T_MA_in, T_MA_out, ...
          p_R_in, p_R_out, p_MA_in, p_MA_out, ...
          dp_tube, residual_max, ...
          h_R_tube_outlet, p_R_tube_outlet] = scanTubes(...
    tubePaths, predecessors_in, ...
    h_R_in, h_R_out, T_MA_in, T_MA_out, ...
    p_R_in, p_R_out, p_MA_in, p_MA_out, ...
    mdot_R, mdot_MA, TCinf, GeoCondition, ...
    CV_num, row, solver_flag, Prop_handle, ...
    h_R_inlet, p_R_inlet, T_MA_inlet, p_MA_inlet, ...
    h_R_tube_inlet, p_R_tube_inlet, ...
    h_R_tube_outlet, p_R_tube_outlet, ...
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
        h_R_tube_inlet(tube) = h_R_inlet;
        p_R_tube_inlet(tube) = p_R_inlet;
    elseif size(tube_upper, 2) == 1
        h_R_tube_inlet(tube) = h_R_tube_outlet(tube_upper);
        p_R_tube_inlet(tube) = p_R_tube_outlet(tube_upper);
    else
        p_R_tube_inlet(tube) = p_R_tube_outlet(tube_upper(1));
        h_R_tube_inlet(tube) = h_R_tube_outlet(tube_upper) * mdot_R(tube_upper) / sum(mdot_R(tube_upper));
    end
    h_R_in_tube(1) = h_R_tube_inlet(tube);
    p_R_in_tube(1) = p_R_tube_inlet(tube);

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

    h_R_tube_outlet(tube) = h_R_out_tube(end);
    p_R_tube_outlet(tube) = p_R_out_tube(end);

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