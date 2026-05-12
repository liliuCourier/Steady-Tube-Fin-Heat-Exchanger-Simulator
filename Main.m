% 主求解器 — 运行前请先执行 PreProcessing 完成流路/几何/边界设置
%%
% if ~exist('Prop_handle','var')
%     refprop_location = 'E:\refprop10\REFPROP';
%     R = 'R134a';
%     Prop_handle = Prop_load(refprop_location,R,1e-3,5.5,80,510,100,25,25);
% end

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

options = optimoptions('fsolve','Display','none',...
    'Algorithm','levenberg-marquardt',...
    'FunctionTolerance',1e-12,...
    'MaxFunctionEvaluations',5e4,...
    'StepTolerance',1e-8,...
    'UseParallel',false,...
    'ScaleProblem','jacobian');

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

% 空气侧
T_MA_in     = T_MA_inlet*ones(CV_num,Tube_num);
T_MA_out    = T_MA_inlet*ones(CV_num,Tube_num);

p_MA_in     = p_MA_inlet*ones(CV_num,Tube_num);
p_MA_out    = p_MA_inlet*ones(CV_num,Tube_num);

% 非均匀的风场请在这设置
mdot_MA     = mdot_MA_inlet/row/CV_num*ones(CV_num,Tube_num);

dp_tube = zeros(Tube_num,1);
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

tic
for i1 = 1:loopmax
    residual_max = 0;
    % 换热路径扫描（广度优先）
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
            residual_max, dp_tube);

        if residual_max < residual_limit
            disp("残差收敛")
            break
        end
    else
        % 有环路：环路压降扫描 + 流量重分配
        [h_R_in, h_R_out, T_MA_in, T_MA_out, ...
         p_R_in, p_R_out, p_MA_in, p_MA_out, ...
         dp_tube, residual_max] = scanTubes(...
            pdropPaths, predecessors_in, ...
            h_R_in, h_R_out, T_MA_in, T_MA_out, ...
            p_R_in, p_R_out, p_MA_in, p_MA_out, ...
            mdot_R, mdot_MA, TCinf, GeoCondition, ...
            CV_num, row, 2, Prop_handle, ...
            h_R_inlet, p_R_inlet, T_MA_inlet, p_MA_inlet, ...
            residual_max, dp_tube);

        % 根据环路压降重新计算流量分配
        R_flow = dp_tube*1e6./(mdot_R.^2);
        u = fsolve(@(u)uF(u,R_flow,N,mdot0),u0,options);
        mdot_R = mdot0 + N*u;
        u0 = u;

        % 再更新一次压力场
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

        % 第二次更新流量场
        R_flow = dp_tube*1e6./(mdot_R.^2);
        u = fsolve(@(u)uF(u,R_flow,N,mdot0),u0,options);
        mdot_R = mdot0 + N*u;

        % 环路压降收敛判断
        if max(abs((dp_tube')*N)) < 1e-6 && (residual_max < residual_limit)
            disp("压降收敛")
            break
        end
        u0 = u;
    end

    % 记录收敛历史
    residual_history(i1) = residual_max;
    if isempty(N)
        dp_loop_history(i1) = 0;
    else
        dp_loop_history(i1) = max(abs((dp_tube')*N));
    end
end
time = toc;

% 导出最终收敛指标
dp_loop_max = dp_loop_history(i1);

fprintf('求解完成，耗时 %.2f s，正在运行后处理...\n', time);

% 自动运行后处理
try
    PostProcessing;
catch ME
    fprintf('后处理运行出错: %s\n', ME.message);
end

%%
function F = uF(u,R_flow,N,mdot0)
 dp_tube = (R_flow).*(mdot0 + N*u).^2;
 dp_loop = (dp_tube')*N;
 F = dp_loop;
end


function  F = alg(x0,BD,InletBD,GeoCondition,CV,N,solver_flag,Prop_handle)

% 迭代上限
loopmax = 100;
% 残差阈值上限：
residual_Energy = 1e-3;
residual_dp = 1e-3;

if solver_flag == 3

    for i = 1:loopmax
        x0_R  = [BD(1);x0(1)];
        BD_R.h_R_inlet      = BD(2);
        BD_R.mdot_R_inlet   = InletBD(1);
        BD_R.p_R_inlet      = x0(2);
        BD_R.Uband_L =      InletBD(2);

        flag = 2;
        out = R_cal_10(x0_R,BD_R,GeoCondition,CV,Prop_handle,flag);
        dp_R        = out;
        if max(abs((dp_R - (x0(2) - x0(1))*1e6)/dp_R))<residual_dp && i<loopmax
            x0(1) = x0(2) - dp_R/1e6;
            break
        elseif i==loopmax
            disp("Uband压降不动点迭代失败")
        else
            x0(1) = x0(2) - dp_R/1e6;
        end
    end

    F = [dp_R,x0(1)];
    return
end

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
    

    for i = 1:loopmax
        out = R_cal_10(x0_R,BD_R,GeoCondition,CV,Prop_handle,[]);

        % 获得计算结果
        %Tin_R        = out{1};
        %Tout_R       = out{2};
        HTC_R        = out{3};
        dEF_R        = out{5};
        T_R          = out{6};

        out_MA = DryA_cal_10(x0_MA,BD_MA,GeoCondition,CV,N,Prop_handle);

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
        else% 根据换热反算焓值
            x0_R(1) = h_R_inlet - Q/mdot_R_inlet/1e3;
            x0_MA(1) = Q/cp_MA/(mdot_MA_inlet) + T_MA_inlet;
        end
    end

    F = [x0_R(1);x0_MA(1)];

elseif solver_flag == 2
    

    x0_R  = [BD(1);x0(1)];
    x0_MA = [BD(2);x0(2)];

    for i = 1:loopmax

        out = R_cal_10(x0_R,BD_R,GeoCondition,CV,Prop_handle,[]);
        dp_R        = out{4};

        out_MA = DryA_cal_10(x0_MA,BD_MA,GeoCondition,CV,N,Prop_handle);

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