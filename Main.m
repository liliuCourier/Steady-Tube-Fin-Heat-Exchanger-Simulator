% 主求解器
%%
if ~exist('Prop_handle','var')
    refprop_location = 'E:\refprop10\REFPROP';
    R = 'R134a';
    %示例：
    % p-Mpa\h-kJ/kg
    % Prop_handle = Prop_load(refprop_location,R,pmin,pmax,hmin,hmax,p_point,u_vap_point,u_liq_point);
    Prop_handle = Prop_load(refprop_location,R,1e-3,5.5,80,510,100,25,25);
end

options = optimoptions('fsolve','Display','none',...
    'Algorithm','levenberg-marquardt',...
    'FunctionTolerance',1e-12,...
    'MaxFunctionEvaluations',5e4,...
    'StepTolerance',1e-8,...
    'UseParallel',false,...
    'ScaleProblem','jacobian');

% % 唤起管路连接程序,建议不要关
% HX_Path_Planner1()

%% 在没有修改管路时，这些不需要重复加载
% 根据管路连接信息生成几何信息
[GeoCondition] = GenerateGeo(TCinf);
[BDCondition] = GenerateBD(GeoCondition,Prop_handle);

[N,u0,mdot0,mdot_R_init,exitflag] = mdot_Initial(BDCondition,GeoCondition,TCinf);

% 生成广度优先路径和环路优先路径
[heatPaths,pdropPaths,predecessors_in,predecessors_out] = buildPath(TCinf.TC_matrix,N);
%%

CV_num = 20;
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

tic
for i1 = 1:loopmax
    residual_max = 0;
    % 换热路径扫描
    for i2 = 1 : length_heat
        % 找到对应管的上游和下游，将边界条件、初始条件对应传入，随后逐控制体计算
        % 根据上游信息修改边界条件和初始条件
        tube = heatPaths(i2);

        p_R_in_tube  = p_R_in(:,tube);
        p_R_out_tube = p_R_out(:,tube);

        h_R_in_tube  = h_R_in(:,tube);
        h_R_out_tube = h_R_out(:,tube);
        
        % 工质侧更新上游混合情况
        tube_upper = predecessors_in{tube};
        if size(tube_upper,2) == 0      % 进口为总进口物性
           h_R_in_tube(1) = h_R_inlet;
           p_R_in_tube(1) = p_R_inlet;
        elseif size(tube_upper,2) == 1  % 内部管道直接传递
           h_R_in_tube(1) = h_R_out(CV_num,tube_upper);
           p_R_in_tube(1) = p_R_out(CV_num,tube_upper);
        else                            % 分管需要混合
           % 压力随便取一条管路计算
           p_R_in_tube(1) = p_R_out(CV_num,tube_upper(1));
           % 焓根据流量进行混合
           h_R_in_tube(1) = h_R_out(CV_num,tube_upper)*mdot_R(tube_upper)/sum(mdot_R(tube_upper));
        end

        mdot_R_tube = mdot_R(tube);
        mdot_MA_tube = mdot_MA(:,tube);

        % 空气侧暂时无混合
        T_MA_in_tube = T_MA_in(:,tube);
        T_MA_out_tube = T_MA_out(:,tube);

        p_MA_in_tube = p_MA_in(:,tube);
        p_MA_out_tube = p_MA_out(:,tube);

        % 为1则和空气控制体同向，为0则是反向
        flowDirection = TCinf.FlowDirection(tube);

        % 逐控制体计算
        % 严格来说，只要内部不报迭代未收敛的问题，就能可以以前后两次之间的差异作为迭代截至的判据
        for i3 = 1:CV_num
            if flowDirection == 1
                x0      = [h_R_out_tube(i3);T_MA_out_tube(i3)];
                P_BD    = [p_R_out_tube(i3);p_MA_out_tube(i3)];
                InletBD = [h_R_in_tube(i3);p_R_in_tube(i3);T_MA_in_tube(i3);p_MA_in_tube(i3);...
                            mdot_R_tube;mdot_MA_tube(i3)];
            else
                x0      = [h_R_out_tube(i3);T_MA_out_tube(CV_num + 1 - i3)];
                P_BD    = [p_R_out_tube(i3);p_MA_out_tube(CV_num + 1 - i3)];
                InletBD = [h_R_in_tube(i3);p_R_in_tube(i3);T_MA_in_tube(CV_num + 1 - i3);p_MA_in_tube(CV_num + 1 - i3);...
                            mdot_R_tube;mdot_MA_tube(CV_num + 1 - i3)];
            end
            
            flag = 1;
            N1 = ceil(tube/row) ;
            xout = alg(x0,P_BD,InletBD,GeoCondition,CV_num,N1,flag,Prop_handle);

            residual_max = max(max(abs(xout - x0)./x0),residual_max);

            h_R_out_tube(i3) = xout(1);
            if flowDirection == 1
               T_MA_out(i3,tube) = xout(2);
            else
               T_MA_out(CV_num + 1- i3,tube) = xout(2);
            end
            if i3 < CV_num
            h_R_in_tube(i3 + 1) = xout(1);
            end
        end
        h_R_out(:,tube) = h_R_out_tube;
        %T_MA_out(:,tube) = T_MA_out_tube;
        h_R_in(:,tube) = h_R_in_tube;
        T_MA_in = [T_MA_inlet*ones(CV_num,row),T_MA_out(:,1:end-row)];
        
    end
    


    % 环路压降算法
    for i2 = 1 : length_dp
        % 找到对应管的上游和下游，将边界条件、初始条件对应传入，随后逐控制体计算
        % 根据上游信息修改边界条件和初始条件
        tube = pdropPaths(i2);

        p_R_in_tube  = p_R_in(:,tube);
        p_R_out_tube = p_R_out(:,tube);

        h_R_in_tube  = h_R_in(:,tube);
        h_R_out_tube = h_R_out(:,tube);
        
        % 工质侧更新上游混合情况
        tube_upper = predecessors_in{tube};
        tube_lower = predecessors_out{tube};
        if size(tube_upper,2) == 0      % 进口为总进口物性
           h_R_in_tube(1) = h_R_inlet;
           p_R_in_tube(1) = p_R_inlet;
        elseif size(tube_upper,2) == 1  % 内部管道直接传递
           h_R_in_tube(1) = h_R_out(CV_num,tube_upper);
           p_R_in_tube(1) = p_R_out(CV_num,tube_upper);
        else                            % 分管需要混合
           % 压力随便取一条管路计算
           p_R_in_tube(1) = p_R_out(CV_num,tube_upper(1));% + dp_Uband(tube_upper(1));
           % 焓根据流量进行混合
           h_R_in_tube(1) = h_R_out(CV_num,tube_upper)*mdot_R(tube_upper)/sum(mdot_R(tube_upper));
        end

        mdot_R_tube = mdot_R(tube);
        mdot_MA_tube = mdot_MA(:,tube);

        % 空气侧暂时无混合
        T_MA_in_tube = T_MA_in(:,tube);
        T_MA_out_tube = T_MA_out(:,tube);

        p_MA_in_tube = p_MA_in(:,tube);
        p_MA_out_tube = p_MA_out(:,tube);

        % 为1则和空气控制体同向，为0则是反向
        flowDirection = TCinf.FlowDirection(tube);

        % 逐控制体计算
        % 严格来说，只要内部不报迭代未收敛的问题，就能可以以前后两次之间的差异作为迭代截至的判据
        for i3 = 1:CV_num
            if flowDirection == 1
                E_BD    = [h_R_out_tube(i3);T_MA_out_tube(i3)];
                x0      = [p_R_out_tube(i3);p_MA_out_tube(i3)];
                InletBD = [h_R_in_tube(i3);p_R_in_tube(i3);T_MA_in_tube(i3);p_MA_in_tube(i3);...
                            mdot_R_tube;mdot_MA_tube(i3)];
            else
                E_BD    = [h_R_out_tube(i3);T_MA_out_tube(CV_num + 1 - i3)];
                x0      = [p_R_out_tube(i3);p_MA_out_tube(CV_num + 1 - i3)];
                InletBD = [h_R_in_tube(i3);p_R_in_tube(i3);T_MA_in_tube(CV_num + 1 - i3);p_MA_in_tube(CV_num + 1 - i3);...
                            mdot_R_tube;mdot_MA_tube(CV_num + 1 - i3)];
            end
            
            flag = 2;
            N1 = ceil(tube/row) ;
            xout = alg(x0,E_BD,InletBD,GeoCondition,CV_num,N1,flag,Prop_handle);

            residual_max = max(max(abs(xout - x0)./x0),residual_max);

            p_R_out_tube(i3) = xout(1);
            p_R_out(i3,tube) = xout(1);

            if flowDirection == 1
               p_MA_out(i3,tube) = xout(2);
            else
               p_MA_out(CV_num + 1- i3,tube) = xout(2);
            end
            %% 这里添加一段U型弯压降计算段，绝热压降段，不考虑换热
            
            if i3 < CV_num
            p_R_in_tube(i3 + 1) = xout(1);
            end

            % if i3 == CV_num
            %     flag = 3;
            %     E_BD    = [h_R_out_tube(i3);h_R_in(1,tube_lower)];
            %     x0      = [p_R_out_tube(i3);p_R_in(1,tube_lower)];
            %     InletBD = [mdot_R_tube,Uband_length(tube)];
            %     % 获取下游的热力条件和压力条件，热力条件不变，获取压力条件
            %     [dp_Uband(tube),p_R_in(1,tube_lower)] = alg(x0,E_BD,InletBD,GeoCondition,CV_num,[],flag,Prop_handle);
            % else
            % 
            % end  
        end
        % 存储环路压降
        dp_tube(tube,1) = p_R_in_tube(1) -  p_R_out_tube(CV_num);
        % 更新进口压力
        p_R_in(:,tube) = p_R_in_tube;
        p_MA_in = [p_MA_inlet*ones(CV_num,row),p_MA_out(:,1:end-row)];
        
    end

    % 单次计算完毕，根据环路压降重新计算流量分配
    
    R_flow = dp_tube*1e6./(mdot_R.^2);
    u = fsolve(@(u)uF(u,R_flow,N,mdot0),u0,options);
    mdot_R = mdot0 + N*u;
    u0 = u;

    % 再更新一次压力
    for i2 = 1 : length_heat
        % 找到对应管的上游和下游，将边界条件、初始条件对应传入，随后逐控制体计算
        % 根据上游信息修改边界条件和初始条件
        tube = heatPaths(i2);

        p_R_in_tube  = p_R_in(:,tube);
        p_R_out_tube = p_R_out(:,tube);

        h_R_in_tube  = h_R_in(:,tube);
        h_R_out_tube = h_R_out(:,tube);
        
        % 工质侧更新上游混合情况
        tube_upper = predecessors_in{tube};
        if size(tube_upper,2) == 0      % 进口为总进口物性
           h_R_in_tube(1) = h_R_inlet;
           p_R_in_tube(1) = p_R_inlet;
        elseif size(tube_upper,2) == 1  % 内部管道直接传递
           h_R_in_tube(1) = h_R_out(CV_num,tube_upper);
           p_R_in_tube(1) = p_R_out(CV_num,tube_upper);
        else                            % 分管需要混合
           % 压力随便取一条管路计算
           p_R_in_tube(1) = p_R_out(CV_num,tube_upper(1));
           % 焓根据流量进行混合
           h_R_in_tube(1) = h_R_out(CV_num,tube_upper)*mdot_R(tube_upper)/sum(mdot_R(tube_upper));
        end

        mdot_R_tube = mdot_R(tube);
        mdot_MA_tube = mdot_MA(:,tube);

        % 空气侧暂时无混合
        T_MA_in_tube = T_MA_in(:,tube);
        T_MA_out_tube = T_MA_out(:,tube);

        p_MA_in_tube = p_MA_in(:,tube);
        p_MA_out_tube = p_MA_out(:,tube);

        % 为1则和空气控制体同向，为0则是反向
        flowDirection = TCinf.FlowDirection(tube);

        % 逐控制体计算
        % 严格来说，只要内部不报迭代未收敛的问题，就能可以以前后两次之间的差异作为迭代截至的判据
        for i3 = 1:CV_num
            if flowDirection == 1
                E_BD    = [h_R_out_tube(i3);T_MA_out_tube(i3)];
                x0      = [p_R_out_tube(i3);p_MA_out_tube(i3)];
                InletBD = [h_R_in_tube(i3);p_R_in_tube(i3);T_MA_in_tube(i3);p_MA_in_tube(i3);...
                            mdot_R_tube;mdot_MA_tube(i3)];
            else
                E_BD    = [h_R_out_tube(i3);T_MA_out_tube(CV_num + 1 - i3)];
                x0      = [p_R_out_tube(i3);p_MA_out_tube(CV_num + 1 - i3)];
                InletBD = [h_R_in_tube(i3);p_R_in_tube(i3);T_MA_in_tube(CV_num + 1 - i3);p_MA_in_tube(CV_num + 1 - i3);...
                            mdot_R_tube;mdot_MA_tube(CV_num + 1 - i3)];
            end
            
            flag = 2;
            N1 = ceil(tube/row) ;
            xout = alg(x0,E_BD,InletBD,GeoCondition,CV_num,N1,flag,Prop_handle);

            residual_max = max(max(abs(xout - x0)./x0),residual_max);

            p_R_out_tube(i3) = xout(1);
            p_R_out(i3,tube) = xout(1);

            if flowDirection == 1
               p_MA_out(i3,tube) = xout(2);
            else
               p_MA_out(CV_num + 1- i3,tube) = xout(2);
            end

            if i3 < CV_num
            p_R_in_tube(i3 + 1) = xout(1);
            end
        end
        % 存储环路压降
        dp_tube(tube,1) = p_R_in_tube(1) -  p_R_out_tube(CV_num);
        p_R_in(:,tube) = p_R_in_tube;
        p_MA_in = [p_MA_inlet*ones(CV_num,row),p_MA_out(:,1:end-row)];
    end

    % 第二次更新流量场
    R_flow = dp_tube*1e6./(mdot_R.^2);
    u = fsolve(@(u)uF(u,R_flow,N,mdot0),u0,options);
    mdot_R = mdot0 + N*u;

    % 这个才是正确的环路压降收敛判断准则——计算的最大环路压降差异<1Pa
    if max(abs((dp_tube')*N)) < 1e-6 && (residual_max < residual_limit)
        disp("压降收敛")
        break
    end
    u0 = u;
end
time = toc;

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