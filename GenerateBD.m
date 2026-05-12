function [BDCondition] = GenerateBD(GeoCondition,Prop_handle)

L        = GeoCondition.L;
L_fin    = GeoCondition.L_fin;

Area_air = L*L_fin;
Ra       = 287.047;

% ==================== 边界条件设置对话框 ====================
prompt = {'工质进口焓 h_R_inlet (kJ/kg):', ...
          '工质进口质量流量 mdot_R_inlet (kg/s):', ...
          '工质进口压力 p_R_inlet (MPa):', ...
          '空气进口温度 T_MA_inlet (K):', ...
          '空气进口压力 p_MA_inlet (MPa):', ...
          '空气进口风速 velocity_MA_in (m/s):', ...
          '空气进口相对湿度 RH_MA_inlet (0~1):'};
dlgtitle = '边界条件设置';
dims = [1 50];
definput = {'450', '0.02', '1', '300', '0.101325', '10', '0.5'};
answer = inputdlg(prompt, dlgtitle, dims, definput);

if isempty(answer)
    error('用户取消了边界条件设置');
end

h_R_inlet       = str2double(answer{1});    % kJ/kg   Inlet enthalpy
mdot_R_inlet    = str2double(answer{2});    % kg/s    Inlet mass flow rate
p_R_inlet       = str2double(answer{3});    % MPa     Inlet Pressure

T_MA_inlet      = str2double(answer{4});    % K       Inlet MA Temperature
p_MA_inlet      = str2double(answer{5});    % MPa     Inlet MA pressure
velocity_MA_in  = str2double(answer{6});    % m/s     Inlet MA velocity
RH_MA_inlet     = str2double(answer{7});    % 1       Relative Humidity

if any(isnan([h_R_inlet, mdot_R_inlet, p_R_inlet, T_MA_inlet, p_MA_inlet, velocity_MA_in, RH_MA_inlet]))
    error('边界条件输入无效，所有值必须为数字');
end
% ========================================================

Density_MA      = p_MA_inlet*1e6/Ra/T_MA_inlet;             % kg/m^3  MA密度
mdot_MA_inlet   = Density_MA*Area_air*velocity_MA_in;       % kg/s    MA质量流量
x_MA_inlet      = RHTox(RH_MA_inlet,T_MA_inlet,p_MA_inlet*1e6,Prop_handle);

% packaging Boundary structure  组装边界条件结构体
BD_MA = struct( ...
    "T_MA_inlet"     ,T_MA_inlet,...
    "mdot_MA_inlet"  ,mdot_MA_inlet,...
    "p_MA_inlet"     ,p_MA_inlet,...
    "x_MA_inlet"     ,x_MA_inlet);

BD_R = struct( ...
    "h_R_inlet"      ,h_R_inlet,...
    "mdot_R_inlet"   ,mdot_R_inlet,...
    "p_R_inlet"      ,p_R_inlet);

BDCondition = struct( ...
    "BD_R"           ,BD_R,...
    "BD_MA"          ,BD_MA);
end


function x = RHTox(RH,T,p,Prop_handle)
if any(RH>1,"all")
    warning('Relative Humidity>1!,相对湿度大于1!');
end
Ra = 287.047;  Rw = 461.523;         % J/K kg
psatvap = Prop_handle.psatvap;
p_w_sat = psatvap(T);  % Pa
p_w = p_w_sat.*RH;     % Pa
W = (Ra/Rw)*(p_w./(p - p_w));
x = W./(1+W);
end

function RH = xToRH(x,T,p,Prop_handle)
Ra = 287.047;  Rw = 461.523;         % J/K kg
psatvap = Prop_handle.psatvap;
p_w_sat = psatvap(T);  % Pa
p_w = p.*x*Rw./(x*Rw + (1-x)*Ra);
RH = p_w./p_w_sat;
if any(RH>1,"all")
    warning('Relative Humidity>1!,计算出了大于1的相对湿度！');
end
end