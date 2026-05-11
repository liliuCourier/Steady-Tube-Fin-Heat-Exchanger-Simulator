function [BDCondition] = GenerateBD(GeoCondition,Prop_handle)

L        = GeoCondition.L;
L_fin    = GeoCondition.L_fin;

Area_air = L*L_fin;
Ra       = 287.047;

% R
h_R_inlet       = 450;    % kJ/kg   Inlet enthalpy          进口焓
mdot_R_inlet    = 20e-3;  % kg/s    Inlet mass flow rate    进口总质量流量
p_R_inlet       = 1;      % MPa     Inlet Pressure          进口压力

% MA
T_MA_inlet      = 300;          % K    Inlet MA Temperature 进口温度
p_MA_inlet      = 0.101325;     % MPa  Inlet MA pressure    进口MA压力

velocity_MA_in  = 5;                                        % m/s     Inlet MA velocity MA进口速度,推荐使用速度
Density_MA      = p_MA_inlet*1e6/Ra/T_MA_inlet;             % kg/m^3  Inlet MA Density  MA进口密度，使用理想气体方程计算

mdot_MA_inlet   = Density_MA*Area_air*velocity_MA_in;       % kg/s    Inlet MA Mass flow rate  换算的MA进口总质量流量，
RH_MA_inlet     = 0.5;                                      % 1       Relative Humidity        MA进口的相对湿度
x_MA_inlet      = RHTox(RH_MA_inlet,T_MA_inlet,p_MA_inlet*1e6,Prop_handle);      % 函数要求输入的压力单位为Pa

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