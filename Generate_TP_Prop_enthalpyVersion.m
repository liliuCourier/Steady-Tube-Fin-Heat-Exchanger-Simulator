function [TP_property] = Generate_TP_Prop_enthalpyVersion(libLoc,R,p_min,p_max,h_min,h_max,p_point,h_point_liq,h_point_vap)

% Example:

% libLoc = 'E:\refprop10\REFPROP';
% 
% % 输入参数
% R = 'R134a';
% p_point = 100;
% u_point_liq = 25;
% u_point_vap = 25;
% 
% u_min = 80;                  % kJ/kg
% u_max = 510;                 % kJ/kg
% p_min = 0.001;               % MPa
% p_max = 5.5;                 % MPa

% Generate_TP_Prop(libLoc,R,p_min,p_max,u_min,u_max,p_point,u_point_liq,u_point_vap)
% Generate_TP_Prop('E:\refprop10\REFPROP','R134a',0.001,5.5,80,510,100,25,25)
% Generate_TP_Prop('E:\refprop10\REFPROP','R290a',0.5,5,150,500,100,25,25)

% 计算三相点的各种值
%T_trip = getFluidProperty(libLoc, "T", "TRIP", 0, "", nan, R, 1, 1, 'MASS BASE SI');
p_trip = getFluidProperty(libLoc, "P", "TRIP", 0, "", nan, R, 1, 1, 'MASS BASE SI');
%u_trip = getFluidProperty(libLoc, "E", "TRIP", 0, "", nan, R, 1, 1, 'MASS BASE SI');

% 临界压力                 Pa
p_critcal_cal = getFluidProperty(libLoc, 'P','CRIT',1,'',123,R, 1, 1, 'MASS BASE SI');     % Pa
h_critcal_cal = getFluidProperty(libLoc, 'H','CRIT',1,'',123,R, 1, 1, 'MASS BASE SI');     % J/kg
p_critcal = p_critcal_cal/1e6;           % MPa

p_series = linspace(log(p_min*1e6),log(p_max*1e6),p_point);
p = exp(p_series);         % Pa
%p = logspace(log10(p_min*1e6),log10(p_max*1e6),p_point);              % Pa

% 压力上下限警告
assert(p(1) > p_trip, ...
    'GenerateProperty:MinPressure %c Greater Than Triple %c', [num2str(p(1)/1e6) ' MPa'], [num2str(p_trip/1e6), ' MPa'])
assert(p(1) < p_critcal_cal, ...
    'GenerateProperty:MinPressure %c Less Than Critical %c', [num2str(p(1)/1e6) ' MPa'], [num2str(p_critcal), ' MPa'])


% 接近临界点的热传输属性（默认）
% 封装变量名为 critical_region, 值为foundation.enum.critical_region.clip（削峰）或foundation.enum.critical_region.no_clip（不削峰）
critical_region = 'foundation.enum.critical_region.clip';

% 临界压力上下削波比例（默认）
% 封装变量名为 p_crit_fraction, 值为常规double
% 如果选了削峰，需要输入值
p_crit_fraction = 0.02;   % 默认削峰值

% 大气压（默认）           MPa
p_atm = 0.101325;
% 倒流的动压阈值（默认）   Pa
p_limit = 0.01;           % Pa

% 归一化液体内能向量
hnorm_liq = linspace(-1,0,h_point_liq);
hnorm_vap = linspace(1,2,h_point_vap);
hnorm = [hnorm_liq , hnorm_vap];

% 饱和液体比内能向量     J/kg
h_sat_liq = zeros(size(p,2),1);
h_sat_vap = zeros(size(p,2),1);
n_sub = p_point - sum(p>p_critcal_cal);

% v\s\T都是直接用归一化内能和压力进行插值
% v = zeros(size(unorm_liq,2) + size(unorm_vap,2),size(p,2));
% s = zeros(size(unorm_liq,2) + size(unorm_vap,2),size(p,2));
% T = zeros(size(unorm_liq,2) + size(unorm_vap,2),size(p,2));

% 以下三个输运参数需要分液相和气相
v_liq = zeros(size(hnorm_liq,2),size(p,2));  v_vap = zeros(size(hnorm_vap,2),size(p,2));
s_liq = zeros(size(hnorm_liq,2),size(p,2));  s_vap = zeros(size(hnorm_vap,2),size(p,2));
T_liq = zeros(size(hnorm_liq,2),size(p,2));  T_vap = zeros(size(hnorm_vap,2),size(p,2));
nu_liq = zeros(size(hnorm_liq,2),size(p,2)); nu_vap = zeros(size(hnorm_vap,2),size(p,2));
k_liq = zeros(size(hnorm_liq,2),size(p,2));  k_vap = zeros(size(hnorm_vap,2),size(p,2));
Pr_liq = zeros(size(hnorm_liq,2),size(p,2)); Pr_vap = zeros(size(hnorm_vap,2),size(p,2));

for j = 1: n_sub
    h_sat_liq(j) = getFluidProperty(libLoc, 'H','P',p(j),'Q',0,R, 1, 1, 'MASS BASE SI'); % J/kg
    h_sat_vap(j) = getFluidProperty(libLoc, 'H','P',p(j),'Q',1,R, 1, 1, 'MASS BASE SI'); % J/kg

    i = h_point_liq;
    [h_sat_liq(j), v_liq(i,j), s_liq(i,j), T_liq(i,j), nu_liq(i,j), k_liq(i,j), Pr_liq(i,j)] ...
        = getFluidProperty(libLoc, 'H,V,S,T,KV,TCX,PRANDTL','P',p(j),'Q',0,R, 1, 1, 'MASS BASE SI'); % J/kg
    i = 1;
    [h_sat_vap(j), v_vap(i,j), s_vap(i,j), T_vap(i,j), nu_vap(i,j), k_vap(i,j), Pr_vap(i,j)] ...
        = getFluidProperty(libLoc, 'H,V,S,T,KV,TCX,PRANDTL','P',p(j),'Q',1,R, 1, 1, 'MASS BASE SI'); % J/kg
end

%%

tolx = 1e-6;
opts = optimset('Display', 'off', 'MaxFunEvals', 200, 'TolX', tolx);
delta_h = max(h_sat_vap(1:n_sub) - h_sat_liq(1:n_sub));
h_last = h_critcal_cal;
for j = n_sub+1 : p_point
    % Midpoint of search interval
    h_mid = h_last;

    % Boundary of search interval, scaled to -1 and +1
    % Ensure that the search boundary is within uRange
    h_scaled_bnds = [-1 1];
    if h_scaled_bnds(1)*delta_h/2 + h_mid < h_min*1e3
        h_scaled_bnds(1) = (h_min*1e3 - h_mid)/(delta_h/2);
    end
    if h_scaled_bnds(2)*delta_h/2 + h_mid > h_max*1e3
        h_scaled_bnds(2) = (h_max*1e3 - h_mid)/(delta_h/2);
    end

    try
        % Search for specific internal energy at peak specific heat
        [h_scaled_check, ~, exitflag] = ...
            fminbnd(@(h_scaled)-getFluidProperty(libLoc, 'CP','P',p(j),'H',h_scaled*delta_h/2 + h_mid,R, 1, 1, 'MASS BASE SI'), ...
            h_scaled_bnds(1), h_scaled_bnds(2), opts);

        % If search fails, ignore this solution
        if exitflag(i) ~= 1
            h_scaled_check = nan;
        end
    catch
        % If search fails, ignore this solution
        %u_scaled_check = nan;
    end

    % Tolerance for checking if the solution is on interval boundary
    tol_bnd = 2*(tolx + 3*abs(h_scaled_check)*sqrt(eps));

    % Check if the solution is on the interval boundaries
    is_left  = h_scaled_check <= h_scaled_bnds(1) + tol_bnd;
    is_right = h_scaled_check >= h_scaled_bnds(2) - tol_bnd;

    % Check if the solution is in the interior
    is_interior = ~(is_left | is_right);

    % Solution is a peak only if it is in interior of interval
    % If peak is not found, then just extend upward from last point
    if is_interior
        h_sat_liq(j) = h_scaled_check*delta_h/2 + h_mid;
        h_sat_vap(j) = h_scaled_check*delta_h/2 + h_mid;
    else
        h_sat_liq(j:n) = h_last;
        h_sat_vap(j:n) = h_last;
        break
    end
    h_last = h_sat_liq(j);
end

%%

% 内能上下限警告
assert(h_min*1e3 < min(h_sat_liq(1:n_sub)), ...
    'GenerateProperty:MinInternalEnergy %c Less Than Saturation %c', [num2str(h_min) ' kJ/kg'], [num2str(min(h_sat_liq(1:n_sub))/1e3) ' kJ/kg'])
assert(h_max*1e3 > max(h_sat_vap(1:n_sub)), ...
    'GenerateProperty:MaxInternalEnergy %c Greater Than Saturation %c', [num2str(h_max) ' kJ/kg'], [num2str(max(h_sat_vap(1:n_sub))/1e3) ' kJ/kg'])

% 通过归一化内能反推内能向量：
% 实际液相区的比内能     J/kg
h_liq = h_min*1e3 + (hnorm_liq + 1).*(h_sat_liq - h_min*1e3) ;
h_vap = h_max*1e3 + (hnorm_vap - 2).*(h_max*1e3 - h_sat_vap) ;
h = [h_liq, h_vap];


% 归一化内能从-1~0 1~2的关于三个状态参数的插值表
for i = n_sub + 1 :p_point
    [v_liq(:,i),s_liq(:,i),T_liq(:,i),nu_liq(:,i),k_liq(:,i),Pr_liq(:,i)] = getFluidProperty(libLoc, 'V,S,T,KV,TCX,PRANDTL','H',h_liq(i,:),'P',p(i),R, 1, 1,'MASS BASE SI');
    [v_vap(:,i),s_vap(:,i),T_vap(:,i),nu_vap(:,i),k_vap(:,i),Pr_vap(:,i)] = getFluidProperty(libLoc, 'V,S,T,KV,TCX,PRANDTL','H',h_vap(i,:),'P',p(i),R, 1, 1,'MASS BASE SI');
    %[v(:,i),s(:,i),T(:,i)] = getFluidProperty(libLoc, 'V,S,T','E',u(i,:),'P',p(i),R, 1, 1,'MASS BASE SI');
end

for j = 1 : p_point
    i = 1 : h_point_liq-1;
    [v_liq(i,j), s_liq(i,j), T_liq(i,j), nu_liq(i,j), k_liq(i,j), Pr_liq(i,j)] ...
        = getFluidProperty(libLoc, 'V,S,T,KV,TCX,PRANDTL','H',h_liq(j,i),'P',p(j),R, 1, 1,'MASS BASE SI');

    i = 2 : h_point_vap;
    [v_vap(i,j), s_vap(i,j), T_vap(i,j), nu_vap(i,j), k_vap(i,j), Pr_vap(i,j)] ...
        = getFluidProperty(libLoc, 'V,S,T,KV,TCX,PRANDTL','H',h_vap(j,i),'P',p(j),R, 1, 1,'MASS BASE SI');

end

% 创建结构体防止复用
TP_property = struct();

% 输入结构体的具体值
TP_property.h_min = h_min;
TP_property.h_max = h_max;
TP_property.p_TLU = p/1e6;
TP_property.p_crit = p_critcal;
TP_property.critical_region = critical_region;
TP_property.p_crit_fraction = p_crit_fraction;
TP_property.p_atm = p_atm;
TP_property.q_rev = p_limit;

TP_property.hnorm_liq = hnorm_liq;
TP_property.v_liq = v_liq;
TP_property.s_liq = s_liq/1e3;

TP_property.T_liq = T_liq;
TP_property.nu_liq = nu_liq*1e4;

TP_property.k_liq = k_liq;
TP_property.Pr_liq = Pr_liq;
TP_property.h_sat_liq = h_sat_liq/1e3;
TP_property.hnorm_vap = hnorm_vap;

TP_property.v_vap = v_vap;
TP_property.s_vap = s_vap/1e3;
TP_property.T_vap = T_vap;
TP_property.nu_vap = nu_vap*1e4;
TP_property.k_vap = k_vap;
TP_property.Pr_vap = Pr_vap;
TP_property.h_sat_vap = h_sat_vap/1e3;



% fluidTables = twoPhaseFluidTables([80,510],[0.001,5.5],25,25,100,'R134a',libLoc);
% modelName = 'Example';
% twoPhaseFluidTables('Example/TP Table',fluidTables)

end