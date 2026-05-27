function  out = R_cal_10(x0,BD_R,GeoCondition,CV_num,Prop_handle,inlet_props,sat_cache)
L =         GeoCondition.L;                     % 单管长
D_inner =   GeoCondition.D_inner;               % 管内径
r =         GeoCondition.r;                     % 表面相对粗糙度
S =         pi*D_inner^2/4;                     % m^2 横截面积，用以计算速度

% 逐控制体的不动点迭代法
hin  =    BD_R.h_R_inlet;             % 进口焓            kJ/kg
mdot =    BD_R.mdot_R_inlet;          % 进口总质量流量    kg/s
pin  =    BD_R.p_R_inlet;             % 进口压力          MPa

hout = x0(1);
pout = x0(2);

h_CV = (hin + hout)/2;
p_CV = (pin + pout)/2;

% 物性调用（入口物性优先使用缓存，避免重复 REFPROP 查询）
if nargin >= 7 && ~isempty(inlet_props)
    Din_CV = inlet_props{1}; Prin_CV = inlet_props{2}; Nuin_CV = inlet_props{3};
    Tin_CV = inlet_props{4}; xin_CV = inlet_props{5}; kin_CV = inlet_props{6};
    vsatliq_in = inlet_props{7}; vsatvap_in = inlet_props{8};
    Prsatliq_in = inlet_props{9}; Prsatvap_in = inlet_props{10};
    Nusatliq_in = inlet_props{11}; Nusatvap_in = inlet_props{12};
    ksatliq_in = inlet_props{13}; ksatvap_in = inlet_props{14};
else
    [Din_CV, Prin_CV, Nuin_CV, Tin_CV, xin_CV, kin_CV, vsatliq_in,...
        vsatvap_in, Prsatliq_in, Prsatvap_in, Nusatliq_in, Nusatvap_in, ksatliq_in, ksatvap_in] = Prop1(pin,hin,Prop_handle);
end

if nargin >= 8 && ~isempty(sat_cache)
    [Dout_CV, Prout_CV, Nuout_CV, Tout_CV, xout_CV, kout_CV,vsatliq_out,...
        vsatvap_out, Prsatliq_out, Prsatvap_out, Nusatliq_out, Nusatvap_out, ksatliq_out, ksatvap_out] = Prop1(pout,hout,Prop_handle,sat_cache);
else
    [Dout_CV, Prout_CV, Nuout_CV, Tout_CV, xout_CV, kout_CV,vsatliq_out,...
        vsatvap_out, Prsatliq_out, Prsatvap_out, Nusatliq_out, Nusatvap_out, ksatliq_out, ksatvap_out] = Prop1(pout,hout,Prop_handle);
end

D_CV = (Din_CV + Dout_CV)/2;
Pr_CV = (Prin_CV + Prout_CV)/2;
Nu_CV = (Nuin_CV + Nuout_CV)/2;
T_CV = (Tin_CV + Tout_CV)/2;
x_CV = (xin_CV + xout_CV)/2;
k_CV = (kin_CV + kout_CV)/2;
if nargin >= 8 && ~isempty(sat_cache)
    vsatliq_CV = (vsatliq_in + sat_cache{1})/2;
    vsatvap_CV = (vsatvap_in + sat_cache{2})/2;
    Prsatliq_CV = (Prsatliq_in + sat_cache{3})/2;
    Prsatvap_CV = (Prsatvap_in + sat_cache{4})/2;
    Nusatliq_CV = (Nusatliq_in + sat_cache{5})/2;
    Nusatvap_CV = (Nusatvap_in + sat_cache{6})/2;
    ksatliq_CV = (ksatliq_in + sat_cache{7})/2;
    ksatvap_CV = (ksatvap_in + sat_cache{8})/2;
else
    vsatliq_CV = (vsatliq_in + vsatliq_out)/2;
    vsatvap_CV = (vsatvap_in + vsatvap_out)/2;
    Prsatliq_CV = (Prsatliq_in + Prsatliq_out)/2;
    Prsatvap_CV = (Prsatvap_in + Prsatvap_out)/2;
    Nusatliq_CV = (Nusatliq_in + Nusatliq_out)/2;
    Nusatvap_CV = (Nusatvap_in + Nusatvap_out)/2;
    ksatliq_CV = (ksatliq_in + ksatliq_out)/2;
    ksatvap_CV = (ksatvap_in + ksatvap_out)/2;
end

% 计算流速
veloctiy_CV = (mdot/D_CV)/S;

%% Correlation area 关联式区域

% 计算单相压降
% 使用单相流量，但是单相物性
Re_go = abs(mdot)*D_inner/(S*Nusatvap_CV*1e-6);
Re_lo = abs(mdot)*D_inner/(S*Nusatliq_CV*1e-6);

Re_lam = 2000;  Re_tur = 3000;

if Re_go <= Re_lam
    f_go = 64 / Re_go;
elseif Re_go > Re_tur
    f_go = 0.25*(log10(150.39/Re_go^0.98865-152.66/Re_go))^(-2);
else
    f_go = (1.1525*Re_go + 895)*1e-5;
end

if Re_lo <= Re_lam
    f_lo = 64 / Re_lo;
elseif Re_go > Re_tur
    f_lo = 0.25*(log10(150.39/Re_lo^0.98865-152.66/Re_lo))^(-2);
else
    f_lo = (1.1525*Re_lo + 895)*1e-5;
end

dpdL_lo = f_lo*(abs(mdot)/S)^2/(2*D_inner/vsatliq_CV);
dpdL_go = f_go*(abs(mdot)/S)^2/(2*D_inner/vsatvap_CV);

% Xu-Fang 2013 冷凝两相摩擦压降(NED 263, 87-96)
% 适用范围: R134a,R22,R410A 等, Dh 0.1–10mm, G 20–800, q 2–55.3
Y = sqrt(dpdL_go / dpdL_lo);
rho_tp = 1 / (x_CV*vsatvap_CV + (1-x_CV)*vsatliq_CV);  % 均相密度
g_acc = 9.81;
Fr_tp = (abs(mdot)/S)^2 / (g_acc * D_inner * rho_tp^2);
sigma = 0.008;  % N/m, R134a 冷凝近似值
We_tp = (abs(mdot)/S)^2 * D_inner / (rho_tp * sigma);
phi2_lo = Y^2 * x_CV^3 + (1 - x_CV^2.59)^0.632 * ...
    (1 + 2*x_CV^1.17*(Y^2 - 1) + 0.00775*x_CV^(-0.475)*Fr_tp^0.535*We_tp^0.188);
dp_f_CV = phi2_lo * dpdL_lo * L/CV_num;

%% Müller-Steinhagen and Heck (1986) — 保留备查
% G_msh = dpdL_lo + 2*(dpdL_go - dpdL_lo)*x_CV;
% dp_f_CV = L*(G_msh*(1-x_CV)^(1/3) + dpdL_go*x_CV^3);

%
% if Re_go > 3000
%     f_go = (-1.8*log10(6.9./Re_go+(r/3.7)^1.11)).^(-2);%0.25*(log(150.39/Re_go^0.98865 - 152.66/Re_go))^(-2);
% elseif Re_go < 2000
%     f_go = 64/Re_go;
% else
%     f_go = 1e-5*(1.1525*Re_go + 895);
% end
%
% if Re_lo > 3000
%     f_lo = 0.25*(log(150.39/(Re_lo^0.98865) - 152.66/Re_lo))^(-2);
% elseif Re_lo < 2000
%     f_lo = 64/Re_lo;
% else
%     f_lo =  (1e-5)*(1.1525*Re_lo + 895);
% end
% G = mdot/S;
% dpdL_lo = f_lo*G^2/(2*D_inner/vsatliq_CV);
% dpdL_go = f_go*G^2/(2*D_inner/vsatvap_CV);
%
% Y = sqrt(dpdL_go/dpdL_lo);
% %Th = (Y^2)*(x_CV^3) + ((1-x_CV)^(1/3))*(1+2*x_CV*(Y^2 -1 ));
% Th = (1 - x_CV)^(1/3) + 3.5 * x_CV^3 + x_CV^3 * (Y^2);
% %L = L + 50*D_inner;
% dp_f_CV = Th*dpdL_lo*L/CV_num;

%L = L + 30*D_inner;
%dp_f_CV =   (f_CV*L/CV_num).*(mdot.^2)./(2*D_CV*D_inner*S^2);    % Pa  摩擦压损
% dp_f_CV = L*(G*(1-x_CV)^(1/3)+ dpdL_go*x_CV^3);  % Müller-Steinhagen and Heck (1986) — 已替换为 Xu-Fang 2013
dp_v_CV =   16*mdot.^2/(pi^2*D_inner^4).*(1./Dout_CV - 1./Din_CV);   % Pa  速度压损
dp =    (dp_f_CV + dp_v_CV);


Re_1P = abs(veloctiy_CV)*D_inner/(Nu_CV*1e-6);
Re_2P = abs(mdot).*(1-x_CV + x_CV.*sqrt(vsatvap_CV./vsatliq_CV))*D_inner/S./(Nusatliq_CV*1e-6./vsatliq_CV);
f_CV =      (-1.8*log10(6.9./Re_1P+(r/3.7)^1.11)).^(-2);
% 单相换热 Gnielinski，层流-湍流 1000~1200 线性光滑过渡
Re_low  = 1000;
Re_high = 1200;
Nu_gn = (f_CV/8.*(Re_1P - 1000).*Pr_CV)./(1+12.7*sqrt(f_CV/8).*(Pr_CV.^(2/3)-1));
if Re_1P <= Re_low
    Nu_1P_CV = 3.66;
elseif Re_1P >= Re_high
    Nu_1P_CV = Nu_gn;
else
    w = (Re_1P - Re_low) / (Re_high - Re_low);
    Nu_1P_CV = 3.66 * (1 - w) + Nu_gn * w;
end
h_1P_CV = k_CV.*Nu_1P_CV/D_inner;

Nu_2P_CV = 0.05*((Re_2P).^0.8).*Prsatliq_CV.^0.33;
h_2P_CV = ksatliq_CV.*Nu_2P_CV/D_inner;

transition_range = 0.05;

% 液相到两相混合
w = min(max(x_CV / transition_range, 0), 1);   % 0→1 的权重
h_mix_1 = h_1P_CV .* (1 - w) + h_2P_CV .* w;  % 线性混合（可改为 Hermite）

% 两相到气相混合
w2 = min(max((x_CV - (1 - transition_range)) / transition_range, 0), 1);
h_R = h_mix_1 .* (1 - w2) + h_1P_CV .* w2;
dEF = mdot.*(hin - hout);

% if h_R == 0
%     pause()
% 
% end

%%
out = cell(7);
out{1} = Tin_CV;
out{2} = Tout_CV;
out{3} = h_R;
out{4} = dp;
out{5} = dEF;
out{6} = T_CV;
out{7} = {Dout_CV, Prout_CV, Nuout_CV, Tout_CV, xout_CV, kout_CV, vsatliq_out, ...
          vsatvap_out, Prsatliq_out, Prsatvap_out, Nusatliq_out, Nusatvap_out, ksatliq_out, ksatvap_out};

end