function  out = R_cal_10(x0,BD_R,GeoCondition,CV_num,Prop_handle,flag)
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

% 物性调用
[Din_CV, Prin_CV, Nuin_CV, Tin_CV, xin_CV, kin_CV, vsatliq_in,...
    vsatvap_in, Prsatliq_in, Prsatvap_in, Nusatliq_in, Nusatvap_in, ksatliq_in, ksatvap_in] = Prop1(pin,hin,Prop_handle);

[Dout_CV, Prout_CV, Nuout_CV, Tout_CV, xout_CV, kout_CV,vsatliq_out,...
    vsatvap_out, Prsatliq_out, Prsatvap_out, Nusatliq_out, Nusatvap_out, ksatliq_out, ksatvap_out] = Prop1(pout,hout,Prop_handle);

D_CV = (Din_CV + Dout_CV)/2;
Pr_CV = (Prin_CV + Prout_CV)/2;
Nu_CV = (Nuin_CV + Nuout_CV)/2;
T_CV = (Tin_CV + Tout_CV)/2;
x_CV = (xin_CV + xout_CV)/2;
k_CV = (kin_CV + kout_CV)/2;
vsatliq_CV = (vsatliq_in + vsatliq_out)/2;
vsatvap_CV = (vsatvap_in + vsatvap_out)/2;
Prsatliq_CV = (Prsatliq_in + Prsatliq_out)/2;
Prsatvap_CV = (Prsatvap_in + Prsatvap_out)/2;
Nusatliq_CV = (Nusatliq_in + Nusatliq_out)/2;
Nusatvap_CV = (Nusatvap_in + Nusatvap_out)/2;
ksatliq_CV = (ksatliq_in + ksatliq_out)/2;
ksatvap_CV = (ksatvap_in + ksatvap_out)/2;

% 物性调用
veloctiy_CV = (mdot/D_CV)/S;

%% Correlation area 关联式区域

Re_1P = abs(veloctiy_CV)*D_inner./(Nu_CV*1e-6);
Re_2P = abs(mdot).*(1-x_CV + x_CV.*sqrt(vsatvap_CV./vsatliq_CV))*D_inner/S./(Nusatliq_CV*1e-6./vsatliq_CV);

f_CV =      (-1.8*log10(6.9./Re_1P+(r/3.7)^1.11)).^(-2);
% 计算压差

if flag == 2
Uband_L = BD_R.Uband_L;

% U型弯压降计算公式：
Re_1P = abs(veloctiy_CV)*D_inner./(Nu_CV*1e-6);
f_CV =      (-1.8*log10(6.9./Re_1P+(r/3.7)^1.11)).^(-2);

dp_f_CV =   (f_CV*Uband_L).*(mdot.^2)./(2*D_CV*D_inner*S^2);    % Pa  摩擦压损
dp_v_CV =   16*mdot.^2/(pi^2*D_inner^4).*(1./Dout_CV - 1./Din_CV);   % Pa  速度压损
dp =    (dp_f_CV + dp_v_CV);

out = dp;
return
end

% 计算单相压降
% 使用单相流量，但是单相物性
% Re_go = abs(mdot)*vsatvap_CV*D_inner./(S*Nusatvap_CV*1e-6);
% Re_lo = abs(mdot)*vsatliq_CV*D_inner./(S*Nusatliq_CV*1e-6);
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

L = L + 30*D_inner;
dp_f_CV =   (f_CV*L/CV_num).*(mdot.^2)./(2*D_CV*D_inner*S^2);    % Pa  摩擦压损
dp_v_CV =   16*mdot.^2/(pi^2*D_inner^4).*(1./Dout_CV - 1./Din_CV);   % Pa  速度压损
dp =    (dp_f_CV + dp_v_CV);

% 单相换热采用Gnielinski公式
Nu_1P_CV = (f_CV/8.*max(Re_1P - 1000,0).*Pr_CV)./(1+12.7*sqrt(f_CV/8).*(Pr_CV.^(2/3)-1));
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

%%
out{1} = Tin_CV;
out{2} = Tout_CV;
out{3} = h_R;
out{4} = dp;
out{5} = dEF;
out{6} = T_CV;

end