function out = DryA_cal_10(x0,BD_MA,GeoCondition,CV_num,N1,Prop_handle,inlet_air_props)
% 取出必要的集合条件，由于翅片的存在，需要的几何条件特别多
D_outer =   GeoCondition.D_outer;       % 管外径
L =         GeoCondition.L;             % 管长
row =       GeoCondition.row;             % 管长

Fin_pitch = GeoCondition.Fin_pitch;     % 翅片间距
H_fin =     GeoCondition.H_fin;         % 翅片宽、深——平行风流动方向——longi length
L_fin =     GeoCondition.L_fin;         % 翅片长度——垂直于风流动方向——vertical length
P_row =     GeoCondition.P_row;         % 管间距，vertical spacing
P_col =     GeoCondition.P_col;         % 管间距，horizontal spacing
dx_fin =    GeoCondition.dx_fin;        % 翅片厚度
A_MA =      GeoCondition.A_MA;          % 空气侧的总换热面积

% 固有物性，空气和水的气体常数
Ra = 287.047;

% 不同于工质侧的物性调用，空气侧的物性调用比较杂乱，这块有优化空间
hair = Prop_handle.hair;
visair =Prop_handle.visair;
kair =Prop_handle.kair;
Prair =Prop_handle.Prair;
cpair = Prop_handle.cpair;

% 取出边界条件
T_in =          BD_MA.T_MA_inlet;         % K    进口温度
mdot    =       BD_MA.mdot_MA_inlet;      % kg/s 列向量，进口质量流量
p_inlet =       BD_MA.p_MA_inlet;         % MPa  进口压力 

pin =       p_inlet;
Tin =       T_in;
vin =   Ra*Tin./(pin*1e6);

% 取出初始条件
Tout =          x0(1);
pout =          x0(2);

%% 计算中间变量，以及必要的物性计算
% 空气物性（入口优先使用缓存，避免重复插值查询）
if nargin >= 7 && ~isempty(inlet_air_props)
    h_air_in   = inlet_air_props{1};
    Pr_air_in  = inlet_air_props{2};
    vis_air_in = inlet_air_props{3};
    k_air_in   = inlet_air_props{4};
    cp_air_in  = inlet_air_props{5};
    vin        = inlet_air_props{6};
else
    h_air_in   = hair(Tin);
    Pr_air_in  = Prair(Tin);
    vis_air_in = visair(Tin);
    k_air_in   = kair(Tin);
    cp_air_in  = cpair(Tin);
    vin        = Ra*Tin./(pin*1e6);
end
h_air_out  = hair(Tout);
Pr_air_out = Prair(Tout);
vis_air_out = visair(Tout);
k_air_out  = kair(Tout);
cp_air_out = cpair(Tout);



% 使用理想气体来求空气的密度
p_CV =  (pin + pout)/2;                 % MPa     控制体压力
T_CV =  (Tin + Tout)/2;                 % K       控制体温度
v_CV =  Ra*T_CV./(p_CV*1e6);            % m^3/kg  控制体比容，采用理想气体方程计算

Pr_CV  =        (Pr_air_in + Pr_air_out)/2;
vis_CV =        (vis_air_in + vis_air_out)/2;
k_CV   =        (k_air_in + k_air_out)/2;
cp_CV   =       (cp_air_in + cp_air_out)/2;

% 出口比容
vout = Ra*Tout./(pout*1e6);

%% Correlation area关联式区域
A_total = L_fin*L;                                                  % m^2   在发生截面收缩之前的总通流面积
velocity_in = abs(mdot).*vin/(A_total/row/CV_num);
velocity_out = abs(mdot).*vout/(A_total/row/CV_num);                  % m/s   对应的进口速度
scale = (P_row/(P_row-D_outer))*(Fin_pitch/(Fin_pitch - dx_fin));   % 1     截面收缩比的倒数
velocity_max = velocity_in*scale;

N = N1;

Dh = 4 * H_fin * A_total/scale/A_MA;
Re = D_outer*velocity_max./(v_CV.*vis_CV);

% 这里要对Re做出限制，不允许很低的雷诺数，但是允许特别高的雷诺数
if any(Re-100<0)
%warning("雷诺数过低，脱离该关联式应用范围")
end

Re = max(Re,100);

P1 = 1.9 - 0.23*log(Re);
P2 = -0.236+0.126*log(Re);
P3 = -0.361 - 0.042*N./log(Re) + 0.158*log(N*(Fin_pitch/D_outer)^0.41);
P4 = -1.224 - 0.076*(P_col/Dh)^1.42./log(Re);
P5 = -0.083 + 0.058*N./log(Re);
P6 = -5.735 + 1.21*log(Re./N);
F1 = -0.764 + 0.739*P_row/P_col + 0.177*Fin_pitch/D_outer - 0.00758./N;
F2 = -15.689 + 64.021./log(Re);
F3 = 1.696 - 15.695./log(Re);

j1 = 0.108*(Re.^(-0.29)).*((P_row/P_col).^P1)*(Fin_pitch/D_outer)^(-1.084)*(Fin_pitch/Dh)^(-0.786).*(Fin_pitch/P_row).^P2;
j2 = 0.086*(Re.^P3).*(N.^P4).*((Fin_pitch/D_outer).^(P5)).*((Fin_pitch/Dh).^(P6))*(Fin_pitch/P_row)^(-0.93);

f = 0.0267*(Re.^F1).*((P_row/P_col).^F2).*((Fin_pitch/D_outer).^F3);
if N1 == 1
    j = j1;
else
    j = j2;
end

h = j.*(1./v_CV).*velocity_max.*cp_CV.*Pr_CV.^(-2/3);                                       
dp_f = f.*(A_MA/(A_total/scale)).*(velocity_max.^2./v_CV/2);
dp = dp_f ;
F_dp = dp_f - (pin - pout)*1e6;

h_in  = h_air_in;
h_out = h_air_out;

H = (P_col - D_outer)/2 + dx_fin/2;
mH = H*sqrt(2*h/220/dx_fin);
n_fin = tanh(mH)./mH;

%% 残差构造

dEF_air = mdot.*(h_in - h_out)  ; 


        out{1} = dEF_air;
        out{2} = n_fin;
        out{3} = h;
        out{4} = cp_CV;
        out{5} = dp;
        out{6} = F_dp;
        out{7} = {h_air_out, Pr_air_out, vis_air_out, k_air_out, cp_air_out, vout};

end