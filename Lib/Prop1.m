function [D,Pr,Nu,T,x,k,vsatliq,vsatvap,Prsatliq,Prsatvap,Nusatliq,Nusatvap,ksatliq,ksatvap] = Prop1(p,h,Prop_handle)

h_sat_liq = Prop_handle.h_sat_liq;
h_sat_vap = Prop_handle.h_sat_vap;
v_vap    = Prop_handle.v_vap;
v_liq    = Prop_handle.v_liq;
hmin     = Prop_handle.hmin;
hmax     = Prop_handle.hmax;
Pr_vap   = Prop_handle.Pr_vap;
Pr_liq   = Prop_handle.Pr_liq;
Nu_vap   = Prop_handle.Nu_vap;
Nu_liq   = Prop_handle.Nu_liq;
T_vap    = Prop_handle.T_vap;
T_liq    = Prop_handle.T_liq;
k_vap    = Prop_handle.k_vap;
k_liq    = Prop_handle.k_liq;

% 做成单输出模式，先检查相态

% 饱和参数计算
hsatliq = h_sat_liq(p);
hsatvap = h_sat_vap(p);

vsatliq = v_liq(0,p);
vsatvap = v_vap(1,p);
Prsatliq = Pr_liq(0,p);
Prsatvap = Pr_vap(1,p);
Nusatliq = Nu_liq(0,p);
Nusatvap = Nu_vap(1,p);
Tsatliq = T_liq(0,p);
ksatliq = k_liq(0,p);
ksatvap = k_vap(1,p);

if h < hsatliq
    hnorm = (h - hmin)/(hsatliq - hmin) - 1;

    v = v_liq(hnorm,p);
    Pr = Pr_liq(hnorm,p);
    Nu = Nu_liq(hnorm,p);
    T = T_liq(hnorm,p);
    k = k_liq(hnorm,p);
    x = 0;
elseif h > hsatvap
    hnorm = (h - hsatvap)./(hmax - hsatvap) + 1;

    v = v_vap(hnorm,p);
    Pr = Pr_vap(hnorm,p);
    Nu = Nu_vap(hnorm,p);
    T = T_vap(hnorm,p);
    k = k_vap(hnorm,p);
    x = 1;
else
    hnorm = (h - hsatliq)./(hsatvap - hsatliq);

    x = hnorm;
    v = vsatliq + hnorm.*(vsatvap - vsatliq);
    Pr = Prsatliq + hnorm.*(Prsatvap - Prsatliq);
    Nu = Nusatliq + hnorm.*(Nusatvap - Nusatliq);
    k = ksatliq + hnorm.*(ksatvap - ksatliq);
    T = Tsatliq;
end

D = 1./v;
end