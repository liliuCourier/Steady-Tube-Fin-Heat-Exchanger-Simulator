% 简单介绍一下这个函数，用来获取工质侧的初始条件，但是大幅减少了变量个数
% 根据节点的连接关系——管道连接矩阵代表了管道流量之间的线性关系——本身的自由变量个数也是有限制的
% 以本程序常用的TCMatric为例，通过管道连接矩阵+进口流量分配方程进行增广，构成的非线性方程约束中
% 真正的自由变量只有一个！
% 也就是说，在质量流量方程中，从原来的求解逻辑中变量为Tube_num个变为了1个！减少了Tube_num - 1个变量

function [N,u0,mdot0,mdot_R_init,exitflag] = mdot_Initial(BDCondition,GeoCondition,TCinf)

% exitflag 输出不为1和4的时候需要调整总初值和管路

% 获取优良的质量流量初始条件
% 从对应的结构体重获取需要的参数
TC_matrix = TCinf.TC_matrix;
inlet_num = TCinf.inlet_num;
con_num   = TCinf.con_num;

Tube_num = GeoCondition.Tube_num;

mdot_R_inlet = BDCondition.BD_R.mdot_R_inlet;

% 构建非线性约束求解流程
b0 = zeros(1,Tube_num);
b0(inlet_num) = 1;
A1 = [TC_matrix;b0];
b = [zeros(con_num,1);mdot_R_inlet];

% 假设已构建好 A1 (r x n) 和 b (r x 1)
epsilon =  mdot_R_inlet / Tube_num;   % 极小正数下限

% 构造线性规划：无目标函数（或设零目标），只需找可行解
f = zeros(Tube_num, 1);                  % 目标函数系数全零
Aeq = A1;
beq = b;
lb = epsilon * ones(Tube_num, 1);        % 下界
ub = [];                          % 无上界

options = optimoptions('linprog', 'Display', 'off');
[mdot_R_init, ~, exitflag] = linprog(f, [], [], Aeq, beq, lb, ub, options);

if exitflag ~= 1
    error('无法找到正流量初始解，请检查管网约束是否矛盾。');
end

% 似乎没有什么好办法，暂且这样处理吧,设置一个好的初值是大前提
% 零空间的列向量，列向量的个数代表了自由变量的个数
% 反过来输出自由变量的初值
mdot0 = A1\b;  % 最小范数初值
N = null(A1,'rational');
u0 = N\(mdot_R_init - mdot0);

end