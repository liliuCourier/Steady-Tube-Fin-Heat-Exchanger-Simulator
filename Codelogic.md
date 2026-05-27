# Code Logic — Steady Tube-Fin Heat Exchanger Simulator (Demo1.14)

## 零、稳态换热器求解策略

### 问题定义

给定管翅式换热器的流路拓扑、几何尺寸、工质和空气的入口边界条件（$T_{in}$, $p_{in}$, $\dot{m}$），求解稳态下：
- 每根管、每个控制容积的温度场、压力场、干度场
- 总换热量、总压降

### 为什么不能直接联立求解？

换热和压降是**双向耦合**的：换热量影响焓变→焓变影响干度→干度影响物性和压降→压降影响压力→压力影响饱和温度→影响换热温差。全联立（Jiang 2003 的非线性 NR+FDM）导致雅可比矩阵过大且需要有限差分。

### 本程序的求解策略：解耦 + 分步迭代

**核心思想**（沿袭 Ding/Liu 2004，本文改进）：把控制方程组拆成两个子问题交替求解。

### 外层迭代（Main.m 主循环）

```
1. 换热扫描（固定压力场，只算换热）
   沿 Kahn 拓扑排序遍历每根管、每段控制容积：
   ├─ 输入: 入口 h, p, mdot
   ├─ 调用 R_cal_10: 不动点迭代解 CV 出口状态
   │   └─ 单相: Gnielinski HTC → 能量平衡收敛
   │   └─ 两相: Cavallini-Zecchin HTC → 能量平衡收敛
   └─ 调用 DryA_cal_10: 空气侧 j/f 因子 → 换热量验证

2. 压力扫描（固定焓场，只算压降）
   同上遍历：
   ├─ 调用 R_cal_10: MSH 两相压降 / Haaland 单相
   ├─ 速度头压降 = ρv² 改变项
   └─ 输出 dp_CV

3. 流量重分配（仅当有环路时）
   环路压降不平衡 → Newton 迭代求解 N'·(R·m^e)=0
   前: mdot → 计算每管 dp → 环路残差 → 修正 mdot

4. 收敛检查
   环路 dp 不平衡 < 1e-6 Pa && 能量残差 < 1e-3 → 退出
```

### 为什么这个策略有效？

关键在于 $\Delta p / P \approx 10^{-3}$ 的量级关系——压力场对换热的影响（③ P→T_sat→ΔT）远弱于换热对压降的影响（② x→ρ,μ→摩擦）。

因此固定压力场算换热带来的误差量级约为 $O(\Delta p/P)$ 而非 $O(1)$，迭代 3-5 轮即可收敛。

### 两层收敛

程序有两层嵌套迭代，外层的有效性**建立在每个内层均成功收敛之上**：如果任何一次内层不动点迭代报"不收敛"，该次外层迭代计算的结果不可信。

#### 内层收敛（每个 CV 的独立计算）

| 迭代对象 | 变量 | 判据 | 上限 |
|---|---|---|---|
| 换热计算 (R_cal_10 + DryA_cal_10) | $h_{out}$, $T_{MA,out}$ | 能量残差相对误差 ≤ **1e-3** | 100 轮 |
| 压降计算 (R_cal_10) | $p_{out}$ | 动量压降残差 ≤ **1e-3** | 100 轮 |
| 弯管压降 (bend_cal) | $p_{out}$ | 不动点 |p_new-p_old| ≤ **1e-6** | 15/100 轮 |

内层任意 CV 的不动点迭代失败（达到上限仍未满足残差）→ 本次外层迭代的所有 CV 结果作废，程序应中止或标记不可信。

#### 外层收敛（整体换热器的迭代）

前提：当前外层迭代中**所有 CV 的内层均已收敛**。

| 判据 | 阈值 |
|---|---|
| **环路压降平衡**：$|(dp_{tube})' \cdot N|$ < **1e-6** Pa | 各环路进出口压降差异 < 1 Pa |
| **前后迭代一致性**：$residual_{max} <$ **1e-3** | 前后两轮外层迭代所有状态量（h, p, T_MA）的最大相对变化 < 1e-3 |
| **上限**：10 轮 | 达到上限仍未收敛 → 中止 |

两条判据**必须同时满足**。外层迭代不使用固定的最大迭代次数作为收敛判据——达到 10 轮仍未收敛则视为求解失败。

### 数据结构层次

```
管进口状态 ──→ CV1入口 ──→ CV1出口 ──→ ... ──→ CVn出口 ──→ 管出口状态
  (tube_inlet)  (h_in,p_in)  (h_out,p_out)       (h_out,p_out)  (tube_outlet)
      └──── 入口弯管 ───────────────────────── 出口弯管 ────┘
```

**弯管处理的设计思路**（Demo1.14 核心抽象）：

单根管在传统意义上仅包含若干控制容积（CV1...CVn）。为处理 U 型弯，本程序将"管"的概念做了一个关键延伸——从只包含控制体扩展到**同时包含其附属的 U 型弯**：

- **串联时**：管包含其下游的 U 型弯（弯管归上游管，作为 tube_outlet 的一部分）
- **分流时**：下游管包含来自分流点的入口 U 型弯（弯管归下游管，在 tube_inlet → CV1 之间）
- **汇流时**：上游管包含去往汇流点的出口 U 型弯（弯管归上游管，在 CVn → tube_outlet 之间）

这样设计的原因是：一根管在两个方向上有且只有一个 U 型弯需要处理（入口或出口），如果将弯管单独列为独立计算单元，则每个 U 型弯会同时属于两个管——入口弯管既属于上游管的出口、又属于下游管的入口——造成归属混乱。将弯管附着到管上，每个弯管只有一个归属，计算路径清晰：管级变量在弯管之前（入口）和之后（出口）做了一步压力变换，CV 循环完全不需要感知弯管的存在。

**最关键的简化**：串联时出口弯管归属上游管的规则，使得所有连接情况（串联、分流、汇流）都能统一处理，无需在流路拓扑层面对弯管做特殊分支判断。

管间传递通过 `tube_inlet(tube) = tube_outlet(upstream)` 实现（串联等焓等压，汇流焓加权平均）。

---

## 目录结构

```
Program_1_AI/
├── PreProcessing.m              ← 入口1: 预处理（流路/几何/边界）
├── Main.m                       ← 入口2: 稳态求解
├── Main_loop_base.m              ← 入口2b: 两阶段环路优先求解器
├── PostProcessing.m             ← 入口3: 后处理（自动触发）
├── readme.md                    ← 版本记录
├── Codelogic.md                 ← 本文档
├── single_phase_pressure_drop_summary.md   ← 单相压降公式总结
├── two_phase_pressure_drop_summary.md      ← 两相压降公式总结
├── .gitignore
├── PreProc/                     ← 预处理模块
│   ├── HX_Path_Planner1.m       GUI 流路设计器
│   ├── GenerateGeo.m            几何参数设置（含 bend_K 常量 1.5）
│   └── GenerateBD.m             边界条件与求解设置
├── Solver/                      ← 求解器核心
│   ├── mdot_Initial.m           流量初始化（线性规划 + 零空间）
│   ├── buildPath.m              广度优先/环路优先路径生成
│   ├── R_cal_10.m               工质侧单控制容积求解
│   ├── DryA_cal_10.m            空气侧单控制容积求解
│   └── Generate_TP_Prop_enthalpyVersion.m  工质物性表生成
├── Lib/                         ← 物性库
│   ├── Prop_load.m              REFPROP 物性插值句柄加载
│   └── Prop1.m                  单点工质物性查询
└── PostProcessing/              ← 后处理模块
    ├── plotAlongPath.m          沿线物性分布曲线
    ├── plotLoopBalance.m        环路压降平衡验证
    ├── summaryTable.m           控制台性能汇总表
    └── exportHxPerf.m           标准化数据导出
```

---

## 一、顶层工作流

```
PreProcessing ──→ Main ──→ PostProcessing (自动)
     │                │
     ▼                ▼
  生成 TCinf        求解热力场/压力场
  生成 GeoCondition  收敛判断
  生成 BDCondition   自动触发后处理
```

### PreProcessing.m（脚本）

**职责：** 依次完成三类预处理，将结果持久化到 base workspace。

**调用链：**
```
PreProcessing
  ├─ Prop_load()                     → Prop_handle (物性句柄)
  ├─ HX_Path_Planner1()              → TCinf (用户交互导出)
  ├─ GenerateGeo(TCinf)              → GeoCondition (对话框)
  └─ GenerateBD(GeoCondition, Prop_handle) → BDCondition (对话框)
```

**依赖路径：** `PreProc/` `Lib/` `Solver/`

---

### Main.m（脚本 + 3 个嵌套子函数）

**职责：** 换热与压降解耦迭代求解。

**主循环（每轮迭代）：**

```
有环路时:                          无环路时:
  ① 广度优先换热扫描                 ① 广度优先换热扫描
  ② 广度优先压力更新                 ② 残差收敛 → 退出
  ③ Newton 流量重分配
  ④ 环路压降收敛 → 退出
```

**调用链：**
```
Main
  ├─ mdot_Initial(BDCondition, GeoCondition, TCinf)
  │     → linprog()  正流量可行解
  │     → null(A1)   零空间基向量 N
  │     输出: N, u0, mdot0, mdot_R_init
  │
  ├─ buildPath(TCinf.TC_matrix, N)
  │     输出: heatPaths, pdropPaths, predecessors_in, predecessors_out
  │
  ├─ scanTubes(heatPaths, ...)        [嵌套子函数]
  │     ├─ upstream passing: tube_inlet(tube) = tube_outlet(upstream)
  │     ├─ bend_in check → bend_cal() 入口弯管压降
  │     ├─ CV loop: alg(x0, BD, ..., solver_flag, Prop_handle)  [嵌套子函数]
  │     │     ├─ solver_flag=1: R_cal_10() + DryA_cal_10()  换热
  │     │     ├─ solver_flag=2: R_cal_10()                   压力
  │     │     └─ solver_flag=3: 等焓不动点                   弯管压降
  │     ├─ CV loop end: tube_outlet(tube) = CV_n outlet
  │     └─ bend_out check → bend_cal() 出口弯管压降
  │
  ├─ Newton 流量更新
  │
  └─ PostProcessing                   自动后处理
```

**嵌套子函数：**

| 函数 | 签名 | 职责 |
|------|------|------|
| `scanTubes` | `function [..., h_R_tube_outlet, p_R_tube_outlet] = scanTubes(paths, pre_in, pre_out, ..., h_R_tube_inlet, p_R_tube_inlet, h_R_tube_outlet, p_R_tube_outlet, ...)` | 沿指定管序遍历，逐管逐 CV 调用 `alg()`，管级变量贯穿 |
| `alg` | `function [F, cache_R_out, cache_MA_out] = alg(x0, BD, InletBD, GeoCondition, CV, N, solver_flag, Prop_handle, cache_R_in, cache_MA_in)` | 单 CV 不动点迭代：solver_flag=1 换热量，=2 压降，=3 弯管压降 |
| `bend_cal` | `function dp_out = bend_cal(h_in, p_in, mdot, L_geom, Geo, Ph)` | 弯管压降（MSH沿程+L_equiv），返回 MPa |
| `bend_len` | `function Lb = bend_len(k1, k2, Gc)` | 管间几何长度：弧长+直线段，利用 P_row/P_col 推导 |

**依赖路径：** `Lib/` `Solver/` `PostProcessing/` `PreProc/`

---

### Main_loop_base.m（脚本 + 6 个嵌套子函数）

**职责：** 两阶段环路优先求解器，旨在通过 Phase1 快速稳定流量分布来减少总迭代数。

**策略：**

```
Phase 1: 恒压热力扫描 + 流量更新（不更新压力场）
  ├─ 全局饱和物性缓存 sat_global（p = p_R_inlet，所有 CV/迭代共用）
  ├─ scanTubes_phase1（仅热力扫描，提取压阻 dp_CV，含 bend 压阻）
  ├─ Newton 流量更新（同 Main）
  └─ 步出判据：前后两次流量相对变化 < 1e-2
       ↓
Phase 2: 完整热力+压力+流量更新（与 Main 一致）
  ├─ scanTubes(solver_flag=1)  热力扫描（含 bend）
  ├─ scanTubes(solver_flag=2)  压力扫描（含 bend）
  ├─ Newton 流量更新
  └─ 收敛判据：同 Main
```

**嵌套子函数：**

| 函数 | 职责 |
|------|------|
| `scanTubes_phase1` | Phase1 专用：仅热力扫描 + 提取各 CV 压阻，压力场不变，全局饱和缓存，含弯管压阻 |
| `alg_phase1` | 与 `alg(solver_flag=1)` 等价，额外返回 dp_CV (Pa)，复用 sat_global；solver_flag=3 为弯管压降 |
| `scanTubes` | Phase2 使用，与 Main 完全一致 |
| `alg` | Phase2 使用，含 solver_flag=3 弯管压降 |
| `bend_cal_phase1` | Phase1 弯管压降（MSH沿程+L_equiv），返回 Pa |
| `bend_cal` | Phase2 弯管压降，返回 MPa |

**与 Main.m 的关键差异：**

| 特性 | Main | Main_loop_base Phase1 |
|------|------|----------------------|
| 压力场 | 每轮迭代更新 | 全程恒定 = p_R_inlet |
| 饱和物性 | 每 CV 独立计算 (p_out) | 全局缓存 (p_inlet)，11 项 |
| 扫描内容 | 热力 + 压力（2 次 scanTubes/iter） | 仅热力（1 次 scanTubes_phase1/iter） |
| 弯管压降 | bend_cal (MPa) | bend_cal_phase1 (Pa，仅压阻) |

---

### PostProcessing.m（脚本）

**职责：** 一键生成全部可视化与性能报告。Main.m 求解完成后自动调用。

**调用链：**
```
PostProcessing
  ├─ plotAlongPath(...)     → Fig: 2×2 沿线物性分布
  ├─ plotLoopBalance(...)   → Fig: 环路压降平衡柱状图
  ├─ summaryTable(...)      → 控制台性能汇总表
  ├─ 收敛历史绘图           → Fig: 能量残差 + 环路残差收敛曲线
  └─ exportHxPerf(...)      → hxPerf 结构体
```

**依赖路径：** `PostProcessing/` `Lib/` `Solver/`

---

## 二、核心数据结构

### 2.1 管级变量（Demo1.14 新增）

| 变量 | 维度 | 物理含义 |
|------|------|----------|
| `h_R_tube_inlet` | 1×Tube_num | 管进口焓（弯管之前） |
| `p_R_tube_inlet` | 1×Tube_num | 管进口压力（弯管之前） |
| `h_R_tube_outlet` | 1×Tube_num | 管出口焓（弯管之后） |
| `p_R_tube_outlet` | 1×Tube_num | 管出口压力（弯管之后） |

**传递规则**（scanTubes 开头）：
```
入口管（无上游）: tube_inlet = 集管入口值
串联（单上游）  : tube_inlet = upstream.tube_outlet
汇流（多上游）  : tube_inlet = weighted_avg(upstream.tube_outlet) by mass flow
```

**管内 CV 链接**：
```
tube_inlet ─→ [bend_in?] ─→ CV₁ ─→ ... ─→ CV_n ─→ [bend_out?] ─→ tube_outlet
```

弯管存在时：CV1 入口压力 = tube_inlet.p - dp_bend_in；CVn 出口后 tube_outlet.p = CVn.p_out - dp_bend_out。焓不变（等焓）。

### 2.2 结构体

| 结构体 | 产生于 | 关键字段 | 消耗于 |
|--------|--------|----------|--------|
| `TCinf` | `HX_Path_Planner1` → `cbExport` | `TC_matrix` `inlet_num` `outlet_num` `FlowDirection` `row` `col` `con_num` `Tube_num` | `GenerateGeo` `mdot_Initial` `buildPath` `plotAlongPath` `summaryTable` `exportHxPerf` |
| `GeoCondition` | `GenerateGeo(TCinf)` | `L` `D_inner` `D_outer` `r` `A_R` `A_MA` `P_row` `P_col` `Fin_pitch` `dx_fin` `Tube_num` `row` `col` | `GenerateBD` `R_cal_10` `DryA_cal_10` `bend_len` `summaryTable` `exportHxPerf` |
| `BDCondition` | `GenerateBD(GeoCondition, Prop_handle)` | `BD_R( h_R_inlet mdot_R_inlet p_R_inlet )` `BD_MA( T_MA_inlet mdot_MA_inlet p_MA_inlet )` `CV_num` | `mdot_Initial` `R_cal_10` `DryA_cal_10` `summaryTable` `exportHxPerf` |
| `Prop_handle` | `Prop_load(libLoc, R, ...)` | ~26 个 `griddedInterpolant` 函数句柄（工质+湿空气） | `Prop1` `GenerateBD` `DryA_cal_10` `R_cal_10` 以及各后处理函数 |
| `N` | `mdot_Initial(...)` | Tube_num × n_loops 矩阵，零空间基向量 | `Main` `buildPath` `plotLoopBalance` `summaryTable` `exportHxPerf` |
| `hxPerf` | `exportHxPerf(...)` | `Geo` `Topo` `Boundary` `Flow` `Thermal` `Performance` | 优化器 |

---

## 三、物性体系

```
Prop_load (REFPROP 初始化)
  ├─ Generate_TP_Prop_enthalpyVersion  ← REFPROP 调用 getFluidProperty()
  │     构建 T-P 性质查找表 (h_sat, v, T, Nu, k, Pr 等)
  ├─ griddedInterpolant()              ← 1D/2D 插值器
  └─ 湿空气物性表 (hair, visair, kair, Prair, cpair, psatvap, ...)
       ↓
Prop1(p, h, Prop_handle)              ← 单点工质查询
  输入: 压力 p (MPa), 焓 h (kJ/kg)
  输出: [D, Pr, Nu, T, x, k, ...]   ← 根据相态 (过冷/两相/过热) 分支处理
```

**关联式：**
- 工质侧单相摩擦系数: Haaland / Blasius / Fang（全 Re 范围，层流→过渡→湍流）
- 工质侧单相换热: Gnielinski（层流-湍流线性过渡）
- 工质侧两相换热: Cavallini and Zecchin（冷凝）
- 工质侧两相压降: Müller-Steinhagen and Heck (MSH)
- 弯管压降: MSH沿程(L_geom + L_equiv) + Paliwoda 两相修正，L_equiv = K·D/(2f)，K=1.5
- 空气侧: Wang-Chi-Chang Plate-Fin (j 因子 + f 因子)

---

## 四、求解收敛判据

| 场景 | 判据 |
|------|------|
| 无环路 | `residual_max < 1e-3`（能量残差） |
| 有环路 | `max(abs((dp_tube')*N)) < 1e-6` **且** `residual_max < 1e-3` |
| 外层迭代上限 | 10 轮 |
| CV 内层不动点迭代上限 | 100 轮 |
| 弯管内层不动点迭代上限 | 15 轮（bend_cal_phase1）/ 100 轮（bend_cal），容差 1e-6 |

---

## 五、调用依赖矩阵

```
                    PreProc/     Solver/       Lib/          PostProc/
调用者              HX_P  GenG  GenB  mdot  bPath  R_cal  DryA  GenTP  PLoad  P1   plotAP  plotLB  sumTbl  expHxP
────────────────────────────────────────────────────────────────────────────────────────────────────────────────
PreProcessing.m      ✓     ✓     ✓                                          ✓
Main.m                              ✓     ✓      ✓     ✓                   (✓)   (✓)    ✓
  └ scanTubes/alg                                      ✓     ✓             (✓)
PostProcessing.m                                                           (✓)    ✓       ✓       ✓       ✓
GenerateGeo.m                                                                     (TCinf 入参)
GenerateBD.m                                                                      ✓
Prop_load.m                                                                            ✓
Prop1.m                                                                                          ✓
plotAlongPath.m                                                                           ✓
summaryTable.m                                                                            ✓
exportHxPerf.m                                                                            ✓
```

`(✓)` = 间接/条件调用，`✓` = 直接调用。

---

## 六、注意事项

1. **`buildPath.m` 函数名不一致：** 文件名为 `buildPath.m`，调用为 `buildPath(...)`，但内部声明为 `function [...] = buildPaths(A, N)`（带 s）。Windows 不区分大小写不报错，跨平台需注意。

2. **`Main.m` / `PreProcessing.m` / `PostProcessing.m` 均为脚本：** 不是函数，直接在 base workspace 操作变量，无显式参数传递。`Main.m` 内含 4 个嵌套子函数（`alg` `scanTubes` `bend_cal` `bend_len`）。

3. **REFPROP 外部依赖：** 整个物性体系依赖 `getFluidProperty()`，需 REFPROP 安装并位于 MATLAB 路径中。

4. **压力单位约定：** 所有 `p_*` 变量为 MPa，`dp_tube` 为 MPa（显示时 ×1e6 转 Pa），`bend_cal_phase1` 返回 Pa，`bend_cal` 返回 MPa。

5. **HX_Path_Planner1 窗口持久化：** 点击关闭按钮隐藏而非销毁，句柄保存在 `hxDesigner`，`figure(hxDesigner)` 可随时唤出。

---

## 七、已知问题与经验

### 7.1 流量初始化解与压降负值问题

**症状**：管压降 `dp_tube` 计算为负值，导致 `R_flow` 和 `dp_Pa` 为负，流量重分配牛顿迭代发散，最终 DryA_cal_10 因温度变复数崩溃。

**根因链路**：

```
mdot_Initial 用 linprog 作流量初始化（epsilon = mdot_inlet/4）
  → 部分支路分配到大流量、部分支路分配到极小流量
  → Xu-Fang 2013 冷凝关联式的摩擦压降 dp_f ∝ G²·L_tot
  → 极小流量管的 dp_f 可低至个位数 Pa
  → 速度压损 dp_v = 16·mdot²/(π²D⁴)·(1/ρ_out − 1/ρ_in)
  → 当工质压力降低同时密度减小（冷凝段：液相→气相膨胀），ρ_out < ρ_in
  → (1/ρ_out − 1/ρ_in) > 0 → dp_v > 0
  → dp_v > dp_f → 总压降 dp = dp_f + dp_v 为负
  → 管出口压力 > 管入口压力（物理上不可能）
  → 求解器发散
```

**修复**（mdot_Initial.m）：调整初始流量场的下界限制，确保每根管的最小流量足以让摩擦压降主导速度压损。具体做法是将 `epsilon` 从 `mdot_R_inlet/4` 收紧或改为按管数均分的初始值，让所有支路的初始流量在一个量级上。
