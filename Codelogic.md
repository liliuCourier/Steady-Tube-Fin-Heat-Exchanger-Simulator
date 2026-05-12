# Code Logic — Steady Tube-Fin Heat Exchanger Simulator (Demo1.02)

## 目录结构

```
Program_1_AI/
├── PreProcessing.m              ← 入口1: 预处理（流路/几何/边界）
├── Main.m                       ← 入口2: 稳态求解
├── PostProcessing.m             ← 入口3: 后处理（自动触发）
├── readme.md
├── .gitignore
├── PreProc/                     ← 预处理模块
│   ├── HX_Path_Planner1.m       GUI 流路设计器
│   ├── GenerateGeo.m            几何参数设置
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
  ② 环路优先压降扫描                 ② 广度优先压力更新
  ③ fsolve 流量重分配                ③ 残差收敛 → 退出
  ④ 广度优先压力更新
  ⑤ fsolve 二次流量重分配
  ⑥ 环路压降收敛 → 退出
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
  │     └─ alg(x0, BD, ..., CV, N, solver_flag, Prop_handle)  [嵌套子函数]
  │           ├─ R_cal_10()  → Prop1() ×3  工质侧换热+压降
  │           └─ DryA_cal_10()             空气侧换热+压降
  │
  ├─ scanTubes(pdropPaths, ...)       环路压降扫描
  ├─ fsolve(@uF)                      环路流量重分配
  ├─ scanTubes(heatPaths, ...)        压力场二次更新
  ├─ fsolve(@uF)                      二次流量重分配
  │
  └─ PostProcessing                   自动后处理
```

**嵌套子函数：**

| 函数 | 签名 | 职责 |
|------|------|------|
| `scanTubes` | `function [h_R_in, h_R_out, T_MA_in, T_MA_out, p_R_in, p_R_out, p_MA_in, p_MA_out, dp_tube, residual_max] = scanTubes(paths, predecessors, ...)` | 沿指定管序遍历，逐管逐 CV 调用 `alg()` |
| `alg` | `function F = alg(x0, BD, InletBD, GeoCondition, CV, N, solver_flag, Prop_handle)` | 单 CV 不动点迭代：solver_flag=1 换热量，=2 压降 |
| `uF` | `function F = uF(u, R_flow, N, mdot0)` | 环路基向量 F(u) = (dp_tube)' × N，供 fsolve 求解流量重分配 |

**依赖路径：** `Lib/` `Solver/` `PostProcessing/` `PreProc/`

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

| 结构体 | 产生于 | 关键字段 | 消耗于 |
|--------|--------|----------|--------|
| `TCinf` | `HX_Path_Planner1` → `cbExport` | `TC_matrix` `inlet_num` `outlet_num` `FlowDirection` `row` `col` `con_num` `Tube_num` | `GenerateGeo` `mdot_Initial` `buildPath` `plotAlongPath` `summaryTable` `exportHxPerf` |
| `GeoCondition` | `GenerateGeo(TCinf)` | `L` `D_inner` `D_outer` `r` `A_R` `A_MA` `P_row` `P_col` `Fin_pitch` `dx_fin` `Tube_num` `row` `col` | `GenerateBD` `R_cal_10` `DryA_cal_10` `summaryTable` `exportHxPerf` |
| `BDCondition` | `GenerateBD(GeoCondition, Prop_handle)` | `BD_R( h_R_inlet mdot_R_inlet p_R_inlet )` `BD_MA( T_MA_inlet mdot_MA_inlet p_MA_inlet x_MA_inlet )` `CV_num` | `mdot_Initial` `R_cal_10` `DryA_cal_10` `summaryTable` `exportHxPerf` |
| `Prop_handle` | `Prop_load(libLoc, R, ...)` | ~26 个 `griddedInterpolant` 函数句柄（工质+湿空气） | `Prop1` `GenerateBD` `DryA_cal_10` `R_cal_10` 以及各后处理函数 |
| `N` | `mdot_Initial(...)` | Tube_num × n_loops 矩阵，零空间基向量 | `Main` `buildPath` `plotLoopBalance` `summaryTable` `exportHxPerf` |
| `hxPerf` | `exportHxPerf(...)` | `Geo` `Topo` `Boundary` `Flow` `Thermal` `Performance` | `fmincon` / `ga` 等优化器 |

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
- 工质侧摩擦系数: Haaland
- 工质侧单相换热: Gnielinski
- 工质侧两相换热: Cavallini and Zecchin
- 空气侧: Wang-Chi-Chang Plate-Fin (j 因子 + f 因子)

---

## 四、求解收敛判据

| 场景 | 判据 |
|------|------|
| 无环路 | `residual_max < 1e-3`（能量残差） |
| 有环路 | `max(abs((dp_tube')*N)) < 1e-6` **且** `residual_max < 1e-3` |
| 外层迭代上限 | 10 轮 |
| CV 内层不动点迭代上限 | 100 轮 |

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

1. **`buildPath.m` 函数名不一致：** 文件名为 `buildPath.m`，`Main.m` 调用为 `buildPath(...)`，但内部声明为 `function [...] = buildPaths(A, N)`（带 s）。Windows 不区分大小写不报错，跨平台需注意。

2. **`Main.m` / `PreProcessing.m` / `PostProcessing.m` 均为脚本：** 不是函数，直接在 base workspace 操作变量，无显式参数传递。`Main.m` 内含 3 个嵌套子函数（`uF` `alg` `scanTubes`），仅在 Main.m 作用域内可见。

3. **REFPROP 外部依赖：** 整个物性体系依赖 `getFluidProperty()`，需 REFPROP 安装并位于 MATLAB 路径中。

4. **压力单位约定：** 所有 `p_*` 变量为 MPa，`dp_tube` 为 MPa（显示时 ×1e6 转 Pa）。

5. **HX_Path_Planner1 窗口持久化：** 点击关闭按钮隐藏而非销毁，句柄保存在 `hxDesigner`，`figure(hxDesigner)` 可随时唤出。
