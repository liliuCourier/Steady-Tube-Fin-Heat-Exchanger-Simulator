# Steady-Tube-Fin-Heat-Exchanger-Simulator

**Steady Tube-Fin Heat Exchanger Simulator**

## 中文版

### 5/11 Demo1.0

第一个适用于管翅式稳态仿真的 MATLAB 程序（工程）。在 Demo1.0 版本中只支持干工况。架构可以自己看，内容不多。

### 5/12 Demo1.1

基于 Demo1.0 全面重构：修复无环路求解 bug、完成后处理六层模块、重构预处理工作流、工程目录模块化。

**架构变更：**
- 根目录仅保留 3 个入口脚本：`PreProcessing.m` `Main.m` `PostProcessing.m`
- 函数按模块分入 `PreProc/` `Solver/` `Lib/` `PostProcessing/` 四个子目录
- 求解与前处理完全分离：先运行 `PreProcessing` 完成三类设置，再运行 `Main` 求解

**Bug 修复：**
- 无环路时 `N` 为空导致 `fsolve` 报错 → 检测 `N` 是否为空，有/无环路分离两套迭代逻辑

**预处理模块（PreProcessing.m + PreProc/）：**
- `HX_Path_Planner1`：流路设计窗口持久化（关闭=隐藏），`hxDesigner` 句柄常驻，增量导出对比
- `GenerateGeo`：交互式几何参数对话框（管长/管径/翅片间距等 8 项）
- `GenerateBD`：交互式边界条件对话框（工质/空气进口参数 + 控制容积数 CV_num）

**求解模块（Main.m + Solver/）：**
- 收敛历史自动记录，求解完成后自动触发 `PostProcessing`
- 有环路：环路压降扫描 → fsolve 流量重分配 → 压力更新 → 二次重分配 → 收敛
- 无环路：广度优先扫描 → 直接压力更新 → 残差收敛

**后处理模块（PostProcessing 自动触发）：**
1. `plotAlongPath.m` — 2×2 沿线物性分布曲线（温度/压力/换热量/累积压降）
2. `plotLoopBalance.m` — 环路压降平衡双图（支路对比 + 不平衡量），无环路自动跳过
3. `summaryTable.m` — 控制台性能汇总表（换热量/压降/温差/各管占比/出口状态矩阵）
4. 收敛历史双图 — 能量残差 + 环路压降残差（对数坐标）
5. `exportHxPerf.m` — 标准化 `hxPerf` 结构体，供优化器直接调用

### 仿真流程

1. 运行 `PreProcessing` 完成三类预处理：
   - **流路设计**：打开 `HX_Path_Planner1` 设计窗口，设置进出口、绘制 U 型弯连线，点击"导出结果"
   - **几何设置**：弹出对话框设置管长、管径、翅片间距等参数
   - **边界条件**：弹出对话框设置工质/空气进口参数
2. 运行 `Main` 开始求解，完成后自动执行 `PostProcessing` 后处理。
3. 如需修改流路或参数，在 `hxDesigner` 窗口中调整后重新导出 TCinf，重新运行 `Main`。

### 求解器

1. 首先使用 `mdot_Initial.m` 分析流路，找到所有环路。若无环路则返回空基向量。
2. 基于流路和环路，`buildPath.m` 生成广度优先路径和环路优先路径。
3. 内部求解采用换热和压降解耦，一次循环求解的顺序如下：
   - **有环路时**：
     - a. 固定压力场，沿广度优先路径更新热力场；
     - b. 固定热力场，沿环路优先路径计算所有环路压降，判断环路压降是否收敛，不收敛时根据环路压阻更新流量；
     - c. 仍然固定热力场，在刚才更新过的流量下沿广度优先路径更新压力场；
     - d. 再一次判断，收敛即退出，否则进入下次循环。
   - **无环路时**：
     - a. 固定压力场，沿广度优先路径更新热力场；
     - b. 直接沿广度优先路径更新压力场，残差收敛即退出。
4. 换热求解和压降求解均采用不动点迭代法，流量更新采用 `fsolve` 求解压阻方程。
5. 在 16 根管、2 进 2 出（环路数为 4）的计算中，流量分配与 CoilDesigner 计算结果基本一致。

### 注意事项

#### ⚠️ 1. 物性数据获取（运行前必读）

本程序使用 **REFPROP** 获取工质与湿空气物性。`Prop_load` 在 `PreProcessing.m` 和 `Main.m` 中均会调用（如未提前加载）。

**必须修改以下配置以匹配你的环境：**

```matlab
% 在 PreProcessing.m 和 Main.m 开头修改：
refprop_location = 'E:\refprop10\REFPROP';   % ← 改为你的 REFPROP 安装路径（Refprop.dll 所在目录）
R = 'R134a';                                   % ← 改为你需要的工质（REFPROP 支持的命名）
% Prop_load 调用格式:
% Prop_load(libLoc, R, pmin, pmax, hmin, hmax, p_point, u_vap_point, u_liq_point)
Prop_handle = Prop_load(refprop_location, R, 1e-3, 5.5, 80, 510, 100, 25, 25);
%                                            ↑     ↑    ↑   ↑    ↑    ↑
%                                          p_min  p_max h_min h_max p_pt h_pt
```

**参数说明：**

| 参数 | 含义 | 当前默认值 | 注意 |
|------|------|-----------|------|
| `pmin / pmax` | 有效压力范围 (MPa) | 1e-3 ~ 5.5 | 必须覆盖工况压力范围 |
| `hmin / hmax` | 有效焓范围 (kJ/kg) | 80 ~ 510 | 必须覆盖工况焓范围 |
| `p_point` | 压力离散点数 | 100 | 越大插值越精确，加载越慢 |
| `u_vap_point` | 气相归一化焓离散点 | 25 | 同上 |
| `u_liq_point` | 液相归一化焓离散点 | 25 | 同上 |

**常见问题：**
- 如果仿真过程中某管压力/焓超出上述范围，`Prop1` 会因插值外推而报错或返回异常值
- 若出现 `getFluidProperty` 未定义，请检查 REFPROP 是否正确安装且 MATLAB 路径包含 REFPROP 目录
- 不同工质（如 R410A）的饱和压力和焓范围差异很大，更换工质时务必调整 `pmin/pmax/hmin/hmax`

#### 2. 关联式选择

不同的关联式会产生不同的结果。具体的关联式请在 `R_cal_10.m`（工质侧求解器）和 `DryA_cal_10.m`（空气侧求解器）中的"关联式区域"自己设置。

**默认关联式：**

- **工质侧：**
  - 压损：采用 Haaland 公式获取 Darcy friction factor（参考 `Simscape_Pipe(2P)`）
  - 换热：
    - 单相：采用 Gnielinski correlation
    - 两相：采用 Cavallini and Zecchin correlation（参考 `Simscape_Pipe(2P)`）

- **空气侧：**
  - 压损：WangChiChangPlateFin
  - 换热：WangChiChangPlateFin

---

## English Version

# Steady-Tube-Fin-Heat-Exchanger-Simulator

**Steady Tube-Fin Heat Exchanger Simulator**

### 5/11 Demo1.0 Release

First MATLAB program (project) for steady-state simulation of tube-fin heat exchangers. Demo1.0 only supports dry-surface conditions.  
You can check the architecture yourself – it's not complicated.

### 5/12 Demo1.1 Release

Comprehensive refactoring from Demo1.0: fixed no-loop solver bug, completed six-layer post-processing suite, restructured preprocessing workflow, and modularized the project directory.

**Architecture changes:**
- Root directory contains only 3 entry scripts: `PreProcessing.m` `Main.m` `PostProcessing.m`
- Functions organized into 4 subdirectories: `PreProc/` `Solver/` `Lib/` `PostProcessing/`
- Solving fully decoupled from preprocessing: run `PreProcessing` first for three setup steps, then `Main` for solving

**Bug fixes:**
- Empty `N` null-space causing `fsolve` failure → two-branch iteration logic for with/without loops

**Preprocessing (PreProcessing.m + PreProc/):**
- `HX_Path_Planner1`: persistent design window (close = hide), `hxDesigner` handle, incremental export diff
- `GenerateGeo`: interactive geometry parameter dialog (8 parameters including tube length/diameters/fin pitch)
- `GenerateBD`: interactive boundary condition dialog (refrigerant/air inlet + CV_num)

**Solver (Main.m + Solver/):**
- Convergence history auto-logged; `PostProcessing` auto-triggered after solving
- With loops: loop pressure scan → fsolve flow redistribution → pressure update → second redistribution → converge
- Without loops: BFS scan → direct pressure update → residual convergence

**Post-processing (auto-triggered from Main.m):**
1. `plotAlongPath.m` — 2×2 parameter distribution curves (temperature/pressure/heat load/cumulative dp)
2. `plotLoopBalance.m` — Loop pressure balance dual-chart (branch comparison + imbalance), auto-skips if no loops
3. `summaryTable.m` — Console performance summary (heat load/dp/temperatures/per-tube breakdown/outlet state matrix)
4. Convergence history dual-plot — energy residual + loop pressure residual (log scale)
5. `exportHxPerf.m` — Standardized `hxPerf` struct for optimizer integration

### Simulation Workflow

1. Run `PreProcessing` to complete three setup steps:
   - **Circuit design**: Opens `HX_Path_Planner1` design window. Set inlets/outlets, draw U-bend connections, click "Export".
   - **Geometry setup**: Dialog for tube length, diameters, fin pitch, and other geometric parameters.
   - **Boundary conditions**: Dialog for refrigerant/air inlet parameters.
2. Run `Main` to start simulation. Post-processing (`PostProcessing`) runs automatically upon completion.
3. To iterate, modify the circuit in the `hxDesigner` window, re-export, and re-run `Main`.

### Solver

1. First, use `mdot_Initial.m` to analyze the circuits and find all loops. If no loops exist, return an empty null-space basis.
2. Based on circuits and loops, `buildPath.m` generates breadth-first paths and loop-first paths.
3. The internal solver decouples heat transfer and pressure drop. The sequence of one iteration is:
   - **With loops**:
     - a. Fix the pressure field, update the thermal field along breadth-first paths.
     - b. Fix the thermal field, calculate pressure drop for all loops along loop-first paths. Check convergence; if not converged, update mass flow rates according to loop resistance.
     - c. Still with the thermal field fixed, update the pressure field along breadth-first paths using the newly updated flow rates.
     - d. Check again; if converged, exit; otherwise proceed to the next iteration.
   - **Without loops**:
     - a. Fix the pressure field, update the thermal field along breadth-first paths.
     - b. Directly update the pressure field along breadth-first paths; exit when residual converges.
4. Both heat transfer and pressure drop solvers use fixed-point iteration. Flow rate update uses `fsolve` to solve the resistance equations.
5. For a case with 16 tubes, 2 inlets and 2 outlets (4 loops), the flow distribution is essentially consistent with CoilDesigner results.

### Important Notes

#### ⚠️ 1. Property data (must read before running)

REFPROP is used for refrigerant and moist-air property data. `Prop_load` is called in both `PreProcessing.m` and `Main.m` (if not already loaded).

**You must modify these settings to match your environment:**

```matlab
% In PreProcessing.m and Main.m:
refprop_location = 'E:\refprop10\REFPROP';   % ← Change to your REFPROP path
R = 'R134a';                                   % ← Change to your desired fluid
% Prop_load signature:
% Prop_load(libLoc, R, pmin, pmax, hmin, hmax, p_point, u_vap_point, u_liq_point)
Prop_handle = Prop_load(refprop_location, R, 1e-3, 5.5, 80, 510, 100, 25, 25);
```

**Parameter guide:**

| Parameter | Meaning | Default | Note |
|-----------|---------|---------|------|
| `pmin / pmax` | Valid pressure range (MPa) | 1e-3 ~ 5.5 | Must span your operating pressures |
| `hmin / hmax` | Valid enthalpy range (kJ/kg) | 80 ~ 510 | Must span your operating enthalpies |
| `p_point` | Pressure discretization points | 100 | More = finer interpolation, slower load |
| `u_vap_point` | Vapor normalized-enthalpy points | 25 | Same tradeoff |
| `u_liq_point` | Liquid normalized-enthalpy points | 25 | Same tradeoff |

**Common issues:**
- If a tube's pressure/enthalpy exceeds the above ranges during simulation, `Prop1` will error or return invalid values due to extrapolation
- If `getFluidProperty` is undefined, check that REFPROP is installed and its directory is on the MATLAB path
- Different fluids (e.g., R410A) have very different saturation pressures and enthalpy ranges — adjust `pmin/pmax/hmin/hmax` accordingly

#### 2. Correlation selection

Different correlations lead to different results. Please set the desired correlations in the "correlation section" of `R_cal_10.m` (refrigerant-side solver) and `DryA_cal_10.m` (air-side solver).

**Default correlations:**

- **Refrigerant side:**
  - Pressure drop: Haaland equation for Darcy friction factor (reference `Simscape_Pipe(2P)`)
  - Heat transfer:
    - Single-phase: Gnielinski correlation
    - Two-phase: Cavallini and Zecchin correlation (reference `Simscape_Pipe(2P)`)

- **Air side:**
  - Pressure drop: WangChiChangPlateFin
  - Heat transfer: WangChiChangPlateFin

---
