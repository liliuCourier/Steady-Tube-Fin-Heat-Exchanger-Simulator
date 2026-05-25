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

### 5/22 Demo1.11

求解器算法更新：压阻模型指数可变化 + 换热/压力统一路径。

**算法变更：**
- 引入 `R_coef = 1.81`（基于 Blasius 摩擦因子关系 $\Delta p \propto m^{1.75}$）替代硬编码平方指数
- `uF` 函数新增 `R_coef` 参数，`fsolve` 压阻方程求解同步参数化
- 移除首次迭代中单独的环路优先路径（`pdropPaths`）压力扫描，换热与压力统一沿 `heatPaths` 广度优先路径计算
- 新增 `tube_cal` 计数器记录总扫描操作次数
- 收敛判据（环路压降 + 残差）统一移至外层迭代末尾
- 添加压阻指数诊断注释 `log(dp_tube*1e6./R_flow)./log(mdot_R)` 供标定参考

**Bug 修复：**
- 环路压降收敛历史图 `dp_loop_history` 缺少 `*1e6` 单位转换（内部 MPa → 显示 Pa）

**为何改指数：** Domanski (1989) 最早采用 $\Delta p = R \cdot m^{1.75}$（湍流 Blasius 标度），Ding (2004) 沿用。$1.81$ 比 $2.0$ 更贴近管内流动的物理标度，且配合 `heatPaths` 统一路径，避免了 `pdropPaths` 额外扫描带来的计算开销。

### 性能优化历史

| 日期 | 提交 | 优化内容 | 物性调用削减 | 分支 |
|------|------|----------|-------------|------|
| 5/25 | `5f39210` | **入口物性缓存**：`alg` 不动点循环前预计算入口 `Prop1`，循环内复用。入口 `(p,h)` 在迭代中不变，原先每轮重复查询 | ~50%（每 CV 迭代从 2 次→1 次） | `perf-property-cache` |
| 5/25 | `3a94197` | **物性算术平均**：取消平均态 `(p_CV,h_CV)` 处的第三次 `Prop1` 调用，14 个物性取进出口算术平均 | ~33%（每 CV 从 3 次→2 次） | `auto-circuit` |

**当前状态**：两项合计，`R_cal_10` 每 CV 每迭代从原先 3 次 Prop1（进口/出口/平均态）降为 1 次（仅出口），REFROP 调用总次数降为原先的 1/3。已知子函数 0.417s 耗时中 0.35s 在物性调用，预计子函数耗时降至 ~0.18s。

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

#### 附加功能依赖：MATLAB Interface for REFPROP and CoolProp

本程序通过 **MATLAB Interface for REFPROP and CoolProp**（官方 MATLAB 接口）调用 REFPROP，该附加功能**必须**由用户先获取并安装。

**系统要求：**

| 依赖 | 版本要求 | 说明 |
|------|----------|------|
| MATLAB | R2020a 或更新（推荐最新版） | 可能兼容更早版本，但不保证 |
| REFPROP | 10.x | `Refprop.dll` 所在目录需在 MATLAB 路径中 |
| CoolProp | 6.6.0 | 通常安装在 `C:\Users\<用户名>\AppData\Roaming\CoolProp` |
| C/C++ 编译器 | 当前 MATLAB 版本支持的编译器 | [编译器支持列表](https://www.mathworks.com/support/requirements/supported-compilers.html) |

**安装步骤：**
1. 安装 REFPROP 10.x 并确认 `Refprop.dll` 路径
2. 安装 CoolProp 6.6.0（通常自动安装至 `%APPDATA%\CoolProp`）
3. 配置 MATLAB 的 C/C++ 编译器：在 MATLAB 中运行 `mex -setup`
4. 将 MATLAB Interface for REFPROP and CoolProp 添加至 MATLAB 路径

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

## Git 提交历史 / Git Commit History

| SHA | 日期 | 说明 |
|-----|------|------|
| `5f39210` | 5/25 | perf: 入口物性缓存 — alg 循环内消除冗余 Prop1 调用 |
| `3a94197` | 5/25 | R_cal_10: 物性平均策略优化 — 取消第三次 Prop1 调用 |
| `7cddae7` | 5/22 | PostProcessing: 修复 SceneNode 警告 — LaTeX 改 TeX |
| `e518286` | 5/22 | Main_loop_base: 两阶段环路优先求解器 (Phase1+Phase2) |
| `9f9cb1e` | 5/22 | 新增 Main_loop_base: 环路优先求解器 |
| `0fdac6e` | 5/22 | 收敛轨迹可视化: 管扫描级 + 迭代级 + 综合指标 |
| `6721b44` | 5/22 | readme: Demo1.11 补充 PostProcessing 单位转换修复记录 |
| `e2f276b` | 5/22 | PostProcessing: 修复 dp_loop_history 缺少 *1e6 单位转换 |
| `586876e` | 5/22 | readme: Demo1.11 中英文版本记录 |
| `903820f` | 5/22 | Main.m: 压阻指数可变化 + 换热/压力统一路径 + tube_cal |
| `ab13412` | 5/12 | 修复 cbExport: 流路尺寸变化时跳过 diff_mat 相减 |
| `38cba0e` | 5/12 | 修复 HX_Path_Planner1 cbResetConn 中 tcbIn 笔误 |
| `f8a344a` | 5/12 | Main.m: 初始化 R_flow=[]，修复无环路时 PostProcessing 报错 |
| `65fe423` | 5/12 | readme: 补充 REFPROP/CoolProp 附加功能依赖说明 |
| `4fdaa4c` | 5/12 | 归档至 Demo1.1: 合并 Demo1.01/1.02，完善 REFPROP 配置指南 |
| `ad648bb` | 5/12 | 添加 Codelogic.md: 完整代码逻辑与数据结构文档 |
| `7c0e3ec` | 5/12 | 重构目录结构：函数分模块存放，根目录仅保留3个入口脚本 |
| `85f733e` | 5/12 | CV_num 移入 GenerateBD 对话框 |
| `5c4594d` | 5/12 | Main.m: 补充 dp_loop_max 变量定义 |
| `2b88e1f` | 5/12 | 删除 plotTubeLayout: 流路拓扑已由 HX_Path_Planner1 覆盖 |
| `c6c3344` | 5/12 | 重构工作流：预处理分离 + GenerateGeo UI + 自动后处理 |
| `688e869` | 5/12 | GenerateBD: 新增边界条件设置对话框 |
| `1f4aebd` | 5/12 | readme: Demo1.02 仿真流程 |
| `22f8db4` | 5/12 | HX_Path_Planner1: 关闭按钮改为隐藏窗口 |
| `232277a` | 5/12 | HX_Path_Planner1: 流路设计窗口持久化与增量导出 |
| `f23b355` | 5/12 | HX_Path_Planner1: 左侧进风方向指示 |
| `7bc367f` | 5/12 | 修复 plotAlongPath 累积压降曲线：有环路时沿前驱路径回溯 |
| `60dc175` | 5/12 | 修复 summaryTable/exportHxPerf 中有环路时压降计算错误 |
| `bfe6e54` | 5/12 | readme: Demo1.02 新增所有生成图像说明 |
| `d4ce3f4` | 5/12 | 简化 plotLoopBalance 子图2：环路不平衡量柱状图 |
| `9600ba7` | 5/12 | 修复 plotLoopBalance: range(ylim) 替换为 diff(ylim) |
| `84c6fad` | 5/12 | 修复后处理中压力单位错误：p_*变量为MPa |
| `96917a0` | 5/12 | readme.md 更新至 Demo1.02 |
| `0dd8013` | 5/12 | 后处理函数统一移至 PostProcessing/ 子目录 |
| `b8ddb72` | 5/12 | Demo1.02: 完成后处理六层模块 |
| `1e518f8` | 5/11 | 添加 readme.md |
| `38c641b` | 5/11 | Demo1.01 初始版本 |

> 分支: `auto-circuit` | 共 36 个提交 | 截至 2026-05-25 未推送 GitHub

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

### 5/22 Demo1.11 Release

Solver algorithm update: variable-exponent pressure-resistance model + unified heat/pressure path.

**Algorithm changes:**
- Introduced `R_coef = 1.81` (based on Blasius friction factor scaling $\Delta p \propto m^{1.75}$) replacing hardcoded square exponent
- `uF` function gains `R_coef` parameter; `fsolve` resistance-equation solving is now parameterized
- Removed separate first-iteration loop-priority path (`pdropPaths`) pressure scan; heat and pressure now both use the same breadth-first path (`heatPaths`)
- Added `tube_cal` counter to track total scan operations
- Convergence criterion (loop pressure drop + residual) unified at outer-loop end
- Added diagnostic comment `log(dp_tube*1e6./R_flow)./log(mdot_R)` for exponent calibration

**Bug fix:**
- Loop pressure-drop convergence history plot: `dp_loop_history` missing `*1e6` unit conversion (internal MPa → display Pa)

**Why change the exponent:** Domanski (1989) first adopted $\Delta p = R \cdot m^{1.75}$ (turbulent Blasius scaling), followed by Ding (2004). $1.81$ is closer to the physical scaling of in-tube flow than $2.0$, and the unified heatPaths scan avoids the extra computational overhead of the separate pdropPaths pass.

### Performance Optimization History

| Date | Commit | Optimization | Prop1 Calls Reduced | Branch |
|------|--------|-------------|-------------------|--------|
| 5/25 | `5f39210` | **Inlet property cache**: Pre-compute inlet `Prop1` before `alg` fixed-point loop; inlet `(p,h)` unchanged during iteration | ~50% (2→1 per CV iteration) | `perf-property-cache` |
| 5/25 | `3a94197` | **Arithmetic averaging**: Eliminated 3rd `Prop1` call at average `(p_CV,h_CV)` state; 14 properties now arithmetic means of inlet/outlet | ~33% (3→2 per CV) | `auto-circuit` |

**Current status**: Combined, `R_cal_10` per-CV per-iteration Prop1 calls reduced from 3 (inlet/outlet/average) to 1 (outlet only). Total REFPROP queries at ~1/3 of original. With 0.35s of 0.417s sub-function time spent on property calls, estimated sub-function time drops to ~0.18s.

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

#### External dependency: MATLAB Interface for REFPROP and CoolProp

This program calls REFPROP through the **MATLAB Interface for REFPROP and CoolProp** (official MATLAB wrapper). This add-on **must** be obtained and installed by the user before running.

**System requirements:**

| Dependency | Version | Notes |
|------------|---------|-------|
| MATLAB | R2020a or later (latest recommended) | May work on older versions, not guaranteed |
| REFPROP | 10.x | Directory containing `Refprop.dll` must be on MATLAB path |
| CoolProp | 6.6.0 | Typically installed at `C:\Users\<user>\AppData\Roaming\CoolProp` |
| C/C++ Compiler | Supported by your MATLAB release | [Compiler support list](https://www.mathworks.com/support/requirements/supported-compilers.html) |

**Installation steps:**
1. Install REFPROP 10.x and note the `Refprop.dll` path
2. Install CoolProp 6.6.0 (usually auto-installed to `%APPDATA%\CoolProp`)
3. Configure MATLAB C/C++ compiler: run `mex -setup` in MATLAB
4. Add the MATLAB Interface for REFPROP and CoolProp to the MATLAB path

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
