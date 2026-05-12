# Steady-Tube-Fin-Heat-Exchanger-Simulator

**Steady Tube-Fin Heat Exchanger Simulator**

## 中文版

### 5/11 Demo1.0

第一个适用于管翅式稳态仿真的 MATLAB 程序（工程）。在 Demo1.0 版本中只支持干工况。架构可以自己看，内容不多。

### 5/12 Demo1.01

基于 Demo1.0，修复了一个重大 bug：在无环路简单流路中，基向量 $N$ 为空，所有管路流量由进口唯一确定，不存在自由变量。此时 `mdot_Initial.m` 中 $u_0 = N \setminus (m_{init} - m_0)$ 会因 $N$ 为空而报错，且 `Main.m` 中 `fsolve` 收到空初值同样报错。

**修复内容：**

1. `mdot_Initial.m`：检测 $N$ 是否为空，为空时直接设 $u_0 = []$。
2. `Main.m`：将求解循环分为有环路 / 无环路两套逻辑：
   - **无环路**：广度优先扫描 → 直接第二次广度优先压力场更新 → 残差收敛即退出，完全跳过环路压降扫描和流量重分配。
   - **有环路**：保持原有逻辑不变（环路压降扫描 → 流量重分配 → 压力场更新 → 二次流量重分配 → 环路压降收敛判断）。

### 5/12 Demo1.02

完成后处理六层模块，提供从流路验证到性能汇总的完整可视化与数据导出能力。所有后处理函数统一放置在 `PostProcessing/` 子目录中，`PostProcessing.m` 作为一键集成入口。

**新增模块与生成图像：**

1. **`PostProcessing/plotTubeLayout.m`** — 流路拓扑可视化
   - 按实际管排布局（$N_{row} \times N_{col}$）绘制管束网格
   - 蓝色箭头 = 正向流动，橙色箭头 = 反向流动
   - 绿色方框标记进口管，红色三角标记出口管
   - 黑线连接同一节点的管路，直观展示分支/汇合关系

2. **`PostProcessing/plotAlongPath.m`** — 沿线物性分布（2×2 子图面板）
   - 子图1：工质出口温度 + 空气平均温度沿 BFS 管序变化
   - 子图2：工质平均压力（左轴）+ 各管压降（右轴柱状图）
   - 子图3：各管换热量柱状图，标注占比百分比
   - 子图4：沿程累积压降曲线

3. **`PostProcessing/plotLoopBalance.m`** — 环路性能对比（双图）
   - 左图：每个环路两支路压降分组柱状图（正向 vs 反向），标注不平衡量
   - 右图：各环路压降不平衡量柱状图（越矮越平衡）
   - 无环路时自动跳过

4. **`PostProcessing/summaryTable.m`** — 控制台整体性能汇总表
   - 基本参数：管排布局、进出口管号、计算耗时
   - 热力性能：总换热量、进出口温度、最小传热温差
   - 流动性能：总流量、总压降、各管压降之和、各管流量
   - 各管换热量占比（含 ASCII 进度条）
   - 各管出口状态矩阵（焓/压力/温度/压降/换热量）

5. **`PostProcessing/exportHxPerf.m`** — 标准化 `hxPerf` 结构体导出，供 `fmincon`、遗传算法等优化器直接调用

6. **`Main.m`** — 新增迭代收敛历史日志（`residual_history`、`dp_loop_history`），`PostProcessing.m` 绘制双图：
   - 左图：能量残差随迭代次数的收敛曲线（对数坐标）
   - 右图：环路压降残差收敛曲线（无环路时显示 N/A）

7. **`PostProcessing.m`** — 一键集成入口，自动调用上述全部模块

8. **`HX_Path_Planner1.m`** — 流路设计窗口持久化与增量导出
   - Figure 添加唯一 `Tag` + 句柄持久化到 `hxDesigner`，随时 `figure(hxDesigner)` 唤出
   - 点击关闭按钮隐藏窗口而非销毁，设计会话不丢失
   - 设计变更后标题栏显示 `[已修改，未导出]`，导出后自动清除
   - 导出时对比上次 `TCinf`，输出拓扑/进出口变化摘要
   - 左侧添加绿色箭头 + `AIR` 标注指示进风方向
   - `Main.m` 求解完成后自动弹回设计窗口

### 仿真流程

1. 运行 `HX_Path_Planner1` 打开设计窗口，设置进出口、绘制 U 型弯连线，点击"导出结果"。
2. 窗口可最小化或关闭（自动隐藏），句柄 `hxDesigner` 始终有效。
3. 在 `GenerateGeo.m` 中进行具体的管道、翅片设计。设计完后不用运行，`Main.m` 求解程序会运行获取信息。
4. 运行 `Main` 开始仿真，求解完成后自动弹回设计窗口。
5. 如需修改流路，直接在弹回的设计窗口中调整，再次导出（控制台显示变更摘要），重新运行 `Main`。
6. 运行 `PostProcessing` 进行后处理可视化与数据导出。

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

#### 1. 物性数据获取

在 `main` 求解主程序中，本程序使用 **Refprop** 获取物性数据，因此用户需要提供其 Refprop 安装位置（`Refprop.dll` 所在位置），并将 `R` 填写为 Refprop 支持的工质。

示例：
```matlab
refprop_location = 'E:\refprop10\REFPROP';
R = 'R134a';
% Prop_handle = Prop_load(refprop_location,R,pmin,pmax,hmin,hmax,p_point,u_vap_point,u_liq_point);
Prop_handle = Prop_load(refprop_location,R,1e-3,5.5,80,510,100,25,25);
```

调用的物性大致需要一个有效压力范围和有效焓范围（求解器基于压力和焓进行求解）。

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

### 5/12 Demo1.01 Release

Based on Demo1.0, this release fixes a critical bug: in simple circuits without loops, the null-space basis $N$ is empty, meaning all tube flow rates are uniquely determined by the inlet flow rate with no free variables. Previously `mdot_Initial.m` would fail computing $u_0 = N \setminus (m_{init} - m_0)$ with an empty $N$, and `Main.m` would fail when `fsolve` received an empty initial point.

**Fixes:**

1. `mdot_Initial.m`: Detect empty $N$ and set $u_0 = []$ directly.
2. `Main.m`: Split the iteration loop into two logic branches:
   - **No loops**: breadth-first scan → direct second breadth-first pressure update → exit on residual convergence, skipping loop pressure drop scan and flow redistribution entirely.
   - **With loops**: original logic preserved.

### 5/12 Demo1.02 Release

Completed a six-layer post-processing suite providing full visualization and data export — from circuit design verification to performance summaries. All post-processing functions reside in the `PostProcessing/` subdirectory, with `PostProcessing.m` as the one-click entry point.

**New modules and generated figures:**

1. **`PostProcessing/plotTubeLayout.m`** — Flow topology visualization
   - Draws tube grid in actual layout ($N_{row} \times N_{col}$)
   - Blue arrows = forward flow, orange arrows = reverse flow
   - Green square marks inlet tubes, red triangle marks outlet tubes
   - Black lines connect tubes sharing a node, showing branching/merging relationships

2. **`PostProcessing/plotAlongPath.m`** — Parameter distribution along BFS path (2×2 subplot panel)
   - Subplot 1: Refrigerant outlet temperature + air average temperature along path
   - Subplot 2: Refrigerant average pressure (left axis) + per-tube pressure drop (right axis bar)
   - Subplot 3: Per-tube heat load bar chart with percentage labels
   - Subplot 4: Cumulative pressure drop along path

3. **`PostProcessing/plotLoopBalance.m`** — Loop balance verification (dual-panel)
   - Left: Grouped bar chart of forward vs. reverse branch pressure drop for each loop, with imbalance annotations
   - Right: Loop pressure imbalance bar chart (shorter = more balanced)
   - Auto-skips when no loops exist

4. **`PostProcessing/summaryTable.m`** — Console overall performance summary table
   - Basic parameters: tube layout, inlet/outlet tube numbers, computation time
   - Thermal performance: total heat load, inlet/outlet temperatures, minimum approach temperature
   - Flow performance: total mass flow, total pressure drop, per-tube flow rates
   - Per-tube heat load percentages (with ASCII bar chart)
   - Per-tube outlet state matrix (enthalpy/pressure/temperature/pressure drop/heat load)

5. **`PostProcessing/exportHxPerf.m`** — Standardized `hxPerf` struct export for direct use by `fmincon`, genetic algorithms, and other optimizers

6. **`Main.m`** — Added iteration convergence history logging (`residual_history`, `dp_loop_history`), with dual-panel plot in `PostProcessing.m`:
   - Left: Energy residual vs. iteration (log scale)
   - Right: Loop pressure residual convergence (N/A when no loops)

7. **`PostProcessing.m`** — One-click entry point that calls all modules above

8. **`HX_Path_Planner1.m`** — Persistent design window + incremental export
   - Unique `Tag` + handle persisted as `hxDesigner` in workspace; call `figure(hxDesigner)` anytime
   - Close button hides the window instead of destroying it; design session survives
   - Title bar shows `[已修改，未导出]` when circuit is modified, cleared on export
   - Incremental export: compares with last `TCinf`, prints diff summary (topology, inlets, outlets)
   - Green arrow + `AIR` label on the left side indicates air inlet direction
   - `Main.m` automatically brings the design window back after solving

### Simulation Workflow

1. Run `HX_Path_Planner1` to open the design window. Set inlets/outlets, draw U-bend connections, click "Export".
2. The window can be minimized or "closed" (auto-hidden); handle `hxDesigner` stays valid.
3. Design tube and fin geometry in `GenerateGeo`. No need to run it; `Main` calls it automatically.
4. Run `Main` to start simulation. The design window pops back up automatically after solving.
5. To iterate, modify the circuit in the returned design window, re-export (console shows change summary), re-run `Main`.
6. Run `PostProcessing` for visualization and data export.

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

#### 1. Property data

In the main solver, Refprop is used to obtain fluid properties. The user must provide the Refprop installation path (the directory containing `Refprop.dll`) and specify the refrigerant `R` supported by Refprop.

Example:
```matlab
refprop_location = 'E:\refprop10\REFPROP';
R = 'R134a';
% Prop_handle = Prop_load(refprop_location,R,pmin,pmax,hmin,hmax,p_point,u_vap_point,u_liq_point);
Prop_handle = Prop_load(refprop_location,R,1e-3,5.5,80,510,100,25,25);
```

The required property calls need a valid pressure range and enthalpy range (the solver works with pressure and enthalpy).

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
