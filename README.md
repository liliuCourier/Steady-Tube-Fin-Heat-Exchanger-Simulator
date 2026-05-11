# Steady-Tube-Fin-Heat-Exchanger-Simulator

**Steady Tube-Fin Heat Exchanger Simulator**

## 中文版

### 5/11 Demo1.0 发布

第一个适用于管翅式稳态仿真的 MATLAB 程序（工程）。在 Demo1.0 版本中只支持干工况。  
架构可以自己看看，内容不多。

### 仿真流程

1. 在仿真时先打开 `HX_Path_Planner1.m` 进行简单的流路设计。**不建议**在导出流路信息 `TCinf` 后就关闭界面，保持界面以实现实时修改流路。
2. 在 `GenerateGeo.m` 中进行具体的管道、翅片设计。设计完后不用运行，`Main.m` 求解程序会运行获取信息。
3. 在 `Main.m` 主程序中开始仿真。
4. 在 `PostProcessing.m` 中使用 MATLAB 进行后处理。

### 求解器

1. 首先使用 `mdot_Initial.m` 分析流路，找到所有环路。
2. 基于流路和环路，`buildPath.m` 生成广度优先路径和环路优先路径。
3. 内部求解采用换热和压降解耦，一次循环求解的顺序如下：
   - a. 固定压力场，沿广度优先路径更新热力场；
   - b. 固定热力场，沿环路优先路径计算所有环路压降，判断环路压降是否收敛，不收敛时根据环路压阻更新流量；
   - c. 仍然固定热力场，在刚才更新过的流量下沿广度优先路径更新压力场；
   - d. 再一次判断，收敛即退出，否则进入下次循环。
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

### Simulation Workflow

1. Run `HX_Path_Planner1.m` first to perform simple circuit design. **It is not recommended** to close the interface after exporting the circuit info `TCinf` – keep it open to enable real-time circuit modification.
2. In `GenerateGeo`, design the specific tube and fin geometry. No need to run it after design; the `Main` solver will obtain the information.
3. Start the simulation in the `Main` program.
4. Use `PostProcessing` in MATLAB for post-processing.

### Solver

1. First, use `mdot_Initial.m` to analyze the circuits and find all loops.
2. Based on circuits and loops, `buildPath.m` generates breadth‑first paths and loop‑first paths.
3. The internal solver decouples heat transfer and pressure drop. The sequence of one iteration is:
   - a. Fix the pressure field, update the thermal field along breadth‑first paths.
   - b. Fix the thermal field, calculate pressure drop for all loops along loop‑first paths. Check convergence; if not converged, update mass flow rates according to loop resistance.
   - c. Still with the thermal field fixed, update the pressure field along breadth‑first paths using the newly updated flow rates.
   - d. Check again; if converged, exit; otherwise proceed to the next iteration.
4. Both heat transfer and pressure drop solvers use fixed‑point iteration. Flow rate update uses `fsolve` to solve the resistance equations.
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

Different correlations lead to different results. Please set the desired correlations in the "correlation section" of `R_cal_10.m` (refrigerant‑side solver) and `DryA_cal_10.m` (air‑side solver).

**Default correlations:**

- **Refrigerant side:**
  - Pressure drop: Haaland equation for Darcy friction factor (reference `Simscape_Pipe(2P)`)
  - Heat transfer:
    - Single‑phase: Gnielinski correlation
    - Two‑phase: Cavallini and Zecchin correlation (reference `Simscape_Pipe(2P)`)

- **Air side:**
  - Pressure drop: WangChiChangPlateFin
  - Heat transfer: WangChiChangPlateFin

---

