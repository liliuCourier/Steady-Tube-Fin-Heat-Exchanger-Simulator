# 工质单相压降关联式 — 文献总结

> 对应关联式：§8 TFRLDP/TFRVDP — 管翅式工质纯液/纯气压降
>
> 整理日期：2026-05-26
>
> 用途：换热器仿真器中工质侧单相摩擦压降计算

---

## 0. 约定：Fanning 摩擦因子 vs Darcy 摩擦因子

所有 Goto 和 Carnavos 关联式均使用 **Fanning 摩擦因子** $f$。Swamee-Jain 和 Churchill 原始论文使用 Darcy $f_D$。**代码统一使用 Fanning $f$**。

| | Fanning $f$ | Darcy $f_D$ |
|---|---|---|
| 层流 | $f = 16/Re$ | $f_D = 64/Re$ |
| 压降公式 | $\Delta p = \dfrac{2fL G^2}{\rho D}$ | $\Delta p = \dfrac{f_D L G^2}{2\rho D}$ |
| 换算 | — | $f = f_D / 4$ |

本文档中所有公式除非特别说明，均以 **Fanning f** 给出。

---

## 1. 层流与过渡区处理

### 1.1 层流：$Re \le 2000$

$$f = \frac{16}{Re}$$

### 1.2 过渡区：$2000 < Re < 4000$

两种方案：

**方案 A — 线性插值**（简单，但不连续）：
```matlab
f_lam = 16 / 2000;
f_tur = turbulent_f(4000);  % 由湍流关联式计算
x = (Re - 2000) / (4000 - 2000);
f = (1-x) * f_lam + x * f_tur;
```

**方案 B — Churchill 公式**（推荐，全范围连续）：
参见 §2。

---

## 2. Churchill (1977) — 全范围连续摩擦因子

> 关联式名称：ChurchillRLDP / ChurchillRVDP / Churchill (MCRLDP)

**文献**：Churchill SW. "Frictional equation spans all fluid flow regimes". *Chemical Engineering*, 84(24), 91-92, 1977. (行业杂志，无 DOI)

**核心公式**（Darcy $f_D$ 形式，Fanning 需除以 4）：

$$f_D = 8\left[\left(\frac{8}{Re}\right)^{12} + \frac{1}{(A+B)^{3/2}}\right]^{1/12}$$

其中：

$$A = \left[2.457\ln\left(\frac{1}{(7/Re)^{0.9} + 0.27\varepsilon/D}\right)\right]^{16}$$

$$B = \left(\frac{37530}{Re}\right)^{16}$$

**特点**：
- 单一公式覆盖层流、过渡区、湍流，无需分段
- 光滑管取 $\varepsilon/D = 0$

**代码实现**（Fanning f）：
```matlab
function f = churchill(Re, eps_D)
    % Churchill 全范围摩擦因子 (Fanning)
    A = (2.457 * log(1 / ((7/Re)^0.9 + 0.27*eps_D)))^16;
    B = (37530 / Re)^16;
    f_D = 8 * ((8/Re)^12 + 1/(A+B)^1.5)^(1/12);
    f = f_D / 4;  % 转换为 Fanning
end
```

---

## 3. Swamee-Jain (1976) — 湍流光管摩擦因子

**文献**：Swamee PK, Jain AK. "Explicit equations for pipe-flow problems". *Journal of the Hydraulics Division, ASCE*, 102(5), 657-664, 1976.

**核心公式**（Darcy $f_D$，转换为 Fanning $f = f_D/4$）：

$$f = \frac{0.25}{4\left[\log_{10}\left(\frac{\varepsilon}{3.7D} + \frac{5.74}{Re^{0.9}}\right)\right]^2}$$

**适用范围**：
- $5000 < Re < 10^8$
- $10^{-6} < \varepsilon/D < 10^{-2}$

### Haaland (1983) 公式 — Swamee-Jain 的替代方案

相比 Swamee-Jain，Haaland 在低 Re 区更稳定，精度略优（偏差 ~0.5% vs ~1.0%）：

$$\frac{1}{\sqrt{f_D}} = -1.8\log_{10}\left[\left(\frac{\varepsilon}{3.7D}\right)^{1.11} + \frac{6.9}{Re}\right]$$

**代码实现**（Fanning f）：
```matlab
% Swamee-Jain
f_D = 0.25 / (log10(eps_D/3.7 + 5.74/Re^0.9))^2;
f = f_D / 4;

% Haaland (替代方案)
f_D = 1 / (-1.8 * log10((eps_D/3.7)^1.11 + 6.9/Re))^2;
f = f_D / 4;
```

**推荐用途**：光管湍流默认选项配合 Churchill 公式。铜管粗糙度 $\varepsilon \approx 1.5\times10^{-6}$ m（近似光滑）。

---

## 4. Carnavos (1980) — 内翅管 Fanning 摩擦因子

**文献**：Carnavos TC. "Heat Transfer Performance of Internally Finned Tubes in Turbulent Flow". *Heat Transfer Engineering*, 1(4), 32-37, 1980.

**测试条件**：11 种内翅管（Noranda Forge-Fin 制造工艺），管内径 6.86–23.8 mm，工质为水（Pr ≈ 6.5）、空气、乙二醇水溶液，Re > 10,000。

**核心公式**（Fanning $f$，±10%）：

$$f = \frac{0.046}{Re^{0.2}} \cdot F^*$$

$$F^* = \left(\frac{A_{fa}}{A_{fn}}\right)^{0.5} (\sec\alpha)^{0.75}$$

其中：
| 符号 | 含义 |
|---|---|
| $A_{fa}$ | 实际流通面积（扣除翅片占据面积） |
| $A_{fn}$ | 标称流通面积（按翅尖直径计算的圆面积） |
| $\alpha$ | 翅片螺旋角（°） |

**几何参数范围**（测试的 11 根管）：
- $A_{fa}/A_{fn}$ = 0.80–0.97
- 翅数 = 10–38
- 螺旋角 = 0°(纵向) 和 6–30°(螺旋)

**注意事项**：
- 螺旋角 $\alpha$ 在 $0^\circ$~$30^\circ$ 范围内验证
- 超出此范围时 $\sec\alpha$ 增长过快，公式可能失效
- 空气数据来自 Carnavos [1]，水数据来自本工作，两者吻合良好

---

## 5. Schlager (1989) — 微翅管增强因子法

**文献**：Schlager LM, Pate MB, Bergles AE. "Heat transfer and pressure drop during evaporation and condensation of R22 in horizontal micro-fin tubes". *International Journal of Refrigeration*, 12, 6-20, 1989.

**测试条件**：
- 3 种微翅管，外径 9.52 mm，最大内径 8.92 mm
- 工质 R22，质量流速 150–500 kg/(m²·s)
- 单相测试雷诺数 18,000–34,000

**核心结果**：

| 指标 | 数值 |
|---|---|
| 单相 HTC 增强因子 | 1.3–1.6（随 Re 增加略有下降）|
| 压降惩罚因子 $\Delta P_M / \Delta P_S$ | ≈ 1.4 |

**⚠️ 注意事项**：
- 单相液相测试数据有限，增强因子 1.4 为近似值
- CoilDesigner 中对该关联式的使用方式是 Gnielinski × 1.9（HTC）配合同一因子计算压降
- 优先使用 Carnavos 或 Goto 的直接几何关联式可获更高精度

---

## 6. Goto (2001) — 螺旋/人字沟管内翅管摩擦因子

**文献**：Goto M, Inoue N, Ishiwatari N. "Condensation and evaporation heat transfer of R410A inside internally grooved horizontal tubes". *International Journal of Refrigeration*, 24, 628-638, 2001.

**测试管参数**：

| 参数 | 螺旋沟管 (C Tube) | 人字沟管 (W Tube) |
|---|---|---|
| 外径 | 8.01 mm | 8.00 mm |
| 平均内径 | 7.30 mm | 7.22 mm |
| 翅高 | 0.17 mm | 0.24 mm |
| 螺旋角 | 18° | 19° |
| 翅数 | 55 | 60 |
| 实际面积/标称面积 | 1.498 | 2.032 |

**单相 Fanning 摩擦因子**（±5%）：

**C 管（螺旋沟）**：
| Re 范围 | $f$ |
|---|---|
| $2000 \leq Re < 6500$ | $f = 1.47 \times 10^{-4} \cdot Re^{0.53}$ |
| $6500 \leq Re < 12700$ | $f = 0.046 \cdot Re^{-0.20}$ |
| $12700 \leq Re \leq 46500$ | $f = 1.23 \times 10^{-2} \cdot Re^{0.21}$ |

**W 管（人字沟）**：
| Re 范围 | $f$ |
|---|---|
| $Re < 3900$ | $f = 2.17 \times 10^{-2} \cdot Re^{-0.08}$ |
| $3900 \leq Re < 11500$ | $f = 1.10 \times 10^{-3} \cdot Re^{0.28}$ |
| $Re \geq 11500$ | $f = 1.53 \times 10^{-2}$（常数）|

---

## 综合对比

| 关联式 | 管型 | 参数输入 | 精度 | 适用 Re |
|---|---|---|---|---|
| Churchill | 全范围(光管/粗糙管) | $\varepsilon/D$ | — | 全 Re 连续 |
| Swamee-Jain / Haaland | 光管湍流 | $\varepsilon, D$ | ~1% | $5\times10^3$–$10^8$ |
| Carnavos | 内翅管 | $A_{fa}/A_{fn}, \alpha$ | ±10% | $>10^4$ |
| Schlager 增强因子 | 微翅管 | 无几何参数 | 近似 | $1.8\times10^4$–$3.4\times10^4$ |
| Goto C-tube | 螺旋沟管 | 管型选择 | ±5% | $2\times10^3$–$4.65\times10^4$ |
| Goto W-tube | 人字沟管 | 管型选择 | ±5% | $<1.15\times10^4$ |

## 代码实现建议

### 统一入口函数

```matlab
function f = single_phase_friction(Re, tube_type, params)
    % tube_type: 'smooth', 'microfin', 'spiral_groove', 'herringbone'
    % params  : 结构体，包含各模型所需参数
    %
    % 返回 Fanning 摩擦因子
    % Δp = 2 * f * L * G^2 / (ρ * D_h)

    switch tube_type
        case 'smooth'
            % Churchill (全范围) 或 Swamee-Jain (仅湍流+过渡区插值)
            f = churchill(Re, params.eps_D);

        case 'microfin'
            % Carnavos 模型
            f = 0.046/Re^0.2 * (params.Afa_Afn)^0.5 * ...
                (1/cosd(params.alpha))^0.75;

        case 'spiral_groove'
            % Goto C-tube (已包含过渡区)
            f = goto_friction(Re, 'C');

        case 'herringbone'
            % Goto W-tube (已包含过渡区)
            f = goto_friction(Re, 'W');
    end
end
```

### 子函数

```matlab
function f = churchill(Re, eps_D)
    % Churchill 全范围摩擦因子 (Fanning)
    A = (2.457 * log(1 / ((7/Re)^0.9 + 0.27*eps_D)))^16;
    B = (37530 / Re)^16;
    f_D = 8 * ((8/Re)^12 + 1/(A+B)^1.5)^(1/12);
    f = f_D / 4;
end

function f = goto_friction(Re, tube_type)
    % Goto 2001 螺旋沟/人字沟管 (Fanning)
    switch tube_type
        case 'C'  % spiral groove
            if Re < 6500
                f = 1.47e-4 * Re^0.53;
            elseif Re < 12700
                f = 0.046 * Re^(-0.20);
            else
                f = 1.23e-2 * Re^0.21;
            end
        case 'W'  % herringbone
            if Re < 3900
                f = 2.17e-2 * Re^(-0.08);
            elseif Re < 11500
                f = 1.10e-3 * Re^0.28;
            else
                f = 1.53e-2;
            end
    end
end
```
