# 工质两相压降关联式 — 文献总结（第 1 批）

> 对应关联式：§10 TFRTPDP — 管翅式工质两相压降
>
> 整理日期：2026-05-26
>
> 第 1 批：3 篇核心 + 2 篇补充
> 补充：Cheng-Thome 2008 CO₂压降 + Xu-Fang 2013 冷凝压降 + Park-Hrnjak 2007 实验验证

---

## 通用符号定义

以下符号在各关联式中通用：

| 符号 | 含义 | 单位 | 来源 |
|---|---|---|---|
| $G$ | 总质量流速 (liquid + vapor) | kg/(m²·s) | $G = \dot{m} / A$，输入 |
| $x$ | 干度（气相质量分数） | — | 能量平衡计算，输入 |
| $\rho_L$ | 液相密度 | kg/m³ | REFPROP/物性，$T_{sat}, P$ |
| $\rho_G$ | 气相密度 | kg/m³ | REFPROP/物性，$T_{sat}, P$ |
| $\mu_L$ | 液相动力粘度 | Pa·s | REFPROP/物性，$T_{sat}, P$ |
| $\mu_G$ | 气相动力粘度 | Pa·s | REFPROP/物性，$T_{sat}, P$ |
| $\sigma$ | 表面张力 | N/m | REFPROP/物性，$T_{sat}, P$ |
| $h_{fg}$ | 汽化潜热 | J/kg | REFPROP/物性，$T_{sat}, P$ |
| $D$ | 管内径 | m | 几何输入 |
| $q$ | 热流密度 | W/m² | $q = Q/A$，输入 |
| $p_r$ | 对比压力 | — | $p_r = P_{sat} / P_{crit}$ |
| $f$ | Fanning 摩擦因子 | — | 由 $Re$ 和单相关联式计算 |

### 基础单相量推导

**全液相（以为全部是液体）摩擦压降梯度**：
$$Re_{Lo} = \frac{G D}{\mu_L}, \quad f_{Lo} = \frac{0.079}{Re_{Lo}^{0.25}} \text{(Blasius, Fanning)}$$
$$\left(\frac{dp}{dz}\right)_{Lo} = \frac{2 f_{Lo} G^2}{\rho_L D}$$

**全气相（以为全部是气体）摩擦压降梯度**：
$$Re_{Go} = \frac{G D}{\mu_G}, \quad f_{Go} = \frac{0.079}{Re_{Go}^{0.25}}$$
$$\left(\frac{dp}{dz}\right)_{Go} = \frac{2 f_{Go} G^2}{\rho_G D}$$

**Martinelli 参数 $X_{tt}$**（两相湍流）：
$$X_{tt} = \left(\frac{1-x}{x}\right)^{0.9} \left(\frac{\rho_G}{\rho_L}\right)^{0.5} \left(\frac{\mu_L}{\mu_G}\right)^{0.1}$$

### 两相摩擦压降梯度通用形式

$$\left(\frac{dp}{dz}\right)_{tp} = \phi^2 \cdot \left(\frac{dp}{dz}\right)_{Lo}$$

其中 $\phi^2$ 为两相乘子，各关联式给出不同形式。

---

## 1. Ould Didi, Kattan, Thome (2002) — 两相压降方法评价

**文献**：Ould Didi MB, Kattan N, Thome JR. "Prediction of two-phase pressure gradients of refrigerants in horizontal tubes". *Int. J. Refrigeration*, 25, 935-947, 2002.

**核心贡献**：5 种制冷剂（R134a, R123, R402A, R404A, R502），管径 10.92/12 mm，质量流速 100–500 kg/m²s，干度 0.04–1.0，评价 7 种压降方法。

### 评价结果排名

| 排名 | 方法 | 特征 |
|---|---|---|
| 1 | **Gronnerud (1979)** | 间歇流/分层波状流最佳 |
| 2 | **Müller-Steinhagen and Heck (1986)** | 环状流最佳，整体最优最稳定 |
| 3 | Friedel (1979) | 被广泛引用但仅排第三 |
| 4 | Lockhart-Martinelli (1949) | |
| 5 | Chisholm (1973) | 明显高估 |

### Müller-Steinhagen and Heck (1986) 公式

**适用范围**：蒸发、冷凝、绝热两相流 **均适用**。在 Ould Didi et al. (2002) 的蒸发实验中环状流最优；在 Xu & Fang (2013) 的冷凝评价中宏通道 MARD 16.3% 为最佳。**MSH 本身可在蒸发和冷凝中直接使用**，Xu-Fang (2012) 和 Xu-Fang (2013) 是在 MSH 基础上的专项改进。

**推导逻辑**：干度 $x=0$ 时退化为全液相值 $a$；$x=1$ 时退化为全气相值 $b$。MSH 用线性插值 $G = a + 2(b-a)x$ 再叠加 $(1-x)^{1/3}$ 和 $x^3$ 权重，构造出平滑过渡的 S 形曲线。

| 符号 | 含义 | 单位 | 推导 |
|---|---|---|---|
| $a$ | 全液相压降梯度 | Pa/m | $a = (dp/dz)_{Lo}$ |
| $b$ | 全气相压降梯度 | Pa/m | $b = (dp/dz)_{Go}$ |
| $G$ | 中间插值函数 | Pa/m | $G = a + 2(b-a)x$ |

**公式**：

$$\left(\frac{dp}{dz}\right)_{tp} = G(1-x)^{1/3} + bx^3$$

**代码实现**：
```matlab
function dpdz = muller_steinhagen_heck(G_mf, x, rho_L, rho_G, mu_L, mu_G, D)
    % === 1. 全液相压降梯度 a ===
    Re_Lo = G_mf * D / mu_L;                % 全液相雷诺数
    f_Lo = 0.079 / Re_Lo^0.25;              % Blasius (Fanning)
    a = 2 * f_Lo * G_mf^2 / (rho_L * D);    % dp/dz_Lo

    % === 2. 全气相压降梯度 b ===
    Re_Go = G_mf * D / mu_G;                % 全气相雷诺数
    f_Go = 0.079 / Re_Go^0.25;
    b = 2 * f_Go * G_mf^2 / (rho_G * D);    % dp/dz_Go

    % === 3. MSH 插值 ===
    G_term = a + 2 * (b - a) * x;
    dpdz = G_term * (1 - x)^(1/3) + b * x^3;

    % 注：x=0 → dpdz=a (全液相); x=1 → dpdz=b (全气相)
end
```

---

## 2. Jung & Radermacher (1989) — 环状流沸腾压降关联式

**文献**：Jung DS, Radermacher R. "Prediction of pressure drop during horizontal annular flow boiling of pure and mixed refrigerants". *Int. J. Heat Mass Transfer*, 32(12), 2435-2446, 1989.

**测试条件**：
- 4 种制冷剂（R22, R114, R12, R152a）纯工质和混合物
- 8m 长管段，内径 9.1 mm，环形流
- 质量流速 230–720 kg/m²s，热流 10–45 kW/m²
- 干度 > 20%，对比压力 0.08–0.16

### 两相压降乘子

**推导逻辑**：Jung & Radermacher 采用分离流模型，假设环形流中气液两相速度不同。他们发现全部数据（纯工质 + 混合物）可以仅用 Martinelli 参数 $X_{tt}$ 关联，无额外的混合物成分依赖。

| 符号 | 含义 | 单位 | 推导 |
|---|---|---|---|
| $\phi_{fo}$ | 两相压降乘子（基于总流量为液相） | — | $\phi_{fo}^2 = \Delta p_{tp} / \Delta p_{Lo}$ |
| $X_{tt}$ | Martinelli 参数（两相湍流） | — | $X_{tt} = \left(\frac{1-x}{x}\right)^{0.9} \left(\frac{\rho_G}{\rho_L}\right)^{0.5} \left(\frac{\mu_L}{\mu_G}\right)^{0.1}$ |
| $p_r$ | 对比压力 | — | $p_r = P_{sat} / P_{crit}$，用于简化版 $X_{tt}$ |

**核心关联式**（回归确定，相关系数 0.985）：

$$\phi_{fo} = 3.58 X_{tt}^{-0.735}$$

**全液相加压降梯度的总两相乘子**（包含干度从 0 到 $x$ 的积分效应，注意此处是总压降乘子而非微分梯度）：

$$\phi_{tp}^2 = \phi_{fo}^2 (1-x)^{1.8} = 12.82 X_{tt}^{-1.47} (1-x)^{1.8}$$

### 对比压力简化版 $X_{tt}$（纯工质专用，$0.06 < p_r < 0.7$，±5%）

Cooper (1980) 发现物性组可以仅用 $p_r$ 关联。Jung & Radermacher 回归得到：

$$\frac{\rho_G}{\rho_L} \left(\frac{\mu_L}{\mu_G}\right)^{0.2} = 0.551 p_r^{-0.011} \approx 0.551 p_r^0$$

代入 $X_{tt}$ 表达式后得：
$$X_{tt} = 0.551 \left(\frac{1-x}{x}\right)^{0.9} p_r^{-0.49}$$

此版本只需干度 $x$ 和对比压力 $p_r$，无需气相/液相物性输入。

**精度**：与实验数据均偏差 **8.4%**

**代码实现**：
```matlab
function dpdz = jung_radermacher(G, x, rho_L, rho_G, mu_L, mu_G, D, p_r)
    % === 1. Martinelli 参数 ===
    % 完整版（推荐）
    Xtt = ((1-x)/x)^0.9 * (rho_G/rho_L)^0.5 * (mu_L/mu_G)^0.1;

    % 简化版（纯工质，只需 x 和 p_r）
    % Xtt = 0.551 * ((1-x)/x)^0.9 * p_r^(-0.49);

    % === 2. 两相乘子 ===
    phi_fo = 3.58 * Xtt^(-0.735);

    % === 3. 全液相基底 ===
    Re_Lo = G * D / mu_L;
    f_Lo = 0.079 / Re_Lo^0.25;
    dpdz_Lo = 2 * f_Lo * G^2 / (rho_L * D);  % Fanning

    % === 4. 两相压降 ===
    dpdz = phi_fo^2 * dpdz_Lo;
end
```

---

## 3. Xu & Fang (2012) — 蒸发压降新关联式

**文献**：Xu Y, Fang X. "A new correlation of two-phase frictional pressure drop for evaporating flow in pipes". *Int. J. Refrigeration*, 35, 2039-2050, 2012.

**数据集**：2622 个数据点，15 种制冷剂，Dh = 0.81–19.1 mm，G = 25.4–1150 kg/m²s

### 评价结论（29 种关联式对比）

| 方法 | MARD (全部) | MARD (宏通道) | MARD (微通道) |
|---|---|---|---|
| Müller-Steinhagen & Heck | 28.5% | **26.9%** | 38.9% |
| Friedel | 29.3% | 29.5% | **28.3%** |
| Xu-Fang (新) | **25.2%** | — | — |

### 新关联式

**推导逻辑**：Xu & Fang 构造了一个三部分组成的两相乘子，物理含义明确：
- **第 1 项**（=1）：$x=0$ 时退化为全液相
- **第 2 项**（$\propto X_{tt}^{-0.5}/Bo^{0.74}$）：沸腾的核态效应，$Bo$ 越大（热流越高）压降越大
- **第 3 项**（$\propto x^{1.5}/We_{tp}^{0.08}$）：高干度区的环状流剪切效应，$We_{tp}$ 越大（惯性越强）压降越大

| 符号 | 含义 | 单位 | 推导 |
|---|---|---|---|
| $\phi_{lo}^2$ | 两相乘子（基于全液相） | — | $\phi_{lo}^2 = (dp/dz)_{tp} / (dp/dz)_{Lo}$ |
| $Bo$ | 沸腾数 | — | $Bo = q / (G h_{fg})$，热流/蒸发潜热通量比 |
| $We_{tp}$ | 两相韦伯数 | — | $We_{tp} = G^2 D / (\rho_{tp} \sigma)$，惯性力/表面张力 |
| $\rho_{tp}$ | 两相均相密度 | kg/m³ | $\rho_{tp} = \left(\frac{x}{\rho_G} + \frac{1-x}{\rho_L}\right)^{-1}$ |

**核心公式**：

$$\phi_{lo}^2 = 1 + \frac{1.28}{Bo^{0.74} X_{tt}^{0.5}} + \frac{0.81 x^{1.5}}{We_{tp}^{0.08}}$$

$$\left(\frac{dp}{dz}\right)_{tp} = \phi_{lo}^2 \cdot \left(\frac{dp}{dz}\right)_{Lo}$$

**代码实现**：
```matlab
function dpdz = xufang_2012(G, x, q, rho_L, rho_G, mu_L, mu_G, sigma, hfg, D)
    % === 1. 两相均相密度 ===
    rho_tp = 1 / (x/rho_G + (1-x)/rho_L);

    % === 2. 全液相基底 ===
    Re_Lo = G * D / mu_L;
    f_Lo = 0.079 / Re_Lo^0.25;             % Blasius (Fanning)
    dpdz_Lo = 2 * f_Lo * G^2 / (rho_L * D);

    % === 3. Martinelli 参数 ===
    Xtt = ((1-x)/x)^0.9 * (rho_G/rho_L)^0.5 * (mu_L/mu_G)^0.1;

    % === 4. 无量纲数 ===
    Bo = q / (G * hfg);                     % 沸腾数
    We_tp = G^2 * D / (rho_tp * sigma);    % 韦伯数

    % === 5. 两相乘子 ===
    phi2_lo = 1 + 1.28 / (Bo^0.74 * Xtt^0.5) + 0.81 * x^1.5 / (We_tp^0.08);

    % === 6. 两相压降 ===
    dpdz = phi2_lo * dpdz_Lo;
end
```

---

## 4. Xu & Fang (2013) — 冷凝压降关联式

**文献**：Xu Y, Fang X. "A new correlation of two-phase frictional pressure drop for condensing flow in pipes". *Nuclear Engineering and Design*, 263, 87-96, 2013.

**数据集**：525 数据点，9 种制冷剂，Dh = 0.1–10.07 mm，G = 20–800 kg/m²s

### 评价结果（29 种关联式对比）

| 方法 | MARD (全部) | MARD (宏通道) | MARD (微通道) |
|---|---|---|---|
| Müller-Steinhagen & Heck | 17.9% | **16.3%** | 23.0% |
| Friedel | 27.5% | 32.5% | 15.3% |
| **Xu-Fang 2013** | **19.4%** | 21.9% | **14.7%** |

### 新关联式

**推导逻辑**：在 MSH 基础上引入 Froude 数和 Weber 数，修正冷凝流的微通道效应。满足 $x=0 \to \phi_{lo}^2=1$（全液）、$x=1 \to \phi_{lo}^2=Y^2$（全气）的边界。

$$\phi_{lo}^2 = Y^2 x^3 + (1 - x^{2.59})^{0.632} \left[1 + 2 x^{1.17} (Y^2 - 1) + 0.00775 x^{-0.475} Fr_{tp}^{0.535} We_{tp}^{0.188}\right]$$

| 符号 | 推导 |
|---|---|
| $Y = \sqrt{(dp/dz)_{Go} / (dp/dz)_{Lo}}$ | 全气相/全液相压降比 |
| $Fr_{tp} = G^2 / (g D \rho_{tp}^2)$ | 两相 Froude 数，惯性/重力 |
| $We_{tp} = G^2 D / (\rho_{tp} \sigma)$ | 两相 Weber 数，惯性/表面张力 |

---

## 5. Cheng-Thome (2008) Part I — CO₂ 流型图压降模型

**文献**：Cheng L, Ribatski G, Moreno Quibén J, Thome JR. "New prediction methods for CO₂ evaporation inside tubes: Part I". *IJHMT*, 51, 111-124, 2008.

**核心结论**：
- 评价 6 种方法对 CO₂ 的预测：**所有方法精度都不够好**
- CO₂ 两相压降远低于常规制冷剂（低压气液密度比 + 低表面张力）
- **推荐工程使用 Müller-Steinhagen & Heck 作为 CO₂ 的近似**

---

## 6. Park & Hrnjak (2007) — CO₂/R410A 低压沸腾实验验证

**文献**：Park CY, Hrnjak PS. "CO₂ and R410A flow boiling heat transfer, pressure drop, and flow pattern at low temperatures in a horizontal smooth tube". *IJR*, 30, 166-178, 2007.

**核心结论**：
- **MSH 和 Friedel 预测 CO₂/R410A 在 ±30% 以内**
- 确认了 MSH 作为通用压降关联式的工程适用性

---

## 七方法对比

| 方法 | 适用场景 | 输入参数 | 精度 | 特点 |
|---|---|---|---|---|
| Müller-Steinhagen & Heck | 通用（推荐默认） | $x, G, \rho, \mu$ | MARD 26.9% | 最简单的插值公式 |
| Jung-Radermacher 1989 | 环状流($x>0.2$) | $x, \rho, \mu, p_r$ | ±8.4% | 可仅用 $p_r$ 参数化 |
| Xu-Fang 2012 (蒸发) | 蒸发流微通道 | $x, G, q, \rho, \mu, \sigma, h_{fg}$ | MARD 25.2% | 含 $Bo$, $We_{tp}$ |
| Xu-Fang 2013 (冷凝) | 冷凝流 | $x, G, \rho, \mu, \sigma$ | MARD 19.4% | 含 $Fr_{tp}$, $We_{tp}$ |
| Friedel 1979 | 通用（$\mu_L/\mu_G<1000$） | $x, G, \rho, \mu, \sigma$ | MARD 29.3% | 40k+ 数据点验证 |
| Cheng-Thome 2008 | CO₂ 专用（流型基础） | — | — | 太复杂，不推荐工程用 |
| Park-Hrnjak 2007 | 实验验证参考 | — | — | 确认 MSH 适用性 |
