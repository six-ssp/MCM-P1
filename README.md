# MCM,美赛A题

---

# **第二问**

# 一、第二问在干什么（承接第一问）

**第一问**做的是“能量守恒式”的 SOC 简化模型：

[
\frac{dSOC}{dt}=-\frac{P(t)}{C_{\text{eff}}}
]

它只回答：**给定功耗曲线 (P(t))**，SOC 怎么掉。

**第二问**在此基础上增加“电池输出端电压”这一层物理机制：
同样的功耗 (P(t)) 下，电流 (I(t)) 会变化，而电流会通过**欧姆内阻**和**RC 极化**产生电压跌落，导致端电压 (V_{\text{term}}) 可能提前触发低压关机。

所以第二问的核心是：

> 在给定使用场景（功耗结构、温度、SOH、亮度、网络等）的情况下，
> 同时模拟 (SOC(t))、(V_{\text{term}}(t))、温度 (T(t))，
> 并用三条终止条件定义手机“电量耗尽时间” (t_{\text{empty}})。

终止条件（代码里就是这三条）：

* **SOC 触底**：(SOC \le 5%)
* **电压截止**：(V_{\text{term}} \le 3.0V)
* **过热保护**：(T \ge 45^\circ C)

---

# 二、代码思路

整套代码就是一个“**输入场景参数 → 生成功耗 → 建立等效电路 → ODE 求解 → 判断终止条件 → 出图/统计**”的流程：

## Step 1：定义场景（S1~S4 + S5敏感性）

每个场景包含：

* 环境温度 (T_{amb})
* 电池健康 SOH
* 初始 SOC
* 屏幕亮度 B、屏幕点亮占空比
* CPU 负载
* 网络类型（4G/5G/WiFi）+ 网络占空比
* GPS 是否开启

这一步对应“文档里各场景设定表”。

## Step 2：由场景生成功耗结构 (P(t))

把总功耗拆成模块功耗之和：

[
P(t)=P_{base}+P_{screen}(t)+P_{cpu}+P_{net}(t)+P_{gps}
]

其中 screen 和 net 用“软方波占空比”模拟“开/关”行为，但比硬方波更平滑。

## Step 3：电池等效电路（Thevenin + 2RC）

电池端电压：

[
V_{\text{term}} = V_{oc}(SOC,T) - I R_0(SOC,T,SOH) - V_1 - V_2
]

并用两条 RC 极化电压微分方程：

[
\dot V_1 = -\frac{V_1}{R_1C_1} + \frac{I}{C_1}
\quad,\quad
\dot V_2 = -\frac{V_2}{R_2C_2} + \frac{I}{C_2}
]

## Step 4：承接第一问的 SOC 动力学

仍然是：

[
\dot{SOC} = -\frac{P(t)}{C_{\text{eff}}}
]

注意做到了“保持两问参数一致”（(C_{\text{eff}}=17Wh)）。

## Step 5：加入热模型

用集总热容 + 对流散热：

* 发热来自：(I^2R_0 + I(V_1+V_2))
* 冷却来自：(k(T-T_{amb}))

得到 (T(t)) 的 ODE。

## Step 6：ODE 求解 + 事件终止

用 `solve_ivp(..., events=...)` 一旦触发三终止条件之一就停止，输出：

* (t_{\text{empty}})
* 触发原因（SOC min / Voltage cutoff / Thermal）

## Step 7：画图 + 汇总表 + 蒙特卡洛

* 每个场景出一套解释图（4张）
* 全场景出对比图（5张）
* 导出 CSV 汇总表
* 蒙特卡洛用参数扰动（±10%），给 t_empty 分布（箱线图）

---

# 三、每个生成文件/图是干什么的（对照文档）

下面按输出文件名讲清楚用途。

---

## A. 每个场景都会生成 4 张图（S1/S2/S3/S4/S5a~e）

### 1）`S*_Q2-1_main.png`

**Main Trajectory: SOC & Terminal Voltage（主轨迹图）**

* 左轴：SOC(%)
* 右轴：端电压 (V_{\text{term}})(V)
* 横线：SOC_min=5%、V_cutoff=3.0V
* 竖线：(t_{\text{empty}})
* 角落标注：终止原因（谁先触发）

 这张图对应文档里**“定义 t_empty & 终止条件”**那一段。
它回答的问题是：

> 在该场景下，手机到底是 SOC 先到底，还是电压先掉到 3.0V，还是过热先触发？

---

### 2）`S*_Q2-2_drop.png`

**Resistance & Voltage Drop Decomposition（电压跌落来源分解）**

画 5 条曲线：

* (R_0(t))：欧姆内阻（温度/SOH影响）
* (I R_0)：欧姆压降（最重要）
* (V_1)、(V_2)：两条 RC 极化压降
* (V_1+V_2)：总极化压降

对应文档里 **“等效电路结构 + 两组 RC 支路意义”**。
它回答的问题是：

> 电压掉到截止值主要是因为：
> 内阻变大？电流变大？还是极化效应（V1/V2）变大？

因为 **V1/V2 是动态状态变量，有时间常数（惯性），不是瞬间变化**，所以它们是平滑的曲线，而不是横平竖直。

---

### 3）`S*_Q2-3_power_stack.png`

**Power Decomposition (Stacked)（功耗结构堆叠）**

堆叠模块功耗：

* Base
* Screen
* CPU
* Network
* GPS

对应文档里 **“功耗结构/各模块耗电”**。
它回答的问题是：

> 这个场景耗电主要来自谁？（屏幕？CPU？网络？）

所以这张图应该更像“连续堆叠面积”，可写报告。

---

### 4）`S*_Q2-4_energy_pie.png`

**Energy Share Pie（能量占比饼图）**

把每个模块对总能量消耗的贡献算出来：

[
E_k=\int P_k(t),dt
]

然后画饼图显示百分比。

 这张图是为了写报告更直观：
堆叠图适合看随时间变化，饼图适合给“总结性结论”。

它对应文档里“功耗分解”那部分的**总结表达**，回答：

> 哪个模块是该场景寿命变短的主因？


## B. 全场景对比图（不分场景）

### 5）`Q2-5_energy_3Dbar.png`

**3D Energy Bar（场景×模块×能量）**

横轴：场景
纵轴：模块
z轴：能量(Wh)

对应文档里“跨场景比较功耗结构差异”。
它回答：

> 哪个场景能量消耗最高？哪个模块在不同场景中变化最大？

---

### 6）`Q2-6_t_empty_comparison.png`

**Battery Life Comparison（寿命对比柱状）**

* x：场景
* y：t_empty (hours)
* 柱子颜色：终止原因（电压/ SOC / 热）

对应文档第二问最核心输出：**各场景电池寿命**。
是论文结论图里最重要的一张。

---

### 7）`Q2-8_shutdown_reason_share.png`

**Shutdown Reason Share（终止原因占比饼图）**

统计所有场景中：

* 电压截止导致关机占比
* SOC最小导致占比
* 过热导致占比

对应文档里“哪个因素更常成为限制”的讨论。

---

### 8）`Q2-9_S5_tornado.png`

**Sensitivity Tornado（敏感性分析）**

S5 系列是对 S1 进行“单因素扰动”：

* 改环境温度
* 改 SOH
* 改亮度
* 改网络等

然后画：

[
\Delta t_{empty}=t_{empty}(S5)-t_{empty}(S1)
]

对应文档里“敏感性分析：哪个参数最影响寿命”。

---

### 9）`Q2-10_MC_boxplot.png`

**Monte Carlo Uncertainty（蒙特卡洛不确定性）**

每个场景做 N=500 次参数扰动（±10%），得到 (t_{empty}) 的分布，画箱线图。

对应文档里“蒙特卡洛”部分。
它回答：

> 在参数不确定的情况下，我们预测的 t_empty 稳不稳？置信区间大概多宽？

---

## C. 表格输出

### 10）`Q2_summary_table.csv`

这是第二问最终的“数据表”。

每一行是一个场景，包含：

* 场景参数（SOC0、T_amb、SOH…）
* 仿真结果（t_empty、终止原因、末端 SOC/V/T）
* 蒙特卡洛区间（p10、p50、p90）

---

# 四、你怎么在论文里写“每图说明”（给你一句话模板）注意：ai给的建议

可以直接这样写图注（我按每图给一句标准图注）：

* **Q2-1**: “Time trajectories of SOC and terminal voltage under scenario S*, with shutdown thresholds and the resulting battery life (t_{empty}).”
* **Q2-2**: “Decomposition of voltage losses into ohmic drop (I R_0) and polarization drops (V_1, V_2), illustrating the physical causes of voltage cutoff.”
* **Q2-3**: “Stacked power decomposition showing time-varying contributions of base, screen, CPU, network, and GPS.”
* **Q2-4**: “Energy share by module, computed via time integration of module power consumption.”
* **Q2-6**: “Comparison of battery life across scenarios; bar colors indicate the shutdown mechanism.”
* **Q2-10**: “Monte Carlo uncertainty of predicted battery life (t_{empty}) under ±10% parameter perturbations.”

---
