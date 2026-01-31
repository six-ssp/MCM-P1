# Problem 2 (Q2) Deliverables Guide (Aligned with the Solution Document)

> 逻辑严格承接 Q1：在给定功耗 \(P(t)\) 情况下，SOC 满足连续时间模型  
> \[
\frac{dSOC}{dt}=-\frac{P(t)}{C_{\text{eff}}},\quad C_{\text{eff}}=17.0\text{ Wh}
\]
> Q2 在此基础上引入 **等效电路(ECM)** + **两 RC 极化** + **热模型**，并用**终止条件**定义手机“电池耗尽时间” \(t_{\text{empty}}\)。

---

## 1. Q2 在问什么（文档对应的核心目标）

Q2 的核心是：在不同使用场景下（S1–S4）以及敏感性场景（S5），基于文档给定的建模假设：

1. 用**功耗分解模型**构造总功耗 \(P(t)\)（Base/Screen/CPU/Network/GPS）。
2. 承接 Q1 用 \(dSOC/dt=-P/C_{\text{eff}}\) 得到 SOC 下降过程。
3. 用 **等效电路**（开路电压 \(V_{oc}(SOC,T)\)、欧姆内阻 \(R_0(SOC,T,SOH)\)、两组并联 RC 极化支路）得到端电压：
   \[
   V_{\text{term}}(t)=V_{oc}(SOC,T)-I(t)R_0 - V_1(t)-V_2(t)
   \]
4. 再引入简化热模型得到电池温度 \(T(t)\)。
5. 定义终止条件（文档阈值）：
   - \(SOC \le 5\%\)
   - \(V_{\text{term}}\le 3.0\text{ V}\)
   - \(T\ge 45^\circ C\)
6. 三个条件**谁先触发谁终止**，对应的终止时刻就是
   \[
   t_{\text{empty}}=\min\{t_{SOC}, t_V, t_T\}
   \]

最终输出：**各场景电池寿命 \(t_{\text{empty}}\)、终止原因、关键轨迹曲线、以及不确定性（蒙特卡洛）**。

---

## 2. 代码总体思路（可以写进“Method”部分）

代码整体分为 5 个模块，完全贴合文档叙事：

### (A) 场景输入（Scenario）
每个场景由一组参数描述（环境温度、健康度 SOH、初始 SOC、亮度、屏幕/网络占空比、CPU负载、GPS 等），例如：
- S1：日常基准
- S2：游戏高压
- S3：导航+5G（含 GPS）
- S4：低温+老化
- S5：基于 S1 的单变量敏感性扰动（温度、SOH、亮度、网络等）

### (B) 功耗结构 → 得到 \(P(t)\)
按文档“模块化功耗结构”构造
\[
P(t)=P_{base}+P_{screen}(t)+P_{cpu}+P_{net}(t)+P_{gps}
\]
其中 Screen/Network 用**平滑占空波**（soft duty wave）模拟“开/关”的负载行为，使图更可读、更像交稿的工程曲线。

### (C) Q1 承接：SOC 动力学
\[
\frac{dSOC}{dt}=-\frac{P(t)}{C_{\text{eff}}}
\]
这一步决定“电量随时间消耗”的主趋势。

### (D) Q2 核心：等效电路 ECM + 2RC 极化 + 热模型
- \(V_{oc}=f(SOC,T)\)：随 SOC 和温度变化的开路电压
- \(R_0=g(T,SOH)\)：温度越低/老化越严重，内阻越高
- 两个 RC 分支电压状态 \(V_1,V_2\)：描述快/慢极化过程
- 端电压 \(V_{\text{term}}\) 与电流 \(I(t)=P(t)/V(t)\)耦合
- 热模型：根据 \(I^2R_0\) 和极化耗散估计升温，并带散热项回到环境温度

### (E) 事件终止（求 \(t_{\text{empty}}\)）
通过 ODE 事件（event）机制检测三条阈值：
- SOC_min、V_cutoff、T_max  
谁先触发，谁就是终止原因，并输出 \(t_{\text{empty}}\)。

### (F) 蒙特卡洛不确定性（文档要求）
对关键参数做 ±10% 随机扰动（例：\(R_0\)基准、温度敏感系数等），重复 N=500 次，得到 \(t_{\text{empty}}\) 的统计分布（P10/P50/P90）。

---

## 3. 每张图是干什么的（逐图解释 + 与文档对应）

下面按文件名解释，报告里可以用“图 X”引用。

---

### **(1) `{S}_Q2-1_main.png`**
**名称：Main Trajectory: SOC & Terminal Voltage (with thresholds & t_empty)**  
**作用：Q2 的主结果图**

- 左轴：SOC (%) 随时间衰减（承接 Q1）
- 右轴：端电压 \(V_{\text{term}}(t)\)（ECM 输出）
- 虚线阈值：
  - \(SOC_{min}=5\%\)
  - \(V_{cutoff}=3.0V\)
- 竖线：\(t_{\text{empty}}\)（最先触发的终止时刻）
- 文本框：打印终止原因（SOC/Voltage/Thermal）

**文档对应：**
- “SOC ODE 动力学 + 终止条件定义 \(t_{empty}\)”
- “端电压模型用于判定电压截止”

---

### **(2) `{S}_Q2-2_drop.png`**
**名称：Voltage Drop Decomposition (ECM)**  
**作用：解释端电压为何会下降、下降由哪些电压损失贡献**

包含 5 条曲线：
- \(R_0(t)\)：欧姆内阻（体现温度/老化影响）
- \(I(t)R_0(t)\)：欧姆压降（直接降低端电压）
- \(V_1(t)\)：快速极化电压
- \(V_2(t)\)：慢速极化电压
- \(V_1(t)+V_2(t)\)：总极化压降

**为什么需要它：**
当某场景是“电压先截止”时，这张图可以说明：
- 是 **I 过大导致 IR0 过大**？
- 还是 **极化累积 V1+V2 变大**？
- 或者低温/老化导致 **R0 增大**？

**文档对应：**
- “等效电路结构：Voc + R0 + 两组RC极化”
- “端电压由电压损失项决定”

---

### **(3) `{S}_Q2-3_power_stack.png`**
**名称：Power Decomposition (Stacked, Readable)**  
**作用：展示功耗结构随时间的组成（解释能量消耗来自哪里）**

堆叠面积图展示：
- Base / Screen / CPU / Network / GPS
- 总高度就是瞬时 \(P(t)\)

**为什么有用：**
- 直接体现“高耗电模块”
- 解释 S2（游戏）为什么寿命短：CPU + Screen 占比大
- 解释 S3 为什么耗电快：Network + GPS 占比明显

**文档对应：**
- “功耗结构模型：按模块叠加得到 P(t)”

---

### **(4) `{S}_Q2-4_energy_pie.png`**
**名称：Energy Share by Module (Pie)**  
**作用：把“功耗堆叠”转换为更直观的“能量占比”总结图**

- 计算各模块能量：
  \[
  E_k=\int P_k(t)\,dt
  \]
- 饼图展示各模块对总能量消耗的贡献比例
- 图下方显示总能量（Wh）

**为什么需要它：**
堆叠图适合看时间结构，但评委更喜欢“一眼结论”。饼图能直接说：
- “屏幕消耗占 XX%”
- “网络消耗占 XX%”
便于写结论段落。

**文档对应：**
- 功耗模型的“能量累积解释”，用于场景差异对比

---

## 4. 全场景对比图（结论与分析用）

### **(5) `Q2-5_energy_3Dbar.png`**
**名称：Cross-Scenario Energy Comparison (3D Bars)**  
**作用：跨场景对比“不同模块消耗的能量”**

- X：场景（S1–S5）
- Y：模块（Base/Screen/CPU/Network/GPS）
- Z：能量（Wh）

**文档对应：**
- 场景对比分析（S1–S4）+ S5 敏感性对比

---

### **(6) `Q2-6_t_empty_comparison.png`**
**名称：Battery Life Comparison: t_empty and Shutdown Reason**  
**作用：Q2 的最终对比结论图**

- 柱高：各场景的 \(t_{\text{empty}}\)
- 柱颜色：终止原因（SOC/Voltage/Thermal）
  - 直观回答“哪个场景最耐用/最危险/最容易电压截止”

**文档对应：**
- “定义并比较 t_empty”
- “终止条件触发的分类结果”

---

### **(7) `Q2-8_shutdown_reason_share.png`**
**名称：Shutdown Reason Share (Pie)**  
**作用：总结所有场景中“主要失败模式”占比**

- 比如电压截止占比更高说明“电压是更常见的限制因素”
- 用于写总结段落：系统最常见的失效机制是什么

**文档对应：**
- “终止机制统计总结”

---

### **(8) `Q2-9_S5_tornado.png`**
**名称：Sensitivity (S5) Tornado: Δt_empty relative to S1**  
**作用：敏感性分析（S5）——哪些因素对寿命影响最大**

- 基准：S1 的 \(t_{\text{empty}}\)
- 横轴：Δt_empty（相对 S1 的变化）
- 每条 bar：一个单变量扰动（温度、SOH、亮度、网络等）

**文档对应：**
- S5 单变量敏感性实验（“改变一项参数，观察寿命变化”）
- 用于回答：寿命对温度/老化/屏幕等哪个更敏感

---

### **(9) `Q2-10_MC_boxplot.png`**
**名称：Monte Carlo Uncertainty: t_empty Distribution**  
**作用：蒙特卡洛不确定性分析（文档要求的鲁棒性）**

- 每个场景一个箱线图，展示 \(t_{\text{empty}}\) 的分布范围
- 说明模型对参数扰动的稳定性
- 对应 CSV 里的 P10/P50/P90

**文档对应：**
- “Monte Carlo N=500，±10% 扰动”要求

---

## 5. 生成的表格文件是干什么的？

### **(10) `Q2_summary_table.csv`**
**作用：Q2 的“可交稿结构化结果表”**，用于：
- 直接放入论文 Table（可复制到 Word）
- 方便计算/排序/绘制更多统计图
- 作为“可复现结果输出”的证明

典型字段含义：
- `t_empty(h)`：该场景电池寿命（小时）
- `ShutdownReason`：终止原因（SOC / Voltage / Thermal）
- `SOC_end(%)`、`Vterm_end(V)`、`T_end(C)`：终止时刻的状态
- `MC_p10/p50/p90(h)`：蒙特卡洛分位数（10%/50%/90%）

---

## 6. 图应该放在报告的哪个位置（建议排版）

建议 Q2 报告按“主线叙事”摆图：

1. **方法描述（Model Section）**
   - 简述功耗结构 + ECM + 终止条件
2. **单场景解释（Results per Scenario）**
   - 放 `S1_Q2-1_main.png`（主图）
   - 紧跟 `S1_Q2-2_drop.png`（解释电压为什么掉）
   - 再放 `S1_Q2-4_energy_pie.png`（能量解释）
   - `S1_Q2-3_power_stack.png`（作为补充/附录也可）
3. **跨场景比较（Cross-Scenario Comparison）**
   - `Q2-6_t_empty_comparison.png`（核心结论）
   - `Q2-5_energy_3Dbar.png`（解释差异来源）
   - `Q2-9_S5_tornado.png`（敏感性结论）
   - `Q2-10_MC_boxplot.png`（鲁棒性/不确定性）
4. **总结段**
   - `Q2-8_shutdown_reason_share.png`（一句话总结失效模式）

---

## 7. 可以直接用的“简短图注模板”（ai分析建议）

- **Fig. Q2-1.** SOC and terminal voltage trajectories under scenario {S}, with cutoff thresholds and computed battery life \(t_{\text{empty}}\).  
- **Fig. Q2-2.** Decomposition of terminal voltage losses into ohmic drop \(I R_0\) and polarization components \(V_1, V_2\).  
- **Fig. Q2-3.** Stacked power consumption by module over time, forming the total load \(P(t)\).  
- **Fig. Q2-4.** Energy consumption share by module obtained by integrating power over time.  
- **Fig. Q2-6.** Cross-scenario comparison of \(t_{\text{empty}}\) and dominant shutdown reason.  
- **Fig. Q2-9.** Tornado plot showing sensitivity of battery life to single-parameter perturbations (S5).  
- **Fig. Q2-10.** Monte Carlo uncertainty of \(t_{\text{empty}}\) under ±10% parameter perturbations (N=500).

---
