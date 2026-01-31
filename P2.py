# -*- coding: utf-8 -*-
"""
Q2 交稿版（整体优化 + 多种图 + 进度条 + legend不遮挡）
- 严格承接Q1：SOC ODE dSOC/dt = -P(t)/C_eff_Wh，且 C_eff_Wh=17.0
- Q2：等效电路 + 两个RC极化 + 热模型 + 三终止条件（SOC/电压/温度）
- 输出图（每场景）：
  Q2-1 main（SOC & Vterm，含阈值线+终止竖线+原因，legend在右侧更靠外）
  Q2-2 drop（R0、I*R0、V1、V2、V1+V2）
  Q2-3 power stack（下采样后更可读）
  Q2-4 energy pie（能量占比）
- 输出图（全场景）：
  Q2-5 3D energy bar
  Q2-6 t_empty comparison（颜色=终止原因）
  Q2-8 shutdown reason share（饼图）
  Q2-9 S5 tornado（敏感性）
  Q2-10 MC boxplot（蒙特卡洛箱线图）
- 输出表：
  Q2_summary_table.csv

你提的三点：
(1) main 图图示往右调：已做（right留白更大 + bbox_to_anchor更靠右）
(2) event_timeline 全挂：默认不画（MAKE_EVENT_TIMELINE=False）；想要可开
(3) MC 加进度条：tqdm，有则显示，无则自动降级
"""

import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import Dict, Tuple, Optional

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from scipy.integrate import solve_ivp

# ---------------------------------------------------------
# 可选：进度条（没装tqdm也不影响运行）
# ---------------------------------------------------------
try:
    from tqdm import tqdm
except Exception:
    tqdm = None

# =========================================================
# 0) 全局开关
# =========================================================
RUN_MC = True
MC_N = 500                 # 文档要求：N=500
MAKE_EVENT_TIMELINE = False  # 你说这图挂了：默认不画（最稳）

# =========================================================
# 1) 画图风格：英文 + Times New Roman + 无网格
# =========================================================
def setup_plot_style():
    mpl.rcParams["font.family"] = "Times New Roman"
    mpl.rcParams["axes.unicode_minus"] = False
    mpl.rcParams["figure.dpi"] = 150
    mpl.rcParams["savefig.dpi"] = 300
    mpl.rcParams["axes.titlesize"] = 13
    mpl.rcParams["axes.labelsize"] = 12
    mpl.rcParams["xtick.labelsize"] = 11
    mpl.rcParams["ytick.labelsize"] = 11
    mpl.rcParams["legend.fontsize"] = 11

setup_plot_style()

def beautify_axes(ax):
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=4, width=1)

def legend_outside_right_fig(fig, handles, labels, x=0.82):
    """
    图级legend放在右侧留白区（更靠右就把x调大）
    """
    fig.legend(handles, labels, loc="center left",
               bbox_to_anchor=(x, 0.5), frameon=False)

# =========================================================
# 2) 承接Q1：容量参数一致
# =========================================================
C_eff_Wh = 17.0

# =========================================================
# 3) 文档阈值：终止条件
# =========================================================
SOC_min = 0.05
V_cutoff = 3.0
T_max = 45.0

# =========================================================
# 4) 场景（按你docx的设定）
# =========================================================
@dataclass
class Scenario:
    name: str
    T_amb: float
    SOH: float
    SOC0: float
    B: float
    screen_duty: float
    U_cpu: float
    network: str
    net_duty: float
    gps_on: bool

def build_scenarios() -> Dict[str, Scenario]:
    S1 = Scenario("S1 (Daily Baseline)", 25.0, 1.0, 0.80, 0.50, 0.40, 0.30, "4G", 0.25, False)
    S2 = Scenario("S2 (Gaming Stress)", 25.0, 1.0, 0.80, 1.00, 1.00, 0.80, "WiFi", 1.00, False)
    S3 = Scenario("S3 (Nav + 5G)", 35.0, 0.90, 0.90, 0.70, 1.00, 0.50, "5G", 1.00, True)
    S4 = Scenario("S4 (Cold + Aged)", 5.0, 0.80, 0.80, 0.50, 0.40, 0.30, "4G", 0.25, False)

    # S5：敏感性（基于S1）
    S5a = Scenario("S5a (T=0C)", 0.0, 1.0, 0.80, 0.50, 0.40, 0.30, "4G", 0.25, False)
    S5b = Scenario("S5b (T=40C)", 40.0, 1.0, 0.80, 0.50, 0.40, 0.30, "4G", 0.25, False)
    S5c = Scenario("S5c (SOH=0.85)", 25.0, 0.85, 0.80, 0.50, 0.40, 0.30, "4G", 0.25, False)
    S5d = Scenario("S5d (B=100%)", 25.0, 1.0, 0.80, 1.00, 0.40, 0.30, "4G", 0.25, False)
    S5e = Scenario("S5e (No Network)", 25.0, 1.0, 0.80, 0.50, 0.40, 0.30, "4G", 0.00, False)

    return {"S1": S1, "S2": S2, "S3": S3, "S4": S4, "S5a": S5a, "S5b": S5b, "S5c": S5c, "S5d": S5d, "S5e": S5e}

SCENARIOS = build_scenarios()

# =========================================================
# 5) 数值积分兼容（避免trapz警告）
# =========================================================
def integrate_trap(y, x):
    """
    兼容写法：新NumPy用 trapezoid，旧NumPy回退 trapz
    """
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(y, x))
    return float(np.trapz(y, x))

# =========================================================
# 6) 功耗结构（文档一致）+ 软方波占空比（交稿更可读）
# =========================================================
P_idle = 0.35
k_screen = 2.6
k_cpu = 3.2
k_gps = 0.8
NET_POWER = {"4G": 1.1, "5G": 1.8, "WiFi": 0.9}

def duty_wave_soft(t_h, duty, period=5.0/60.0, smooth=0.03):
    """
    平滑占空比：避免硬方波带来的锯齿误解（不改变趋势/终止原因）
    """
    if duty <= 0.0:
        return 0.0
    if duty >= 1.0:
        return 1.0
    phase = (t_h % period) / period
    k = 1.0 / max(smooth, 1e-3)
    up = 0.5 * (1 + np.tanh(k * (phase - 0.0)))
    down = 0.5 * (1 - np.tanh(k * (phase - duty)))
    return float(np.clip(up * down, 0.0, 1.0))

def power_components(t_h, sc: Scenario):
    I_screen = duty_wave_soft(t_h, sc.screen_duty)
    I_net = duty_wave_soft(t_h, sc.net_duty)
    return {
        "Base": P_idle,
        "Screen": k_screen * sc.B * I_screen,
        "CPU": k_cpu * sc.U_cpu,
        "Network": NET_POWER.get(sc.network, 1.1) * I_net,
        "GPS": k_gps if sc.gps_on else 0.0,
    }

# =========================================================
# 7) 等效电路：Voc、R0、两RC
# =========================================================
def Voc(SOC, T_c):
    soc = np.clip(SOC, 0.0, 1.0)
    v = 3.0 + 1.2*(1/(1+np.exp(-8*(soc-0.35))))
    v += 0.002*(T_c - 25.0)
    return float(v)

def R0(T_c, SOH, R0_ref=0.08, Ea_R=1800.0):
    T_k = (T_c + 273.15)
    Tref_k = 298.15
    temp_factor = np.exp(Ea_R*(1/T_k - 1/Tref_k))
    return float(R0_ref * temp_factor / max(SOH, 0.2))

def rc_params(T_c):
    # 参考参数 + 温度调节（保持你之前文档思路）
    R1, C1 = 0.015, 2000.0
    R2, C2 = 0.035, 6000.0
    factor = 1.0 + 0.01*(25.0 - T_c)
    factor = np.clip(factor, 0.8, 1.6)
    return float(R1*factor), float(C1), float(R2*factor), float(C2)

# =========================================================
# 8) 热模型（集总）
# =========================================================
C_th = 180.0
k_cool = 8.0

# =========================================================
# 9) ODE系统 + 事件终止
# =========================================================
def rhs(t_h, y, sc: Scenario, params: dict):
    SOC, V1, V2, T_c, SOH = y
    comp = power_components(t_h, sc)
    P = float(sum(comp.values()))

    voc = Voc(SOC, T_c)
    r0 = R0(T_c, SOH, R0_ref=params["R0_ref"], Ea_R=params["Ea_R"])
    R1, C1, R2, C2 = rc_params(T_c)

    v_est = max(voc - V1 - V2, 3.0)
    I = P / v_est

    # 承接Q1：SOC ODE
    dSOC_dt = -P / C_eff_Wh

    # 两RC极化（秒->小时）
    tau1_h = (R1*C1) / 3600.0
    tau2_h = (R2*C2) / 3600.0
    dV1_dt = -(V1 / tau1_h) + (I / C1) * 3600.0
    dV2_dt = -(V2 / tau2_h) + (I / C2) * 3600.0

    # 热：欧姆热 + 极化耗散 - 对流散热
    P_heat = (I**2)*r0 + I*V1 + I*V2
    dT_dt = (P_heat / C_th)*3600.0 - (k_cool*(T_c - sc.T_amb)/C_th)*3600.0

    dSOH_dt = 0.0
    return [dSOC_dt, dV1_dt, dV2_dt, dT_dt, dSOH_dt]

def make_events(sc: Scenario, params: dict):
    def ev_soc(t_h, y): return y[0] - SOC_min
    ev_soc.terminal = True; ev_soc.direction = -1

    def ev_vcut(t_h, y):
        SOC, V1, V2, T_c, SOH = y
        voc = Voc(SOC, T_c)
        r0 = R0(T_c, SOH, R0_ref=params["R0_ref"], Ea_R=params["Ea_R"])
        P = float(sum(power_components(t_h, sc).values()))
        v_est = max(voc - V1 - V2, 3.0)
        I = P / v_est
        Vterm = voc - I*r0 - V1 - V2
        return Vterm - V_cutoff
    ev_vcut.terminal = True; ev_vcut.direction = -1

    def ev_th(t_h, y): return T_max - y[3]
    ev_th.terminal = True; ev_th.direction = -1

    return [ev_soc, ev_vcut, ev_th]

def simulate_one(sc: Scenario, t_max_h=24.0, n_eval=2200, params=None):
    if params is None:
        params = {"R0_ref": 0.08, "Ea_R": 1800.0}

    y0 = [sc.SOC0, 0.0, 0.0, sc.T_amb, sc.SOH]
    t_eval = np.linspace(0, t_max_h, n_eval)

    sol = solve_ivp(
        fun=lambda tt, yy: rhs(tt, yy, sc, params),
        t_span=(0, t_max_h),
        y0=y0,
        t_eval=t_eval,
        method="BDF",
        rtol=1e-6,
        atol=1e-8,
        events=make_events(sc, params),
    )

    t_end = float(sol.t[-1])
    reason = "Reached t_max"
    if sol.status == 1 and sol.t_events:
        te = [ev[0] if len(ev) else np.inf for ev in sol.t_events]
        k = int(np.argmin(te))
        reason = ["SOC_min (5%)", "Voltage cutoff (3.0V)", "Thermal shutdown (45C)"][k]

    SOC, V1, V2, T_c, SOH = sol.y
    t = sol.t

    # 后处理：I, Vterm, 以及功耗分解序列
    comp_series = {k: [] for k in ["Base", "Screen", "CPU", "Network", "GPS"]}
    I_list, Vt_list, R0_list = [], [], []
    V1_list, V2_list = [], []

    for tt, soc, v1, v2, temp, soh in zip(t, SOC, V1, V2, T_c, SOH):
        comp = power_components(tt, sc)
        P = float(sum(comp.values()))
        voc = Voc(soc, temp)
        r0 = R0(temp, soh, R0_ref=params["R0_ref"], Ea_R=params["Ea_R"])
        v_est = max(voc - v1 - v2, 3.0)
        I = P / v_est
        Vterm = voc - I*r0 - v1 - v2

        for kk in comp_series.keys():
            comp_series[kk].append(comp[kk])

        I_list.append(I); Vt_list.append(Vterm); R0_list.append(r0)
        V1_list.append(v1); V2_list.append(v2)

    return {
        "t": t,
        "SOC": SOC,
        "T": T_c,
        "SOH": SOH,
        "I": np.array(I_list),
        "Vterm": np.array(Vt_list),
        "R0": np.array(R0_list),
        "V1": np.array(V1_list),
        "V2": np.array(V2_list),
        "comp": {k: np.array(v) for k, v in comp_series.items()},
        "t_empty": t_end,
        "reason": reason,
    }

# =========================================================
# 10) 平滑/下采样（让图可读）
# =========================================================
def moving_average(x: np.ndarray, win: int) -> np.ndarray:
    if win <= 1:
        return x
    kernel = np.ones(int(win)) / float(win)
    return np.convolve(x, kernel, mode="same")

# =========================================================
# 11) 单场景出图（Q2-1~Q2-4）
# =========================================================
def plot_bundle(res, sc_key: str):
    t = res["t"]
    SOC = res["SOC"]
    Vt = res["Vterm"]
    comp = res["comp"]

    # ---------------- Q2-1 Main ----------------
    fig, ax1 = plt.subplots(figsize=(10.8, 4.9))

    # 关键：留更大右边空白（你要“往右边调”就调这个）
    fig.subplots_adjust(right=0.68)  # 原来0.72更紧；改小=留白更大

    ln_soc, = ax1.plot(t, 100*SOC, linewidth=3.0, label="SOC (%)")
    ln_soc_min = ax1.axhline(100*SOC_min, color="#666666", linestyle="--", linewidth=1.2, label="SOC_min = 5%")
    ax1.set_xlabel("Time t (hours)")
    ax1.set_ylabel("SOC (%)")
    ax1.set_ylim(0, 105)
    beautify_axes(ax1)

    ax2 = ax1.twinx()
    Vt_s = moving_average(Vt, win=max(3, len(Vt)//400))
    ln_vt, = ax2.plot(t, Vt_s, linewidth=2.8, label="V_term (V)")
    ln_vcut = ax2.axhline(V_cutoff, color="#666666", linestyle="--", linewidth=1.2, label="V_cutoff = 3.0V")
    ax2.set_ylabel("")  # 避免和legend冲突
    ax2.set_ylim(2.6, 4.4)
    ax2.grid(False)

    # 终止竖线 + 原因标注（写论文很好用）
    ax1.axvline(res["t_empty"], color="#444444", linestyle=":", linewidth=1.6)
    ax1.text(0.02, 0.05,
             f"t_empty = {res['t_empty']:.2f} h\nReason: {res['reason']}",
             transform=ax1.transAxes, fontsize=11, va="bottom")
    ax1.text(0.98, 0.98, "Right axis: V_term (V)", transform=ax1.transAxes,
             ha="right", va="top", fontsize=11)

    ax1.set_title(f"{sc_key} Q2-1 Main Trajectory (Thresholds & t_empty)")

    handles = [ln_soc, ln_soc_min, ln_vt, ln_vcut]
    labels = [h.get_label() for h in handles]
    # 关键：legend更往右（你要再右就把0.82改成0.84/0.86）
    legend_outside_right_fig(fig, handles, labels, x=0.80)

    fig.savefig(f"{sc_key}_Q2-1_main.png")
    plt.close(fig)

    # ---------------- Q2-2 Drop ----------------
    I = res["I"]
    R0_arr = res["R0"]
    V1 = res["V1"]
    V2 = res["V2"]
    IR0 = I * R0_arr
    Vpol = V1 + V2

    win = max(3, len(t)//400)
    IR0_s = moving_average(IR0, win)
    V1_s = moving_average(V1, win)
    V2_s = moving_average(V2, win)
    Vpol_s = moving_average(Vpol, win)

    fig, ax = plt.subplots(figsize=(10.8, 4.9))
    fig.subplots_adjust(right=0.78)
    ax.plot(t, R0_arr, linewidth=2.3, label="R0 (Ohm)")
    ax.plot(t, IR0_s, linewidth=2.3, label="I*R0 (V)")
    ax.plot(t, V1_s, linewidth=2.1, label="V1 (V)")
    ax.plot(t, V2_s, linewidth=2.1, label="V2 (V)")
    ax.plot(t, Vpol_s, linewidth=2.7, label="V1+V2 (V)")
    ax.set_title(f"{sc_key} Q2-2 Voltage Drop Decomposition (ECM)")
    ax.set_xlabel("Time t (hours)")
    ax.set_ylabel("Magnitude")
    beautify_axes(ax)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    fig.savefig(f"{sc_key}_Q2-2_drop.png")
    plt.close(fig)

    # ---------------- Q2-3 Power Stack（下采样避免竖条） ----------------
    step = max(1, len(t)//600)
    t_ds = t[::step]
    comp_ds = {k: comp[k][::step] for k in comp.keys()}

    fig, ax = plt.subplots(figsize=(10.8, 4.9))
    fig.subplots_adjust(right=0.78)
    ax.stackplot(
        t_ds,
        comp_ds["Base"], comp_ds["Screen"], comp_ds["CPU"], comp_ds["Network"], comp_ds["GPS"],
        labels=["Base", "Screen", "CPU", "Network", "GPS"],
        alpha=0.90
    )
    ax.set_title(f"{sc_key} Q2-3 Power Decomposition (Stacked, Readable)")
    ax.set_xlabel("Time t (hours)")
    ax.set_ylabel("Power (W)")
    beautify_axes(ax)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    fig.savefig(f"{sc_key}_Q2-3_power_stack.png")
    plt.close(fig)

    # ---------------- Q2-4 Energy Pie ----------------
    E = {k: integrate_trap(comp[k], t) for k in comp.keys()}
    labels_p = list(E.keys())
    values = np.array([E[k] for k in labels_p], dtype=float)
    total = float(values.sum()) if values.sum() > 0 else 1.0

    fig, ax = plt.subplots(figsize=(7.4, 5.2))
    ax.set_title(f"{sc_key} Q2-4 Energy Share by Module")
    ax.pie(values, labels=labels_p,
           autopct=lambda p: f"{p:.1f}%", startangle=90, counterclock=False)
    ax.text(0.0, -1.25, f"Total energy ≈ {total:.2f} Wh", ha="center", fontsize=11)
    fig.savefig(f"{sc_key}_Q2-4_energy_pie.png")
    plt.close(fig)

    # 可选：event timeline（默认不画）
    if MAKE_EVENT_TIMELINE:
        plot_event_timeline_safe(res, sc_key)

    return E

def plot_event_timeline_safe(res, sc_key: str):
    """
    修复版 event timeline：
    - NaN 改 0 + 显示N/A
    - 不用 bbox_inches='tight'，避免裁成细线
    """
    t = res["t"]
    SOC = res["SOC"]
    Vt = res["Vterm"]
    T = res["T"]

    soc_dist = SOC - SOC_min
    v_dist = Vt - V_cutoff
    t_dist = T_max - T

    def first_cross_time(dist):
        idx = np.where(dist <= 0)[0]
        return float(t[idx[0]]) if len(idx) else np.nan

    tsoc = first_cross_time(soc_dist)
    tv = first_cross_time(v_dist)
    tt = first_cross_time(t_dist)

    items = ["SOC_min", "V_cutoff", "T_max"]
    vals_raw = np.array([tsoc, tv, tt], dtype=float)
    vals_plot = np.nan_to_num(vals_raw, nan=0.0)

    fig, ax = plt.subplots(figsize=(7.8, 4.2))
    bars = ax.bar(items, vals_plot)
    ymax = np.nanmax(vals_raw) if np.isfinite(vals_raw).any() else 0.0
    ymax = max(ymax, res["t_empty"], 1.0)
    ax.set_ylim(0, 1.15*ymax)

    ax.set_title(f"{sc_key} Q2-7 Event Timeline (Which Hits First?)")
    ax.set_ylabel("First Crossing Time (hours)")
    beautify_axes(ax)

    for b, v in zip(bars, vals_raw):
        if np.isfinite(v):
            ax.text(b.get_x()+b.get_width()/2, b.get_height(), f"{v:.2f}",
                    ha="center", va="bottom", fontsize=11)
        else:
            ax.text(b.get_x()+b.get_width()/2, 0.02*ymax, "N/A",
                    ha="center", va="bottom", fontsize=11)

    fig.savefig(f"{sc_key}_Q2-7_event_timeline.png")
    plt.close(fig)

# =========================================================
# 12) 蒙特卡洛：加进度条
# =========================================================
def monte_carlo_t_empty(sc: Scenario, N=500, seed=7):
    rng = np.random.default_rng(seed)
    samples = []

    iterator = range(N)
    if tqdm is not None:
        iterator = tqdm(iterator, desc=f"MC {sc.name}", ncols=90)

    for _ in iterator:
        Ea_R = 1800.0 * (1 + rng.normal(0, 0.10))
        R0_ref = 0.08 * (1 + rng.normal(0, 0.10))
        SOH_pert = sc.SOH * (1 + rng.normal(0, 0.10))

        Ea_R = max(Ea_R, 100.0)
        R0_ref = max(R0_ref, 1e-3)
        SOH_pert = float(np.clip(SOH_pert, 0.3, 1.0))

        sc_mc = Scenario(sc.name, sc.T_amb, SOH_pert, sc.SOC0, sc.B,
                         sc.screen_duty, sc.U_cpu, sc.network, sc.net_duty, sc.gps_on)

        res = simulate_one(sc_mc, t_max_h=24.0, n_eval=500,
                           params={"R0_ref": R0_ref, "Ea_R": Ea_R})
        samples.append(res["t_empty"])

    arr = np.array(samples, dtype=float)
    p10, p50, p90 = np.percentile(arr, [10, 50, 90])
    return {"samples": arr, "p10": float(p10), "p50": float(p50), "p90": float(p90)}

# =========================================================
# 13) 全场景图：3D能量柱、t_empty柱、原因饼、S5 tornado、MC箱线
# =========================================================
def plot_energy_3d(all_energy: Dict[str, Dict[str, float]]):
    keys = list(all_energy.keys())
    modules = ["Base", "Screen", "CPU", "Network", "GPS"]

    Z = np.zeros((len(keys), len(modules)))
    for i, k in enumerate(keys):
        for j, m in enumerate(modules):
            Z[i, j] = all_energy[k].get(m, 0.0)

    fig = plt.figure(figsize=(12.0, 6.8))
    ax = fig.add_subplot(111, projection="3d")

    xpos, ypos = np.meshgrid(np.arange(len(keys)), np.arange(len(modules)), indexing="ij")
    xpos = xpos.ravel(); ypos = ypos.ravel()
    zpos = np.zeros_like(xpos)

    dx = 0.5*np.ones_like(xpos)
    dy = 0.5*np.ones_like(ypos)
    dz = Z.ravel()

    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, shade=True)
    ax.set_title("Q2-5 Cross-Scenario Energy Comparison (3D Bars)")
    ax.set_xlabel("Scenario")
    ax.set_ylabel("Module")
    ax.set_zlabel("Energy (Wh)")
    ax.set_xticks(np.arange(len(keys)) + 0.25)
    ax.set_xticklabels(keys)
    ax.set_yticks(np.arange(len(modules)) + 0.25)
    ax.set_yticklabels(modules)

    fig.savefig("Q2-5_energy_3Dbar.png")
    plt.close(fig)

def plot_tempty_bar(df: pd.DataFrame):
    reason_color = {
        "SOC_min (5%)": "#1f77b4",
        "Voltage cutoff (3.0V)": "#ff7f0e",
        "Thermal shutdown (45C)": "#d62728",
        "Reached t_max": "#7f7f7f",
    }

    keys = df["ScenarioKey"].tolist()
    tvals = df["t_empty(h)"].to_numpy(dtype=float)
    reasons = df["ShutdownReason"].tolist()
    colors = [reason_color.get(r, "#7f7f7f") for r in reasons]

    fig, ax = plt.subplots(figsize=(11.8, 5.2))
    bars = ax.bar(keys, tvals)
    for i, b in enumerate(bars):
        b.set_facecolor(colors[i])

    ax.set_title("Q2-6 Battery Life Comparison: t_empty and Shutdown Reason")
    ax.set_xlabel("Scenario")
    ax.set_ylabel("t_empty (hours)")
    beautify_axes(ax)

    used = []
    handles = []
    labels = []
    for r in reasons:
        if r in used:
            continue
        used.append(r)
        handles.append(plt.Line2D([0], [0], color=reason_color.get(r, "#7f7f7f"), lw=8))
        labels.append(r)
    ax.legend(handles, labels, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)

    fig.savefig("Q2-6_t_empty_comparison.png")
    plt.close(fig)

def plot_reason_share(df: pd.DataFrame):
    counts = df["ShutdownReason"].value_counts()
    fig, ax = plt.subplots(figsize=(7.4, 5.2))
    ax.set_title("Q2-8 Shutdown Reason Share (All Scenarios)")
    ax.pie(counts.values, labels=counts.index,
           autopct=lambda p: f"{p:.1f}%", startangle=90, counterclock=False)
    fig.savefig("Q2-8_shutdown_reason_share.png")
    plt.close(fig)

def plot_s5_tornado(df: pd.DataFrame):
    base_vals = df.loc[df["ScenarioKey"] == "S1", "t_empty(h)"].values
    if len(base_vals) == 0:
        return
    base = float(base_vals[0])

    s5 = df[df["ScenarioKey"].str.startswith("S5")].copy()
    if len(s5) == 0:
        return
    s5["Delta(h)"] = s5["t_empty(h)"] - base
    s5 = s5.sort_values("Delta(h)")

    fig, ax = plt.subplots(figsize=(10.0, 4.8))
    ax.barh(s5["ScenarioKey"], s5["Delta(h)"])
    ax.axvline(0, color="#555555", linewidth=1.2)
    ax.set_title("Q2-9 Sensitivity (S5) Tornado: Δt_empty relative to S1")
    ax.set_xlabel("Δt_empty (hours)")
    ax.set_ylabel("Scenario")
    beautify_axes(ax)
    fig.savefig("Q2-9_S5_tornado.png")
    plt.close(fig)

def plot_mc_boxplot(mc_samples: Dict[str, np.ndarray]):
    keys = list(mc_samples.keys())
    data = [mc_samples[k] for k in keys]

    fig, ax = plt.subplots(figsize=(11.8, 5.2))
    ax.boxplot(data, labels=keys, showfliers=False)
    ax.set_title("Q2-10 Monte Carlo Uncertainty: t_empty Distribution")
    ax.set_xlabel("Scenario")
    ax.set_ylabel("t_empty (hours)")
    beautify_axes(ax)
    fig.savefig("Q2-10_MC_boxplot.png")
    plt.close(fig)

# =========================================================
# 14) 主流程
# =========================================================
def main():
    rows = []
    all_energy = {}
    mc_samples = {}

    for key, sc in SCENARIOS.items():
        # 名义参数仿真（出图）
        res = simulate_one(sc, t_max_h=24.0, n_eval=2200,
                           params={"R0_ref": 0.08, "Ea_R": 1800.0})

        # 每场景出图 + 能量汇总
        E = plot_bundle(res, key)
        all_energy[key] = E

        row = {
            "ScenarioKey": key,
            "ScenarioName": sc.name,
            "SOC0": sc.SOC0,
            "T_amb(C)": sc.T_amb,
            "SOH": sc.SOH,
            "t_empty(h)": res["t_empty"],
            "ShutdownReason": res["reason"],
            "SOC_end(%)": float(100 * res["SOC"][-1]),
            "Vterm_end(V)": float(res["Vterm"][-1]),
            "T_end(C)": float(res["T"][-1]),
        }

        # MC + 进度条
        if RUN_MC:
            mc = monte_carlo_t_empty(sc, N=MC_N, seed=7)
            row.update({"MC_p10(h)": mc["p10"], "MC_p50(h)": mc["p50"], "MC_p90(h)": mc["p90"]})
            mc_samples[key] = mc["samples"]

        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv("Q2_summary_table.csv", index=False, encoding="utf-8-sig")

    # 全局对比图
    plot_energy_3d(all_energy)
    plot_tempty_bar(df)
    plot_reason_share(df)
    plot_s5_tornado(df)
    if RUN_MC and len(mc_samples) > 0:
        plot_mc_boxplot(mc_samples)

    print("Done.")
    print("Saved CSV: Q2_summary_table.csv")

if __name__ == "__main__":
    main()
