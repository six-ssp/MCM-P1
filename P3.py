"""
Q3 模型检验：局部灵敏度（Tornado） + Sobol 全局灵敏度
严格按《第三问思路.docx》：
1) 局部灵敏度：弹性系数 S_i（归一化偏导）→ Tornado 图
2) 全局灵敏度：Sobol 一阶 / 总阶 → 条形图
3) 综合对比表：Local排名 + Sobol(Si, STi) → CSV

- 模型承接Q1/Q2：SOC ODE + ECM + 热耦合 + 终止条件
"""

import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt

from dataclasses import dataclass
from typing import Dict, List, Tuple

from scipy.integrate import solve_ivp
from scipy.stats import qmc


try:
    from tqdm import tqdm
except Exception:
    tqdm = None



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



SOC_min = 0.05
V_cutoff = 3.0
T_max = 45.0

C_th = 180.0
k_cool = 8.0

P_idle = 0.35
NET_POWER = {"4G": 1.1, "5G": 1.8, "WiFi": 0.9}

K_SCREEN_NOM = 2.6
K_CPU_NOM = 3.2

R0_REF_NOM = 0.08
EA_R_NOM = 1800.0

C_EFF_NOM = 17.0  # Wh


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

def scenario_S1() -> Scenario:

    return Scenario(
        name="S1 (Daily Baseline)",
        T_amb=25.0,
        SOH=1.0,
        SOC0=0.80,
        B=0.50,
        screen_duty=0.40,
        U_cpu=0.30,
        network="4G",
        net_duty=0.25,
        gps_on=False
    )


def duty_wave_soft(t_h, duty, period=5.0/60.0, smooth=0.03):
    if duty <= 0.0:
        return 0.0
    if duty >= 1.0:
        return 1.0
    phase = (t_h % period) / period
    k = 1.0 / max(smooth, 1e-3)
    up = 0.5 * (1 + np.tanh(k * (phase - 0.0)))
    down = 0.5 * (1 - np.tanh(k * (phase - duty)))
    return float(np.clip(up * down, 0.0, 1.0))

def power_components(t_h, sc: Scenario, k_screen: float, k_cpu: float, net_scale: float):
    I_screen = duty_wave_soft(t_h, sc.screen_duty)
    I_net = duty_wave_soft(t_h, sc.net_duty)

    return {
        "Base": P_idle,
        "Screen": k_screen * sc.B * I_screen,
        "CPU": k_cpu * sc.U_cpu,
        "Network": net_scale * NET_POWER.get(sc.network, 1.1) * I_net,
        "GPS": 0.8 if sc.gps_on else 0.0,
    }



def Voc(SOC, T_c):
    soc = np.clip(SOC, 0.0, 1.0)
    v = 3.0 + 1.2*(1/(1+np.exp(-8*(soc-0.35))))
    v += 0.002*(T_c - 25.0)
    return float(v)

def R0(T_c, SOH, R0_ref, Ea_R):
    T_k = (T_c + 273.15)
    Tref_k = 298.15
    temp_factor = np.exp(Ea_R*(1/T_k - 1/Tref_k))
    return float(R0_ref * temp_factor / max(SOH, 0.2))

def rc_params(T_c):
    R1, C1 = 0.015, 2000.0
    R2, C2 = 0.035, 6000.0
    factor = 1.0 + 0.01*(25.0 - T_c)
    factor = np.clip(factor, 0.8, 1.6)
    return float(R1*factor), float(C1), float(R2*factor), float(C2)



def rhs(t_h, y, sc: Scenario, C_eff, R0_ref, Ea_R, k_screen, k_cpu, net_scale):

    SOC, V1, V2, T_c, SOH = y

    comp = power_components(t_h, sc, k_screen, k_cpu, net_scale)
    P = float(sum(comp.values()))

    voc = Voc(SOC, T_c)
    r0 = R0(T_c, SOH, R0_ref=R0_ref, Ea_R=Ea_R)
    R1, C1, R2, C2 = rc_params(T_c)

    v_est = max(voc - V1 - V2, 3.0)
    I = P / v_est

    dSOC_dt = -P / C_eff

    tau1_h = (R1*C1) / 3600.0
    tau2_h = (R2*C2) / 3600.0
    dV1_dt = -(V1 / tau1_h) + (I / C1) * 3600.0
    dV2_dt = -(V2 / tau2_h) + (I / C2) * 3600.0

    P_heat = (I**2)*r0 + I*V1 + I*V2
    dT_dt = (P_heat / C_th)*3600.0 - (k_cool*(T_c - sc.T_amb)/C_th)*3600.0

    dSOH_dt = 0.0  
    return [dSOC_dt, dV1_dt, dV2_dt, dT_dt, dSOH_dt]

def make_events(sc: Scenario, C_eff, R0_ref, Ea_R, k_screen, k_cpu, net_scale):
    def ev_soc(t_h, y): return y[0] - SOC_min
    ev_soc.terminal = True; ev_soc.direction = -1

    def ev_vcut(t_h, y):
        SOC, V1, V2, T_c, SOH = y
        comp = power_components(t_h, sc, k_screen, k_cpu, net_scale)
        P = float(sum(comp.values()))
        voc = Voc(SOC, T_c)
        r0 = R0(T_c, SOH, R0_ref=R0_ref, Ea_R=Ea_R)
        v_est = max(voc - V1 - V2, 3.0)
        I = P / v_est
        Vterm = voc - I*r0 - V1 - V2
        return Vterm - V_cutoff
    ev_vcut.terminal = True; ev_vcut.direction = -1

    def ev_th(t_h, y): return T_max - y[3]
    ev_th.terminal = True; ev_th.direction = -1

    return [ev_soc, ev_vcut, ev_th]

def simulate_t_empty(sc: Scenario,
                     C_eff=C_EFF_NOM,
                     R0_ref=R0_REF_NOM,
                     Ea_R=EA_R_NOM,
                     k_screen=K_SCREEN_NOM,
                     k_cpu=K_CPU_NOM,
                     net_scale=1.0,
                     t_max_h=24.0,
                     n_eval=450):
    y0 = [sc.SOC0, 0.0, 0.0, sc.T_amb, sc.SOH]
    t_eval = np.linspace(0, t_max_h, n_eval)

    sol = solve_ivp(
        fun=lambda tt, yy: rhs(tt, yy, sc, C_eff, R0_ref, Ea_R, k_screen, k_cpu, net_scale),
        t_span=(0, t_max_h),
        y0=y0,
        t_eval=t_eval,
        method="BDF",
        rtol=1e-6,
        atol=1e-8,
        events=make_events(sc, C_eff, R0_ref, Ea_R, k_screen, k_cpu, net_scale),
    )
    return float(sol.t[-1])



PARAMS = [
    ("C_eff (Wh)",       (0.90*C_EFF_NOM,   1.10*C_EFF_NOM)),     # 等效容量
    ("R0_ref (Ohm)",     (0.90*R0_REF_NOM,  1.10*R0_REF_NOM)),    # 内阻参考
    ("Ea_R",             (0.90*EA_R_NOM,    1.10*EA_R_NOM)),      # Arrhenius 指数
    ("SOH",              (0.90*1.0,         1.10*1.0)),           # 健康度（基准S1=1.0）
    ("T_amb (C)",        (25.0-5.0,         25.0+5.0)),           # 外部环境温度（示例±5°C）
    ("B (brightness)",   (0.90*0.50,        1.10*0.50)),          # 屏幕亮度
    ("U_cpu",            (0.90*0.30,        1.10*0.30)),          # CPU 负载
    ("Net scale",        (0.90*1.0,         1.10*1.0)),           # 网络功耗缩放（模型简化/标定误差）
]

def eval_model_from_x(base_sc: Scenario, x: np.ndarray) -> float:

    sc = Scenario(**base_sc.__dict__)

    C_eff = float(x[0])
    R0_ref = float(x[1])
    Ea_R = float(x[2])
    sc.SOH = float(x[3])
    sc.T_amb = float(x[4])
    sc.B = float(x[5])
    sc.U_cpu = float(x[6])
    net_scale = float(x[7])

    return simulate_t_empty(sc,
                            C_eff=C_eff,
                            R0_ref=R0_ref,
                            Ea_R=Ea_R,
                            k_screen=K_SCREEN_NOM,
                            k_cpu=K_CPU_NOM,
                            net_scale=net_scale)


def local_sensitivity_tornado(base_sc: Scenario,
                              delta=0.10,
                              save_png="Q3-1_local_tornado.png") -> pd.DataFrame:

    x0 = np.array([(lo + hi) / 2 for _, (lo, hi) in PARAMS], dtype=float)
    t0 = eval_model_from_x(base_sc, x0)

    records = []
    for i, (name, (lo, hi)) in enumerate(PARAMS):
        xi0 = x0[i]
        x_plus = x0.copy();  x_minus = x0.copy()
        x_plus[i] = float(np.clip(xi0*(1+delta), lo, hi))
        x_minus[i] = float(np.clip(xi0*(1-delta), lo, hi))

        t_plus = eval_model_from_x(base_sc, x_plus)
        t_minus = eval_model_from_x(base_sc, x_minus)

        Si = (t_plus - t_minus) / (2*delta*t0)  
        records.append((name, Si, t_plus, t_minus))

    df = pd.DataFrame(records, columns=["Parameter", "Local_Si", "t_plus", "t_minus"])
    df["abs_Si"] = df["Local_Si"].abs()
    df = df.sort_values("abs_Si", ascending=True).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(10.6, 5.2))
    y = np.arange(len(df))
    ax.barh(y, df["Local_Si"].values)
    ax.axvline(0, color="#444444", linewidth=1.2)

    ax.set_yticks(y)
    ax.set_yticklabels(df["Parameter"].tolist())
    ax.set_xlabel("Local Sensitivity (Elasticity)  $S_i$")
    ax.set_title("Q3-1 Local Sensitivity (Tornado): $t_{empty}$ vs. Parameters")
    beautify_axes(ax)

    fig.savefig(save_png)
    plt.close(fig)
    return df

def sobol_global_indices(base_sc: Scenario,
                         N=256,
                         seed=7,
                         save_png="Q3-2_sobol_indices.png") -> pd.DataFrame:
  
    d = len(PARAMS)
    bounds = np.array([b for _, b in PARAMS], dtype=float)  
    sampler = qmc.Sobol(d=d, scramble=True, seed=seed)
    U = sampler.random(n=2*N)  # [0,1)
    X = qmc.scale(U, bounds[:, 0], bounds[:, 1])  

    A = X[:N, :]
    B = X[N:, :]

    def batch_eval(mat: np.ndarray, desc: str):
        out = np.zeros((mat.shape[0],), dtype=float)
        it = range(mat.shape[0])
        if tqdm is not None:
            it = tqdm(it, desc=desc, ncols=90)
        for k in it:
            out[k] = eval_model_from_x(base_sc, mat[k, :])
        return out

    YA = batch_eval(A, "Sobol YA")
    YB = batch_eval(B, "Sobol YB")

    VarY = np.var(np.concatenate([YA, YB]), ddof=1)
    VarY = float(VarY if VarY > 1e-12 else 1e-12)

    Si = np.zeros((d,), dtype=float)
    STi = np.zeros((d,), dtype=float)

    for i in range(d):
        ABi = A.copy()
        ABi[:, i] = B[:, i]
        YABi = batch_eval(ABi, f"Sobol YAB[{i+1}/{d}]")

        Si[i] = np.mean(YB * (YABi - YA)) / VarY
        STi[i] = 0.5 * np.mean((YA - YABi)**2) / VarY

    df = pd.DataFrame({
        "Parameter": [p[0] for p in PARAMS],
        "Sobol_Si": Si,
        "Sobol_STi": STi
    })

    fig, ax = plt.subplots(figsize=(12.0, 5.0))
    x = np.arange(d)
    w = 0.38
    ax.bar(x - w/2, df["Sobol_Si"].values, width=w, label="First-order $S_i$")
    ax.bar(x + w/2, df["Sobol_STi"].values, width=w, label="Total-order $S_{T_i}$")

    ax.set_xticks(x)
    ax.set_xticklabels(df["Parameter"].tolist(), rotation=20, ha="right")
    ax.set_ylabel("Sobol Index (Variance Contribution)")
    ax.set_title("Q3-2 Global Sensitivity (Sobol): First-order vs Total-order")
    beautify_axes(ax)

    fig.subplots_adjust(right=0.78)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)

    fig.savefig(save_png)
    plt.close(fig)
    return df


def make_comparison_table(local_df: pd.DataFrame,
                          sobol_df: pd.DataFrame,
                          save_csv="Q3_comparison_table.csv") -> pd.DataFrame:
    tmp = local_df.copy()
    tmp = tmp.sort_values("abs_Si", ascending=False).reset_index(drop=True)
    tmp["LocalRank"] = np.arange(1, len(tmp) + 1)
    tmp = tmp[["Parameter", "Local_Si", "LocalRank"]]

    df = pd.merge(tmp, sobol_df, on="Parameter", how="left")
    df.to_csv(save_csv, index=False, encoding="utf-8-sig")
    return df


def main():
    base = scenario_S1()

    # 3.1 局部灵敏度（Tornado）
    local_df = local_sensitivity_tornado(base, delta=0.10,
                                         save_png="Q3-1_local_tornado.png")

    # 3.2 全局灵敏度（Sobol）
    sobol_df = sobol_global_indices(base, N=256, seed=7,
                                    save_png="Q3-2_sobol_indices.png")

    # 3.3 综合对比表
    comp = make_comparison_table(local_df, sobol_df,
                                 save_csv="Q3_comparison_table.csv")

    print("Done.")
    print("Saved: Q3-1_local_tornado.png")
    print("Saved: Q3-2_sobol_indices.png")
    print("Saved: Q3_comparison_table.csv")
    print(comp)

if __name__ == "__main__":
    main()
