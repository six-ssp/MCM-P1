#P1 连续时间 SOC 最简动力学模型
#图左上 模拟的是系统状态 SOC(t) 的演化
#图右上 模拟的是功耗函数 P(t) 的时间变化
#图下 具体三种不同程度的使用情况屏幕亮起与时间关系
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def setup_plot_style():
    mpl.rcParams["font.family"] = "Times New Roman"
    mpl.rcParams["axes.unicode_minus"] = False  
    mpl.rcParams["figure.dpi"] = 140
    mpl.rcParams["savefig.dpi"] = 220

setup_plot_style()

#  1) 参数，随便的常见数值

# 电池有效容量（Wh），取合理量级
C_eff_Wh = 17.0
# 基线功耗（W）：待机/后台
P_base_W = 0.35
# 屏幕亮度功耗系数（W）：亮度=1、亮屏时增加的功耗量级
alpha_screen_W = 2.2

#  2) 构造使用情景（输入函数 B(t), I_screen(t)）

def piecewise_schedule(t, segments):
    """
    简单分段函数：
    segments: list of (t0, t1, value), t 单位：小时
    """
    for (t0, t1, val) in segments:
        if t0 <= t < t1:
            return val
    return segments[-1][2] if segments else 0.0


def make_scenario(name="Light"):
    """
    返回一个情景字典，包含：
    - I_screen(t): 屏幕开关门控（0/1）
    - B(t): 亮度（0~1）
    """
    if name == "Light":
        # 轻度：主要待机，偶尔亮屏
        screen_on = [
            (0.0, 1.5, 0),
            (1.5, 1.7, 1),
            (1.7, 4.0, 0),
            (4.0, 4.2, 1),
            (4.2, 8.0, 0),
        ]
        brightness = [(0.0, 8.0, 0.35)]
    elif name == "Moderate":
        # 中度：亮屏更频繁，亮度有变化
        screen_on = [
            (0.0, 0.5, 0),
            (0.5, 1.2, 1),
            (1.2, 2.0, 0),
            (2.0, 3.0, 1),
            (3.0, 3.7, 0),
            (3.7, 5.3, 1),
            (5.3, 8.0, 0),
        ]
        brightness = [
            (0.0, 2.5, 0.45),
            (2.5, 5.5, 0.60),
            (5.5, 8.0, 0.40),
        ]
    elif name == "Heavy":
        # 重度：长时间高亮度亮屏
        screen_on = [
            (0.0, 0.3, 0),
            (0.3, 7.5, 1),
            (7.5, 8.0, 0),
        ]
        brightness = [(0.0, 8.0, 0.85)]
    else:
        raise ValueError("Unknown scenario name. Use: Light / Moderate / Heavy")

    def I_screen(t):
        return piecewise_schedule(t, screen_on)

    def B(t):
        return piecewise_schedule(t, brightness)

    return {"name": name, "I_screen": I_screen, "B": B}


def P_total(t, scenario):
    """
    最简功耗结构：
        P(t) = P_base + alpha * B(t) * I_screen(t)
    """
    return P_base_W + alpha_screen_W * scenario["B"](t) * scenario["I_screen"](t)

#  3) 连续时间 SOC 微分方程与数值求解 

def dSOC_dt(t, soc, scenario):
    """
    ODE: dSOC/dt = -P(t) / C_eff
    单位核对：
        P (W) = Wh / h
        C_eff (Wh)
        => P/C_eff 是 1/h
    """
    return -P_total(t, scenario) / C_eff_Wh


def simulate_soc(scenario, soc0=1.0, t_end=8.0, n=2001):
    """
    求解 SOC(t)，并返回 t, SOC, P(t), I_screen(t), B(t)
    """
    t_eval = np.linspace(0, t_end, n)
    sol = solve_ivp(
        fun=lambda t, y: dSOC_dt(t, y, scenario),
        t_span=(0, t_end),
        y0=[soc0],
        t_eval=t_eval,
        rtol=1e-7,
        atol=1e-9
    )

    t = sol.t
    soc = sol.y[0]

    P = np.array([P_total(tt, scenario) for tt in t])
    I = np.array([scenario["I_screen"](tt) for tt in t])
    B = np.array([scenario["B"](tt) for tt in t])

    # 估计 time-to-empty（若未耗尽则 None）
    idx = np.where(soc <= 0)[0]
    t_empty = float(t[idx[0]]) if len(idx) > 0 else None

    return t, soc, P, I, B, t_empty

# 4) 绘图

def place_legend_outside_right(ax, **kwargs):
    """
    将 legend 放到图外右侧，避免盖住曲线
    """
    return ax.legend(
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        borderaxespad=0.0,
        frameon=True,
        **kwargs
    )


def plot_results(results, t_end=8.0, soc0=1.0):
    """
    results: list of dict {name, t, soc, P, I, B, t_empty}
    """
    fig = plt.figure(figsize=(12.8, 8.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 0.95], hspace=0.28, wspace=0.6)
    fig.subplots_adjust(right=0.82)  

    ax1 = fig.add_subplot(gs[0, 0])  # SOC
    ax2 = fig.add_subplot(gs[0, 1])  # P(t)
    ax3 = fig.add_subplot(gs[1, :])  # inputs (Moderate)

    # --- (1) SOC vs time ---
    for r in results:
        ax1.plot(r["t"], 100 * r["soc"], linewidth=2.2, label=r["name"])

        if r["t_empty"] is not None:
            ax1.axvline(r["t_empty"], linestyle="--", linewidth=1.2)

    ax1.set_title("SOC Trajectory (Continuous-Time ODE Model)")
    ax1.set_xlabel("Time t (hours)")
    ax1.set_ylabel("SOC (%)")
    ax1.set_xlim(0, t_end)
    ax1.set_ylim(0, 105)
    ax1.grid(True, alpha=0.25)
    place_legend_outside_right(ax1)

    # --- (2) Power consumption P(t) ---
    for r in results:
        ax2.plot(r["t"], r["P"], linewidth=2.0, label=r["name"])

    ax2.set_title("Instantaneous Power P(t) (Minimal Structure)")
    ax2.set_xlabel("Time t (hours)")
    ax2.set_ylabel("Power P(t) (W)")
    ax2.set_xlim(0, t_end)
    ax2.grid(True, alpha=0.25)
    place_legend_outside_right(ax2)

    # --- (3) Inputs example: choose Moderate if exists ---
    mid = next((r for r in results if r["name"] == "Moderate"), results[0])
    P_screen = alpha_screen_W * mid["B"] * mid["I"]

    ax3.plot(mid["t"], mid["B"], linewidth=2.2, label="Brightness B(t) (0-1)")
    ax3.step(mid["t"], mid["I"], where="post", linewidth=2.2, label="Screen-on Gate I_screen(t) (0/1)")
    ax3.plot(mid["t"], P_screen, linewidth=2.0, label="Screen Term α·B·I (W)")

    ax3.set_title("Input / Gating Functions Example (Moderate Scenario)")
    ax3.set_xlabel("Time t (hours)")
    ax3.set_ylabel("Value")
    ax3.set_xlim(0, t_end)
    ax3.grid(True, alpha=0.25)
    place_legend_outside_right(ax3)

    # --- super title ---
    fig.suptitle(
        f"Q1 Minimal Continuous-Time Model: dSOC/dt = -P(t)/C_eff   "
        f"(C_eff={C_eff_Wh:.1f} Wh, SOC0={soc0*100:.0f}%)",
        y=0.98,
        fontsize=12
    )

    plt.show()
    print("\n=== Summary (8 hours) ===")
    for r in results:
        soc_end = 100 * r["soc"][-1]
        te = r["t_empty"]
        te_str = f"{te:.2f} h" if te is not None else "Not depleted"
        print(f"{r['name']:<9s} | SOC(8h) = {soc_end:6.2f}% | Time-to-empty = {te_str}")


#5) 主程序入口

if __name__ == "__main__":
    t_end = 8.0
    soc0 = 1.0  # 初始 SOC，可改 0.8 / 0.5 做对比

    scenarios = [make_scenario("Light"), make_scenario("Moderate"), make_scenario("Heavy")]

    results = []
    for sc in scenarios:
        t, soc, P, I, B, t_empty = simulate_soc(sc, soc0=soc0, t_end=t_end)
        results.append({
            "name": sc["name"],
            "t": t,
            "soc": soc,
            "P": P,
            "I": I,
            "B": B,
            "t_empty": t_empty
        })

    plot_results(results, t_end=t_end, soc0=soc0)
