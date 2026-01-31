#pip install -U plotly
#pip install -U kaleido

import numpy as np
import plotly.graph_objects as go

FONT_FAMILY = "Times New Roman"   
T_END = 8.0
N = 2501
SOC0 = 1.0

C_EFF_WH = 17.0
P_BASE_W = 0.35
ALPHA_SCREEN_W = 2.2

t = np.linspace(0, T_END, N)


def piecewise_array(t_arr, segments):
    """segments: list of (t0, t1, value)"""
    y = np.zeros_like(t_arr, dtype=float)
    for t0, t1, v in segments:
        y[(t_arr >= t0) & (t_arr < t1)] = v
    return y

def scenario_segments(name):
    if name == "Light":
        I_seg = [(0.0,1.5,0),(1.5,1.7,1),(1.7,4.0,0),(4.0,4.2,1),(4.2,8.0,0)]
        B_seg = [(0.0,8.0,0.35)]
    elif name == "Moderate":
        I_seg = [(0.0,0.5,0),(0.5,1.2,1),(1.2,2.0,0),
                 (2.0,3.0,1),(3.0,3.7,0),(3.7,5.3,1),(5.3,8.0,0)]
        B_seg = [(0.0,2.5,0.45),(2.5,5.5,0.60),(5.5,8.0,0.40)]
    elif name == "Heavy":
        I_seg = [(0.0,0.3,0),(0.3,7.5,1),(7.5,8.0,0)]
        B_seg = [(0.0,8.0,0.85)]
    else:
        raise ValueError("Unknown scenario. Use: Light / Moderate / Heavy")
    return I_seg, B_seg


def integrate_soc(t_arr, P_arr, C_eff_wh, soc0=1.0):
    soc = np.zeros_like(P_arr, dtype=float)
    soc[0] = soc0
    for k in range(1， len(t_arr)):
        dt = t_arr[k] - t_arr[k-1]
        soc[k] = soc[k-1] - (P_arr[k-1] / C_eff_wh) * dt
    return soc

SCENARIOS = ["Light", "Moderate", "Heavy"]
data = {}

for s in SCENARIOS:
    I_seg, B_seg = scenario_segments(s)
    I = piecewise_array(t, I_seg)
    B = piecewise_array(t, B_seg)
    P = P_BASE_W + ALPHA_SCREEN_W * B * I
    SOC = integrate_soc(t, P, C_EFF_WH, soc0=SOC0)
    data[s] = {"I": I, "B": B, "P": P, "SOC": SOC}

def apply_base_layout(fig, title, x_title, y_title, y_range=None, extra_right_margin=0):
    """统一设置：轴线、无网格、外置legend、Times字体、足够右边距"""
    right_margin = 280 + extra_right_margin  # 给 legend + 长文本留空间
    fig.update_layout(
        title=title,
        template="plotly_white",
        font=dict(family=FONT_FAMILY, size=16),
        margin=dict(l=85, r=right_margin, t=85, b=70),
        legend=dict(
            x=1.02, y=0.5,
            xanchor="left", yanchor="middle",
            bgcolor="rgba(0,0,0,0)",
            bordercolor="rgba(0,0,0,0)",
            font=dict(size=14),
        ),
    )

    fig.update_xaxes(
        title=x_title,
        showline=True, linewidth=2, linecolor="black",
        showgrid=False, zeroline=False,
        ticks="outside", tickwidth=1.5, ticklen=6
    )
    fig.update_yaxes(
        title=y_title,
        showline=True, linewidth=2, linecolor="black",
        showgrid=False, zeroline=False,
        ticks="outside", tickwidth=1.5, ticklen=6
    )
    if y_range is not None:
        fig.update_yaxes(range=y_range)

def add_axis_arrows(fig):
    # x-axis arrow (→)
    fig.add_annotation(
        xref="paper", yref="paper",
        x=1.02, y=0.0,
        ax=1.0, ay=0.0,
        showarrow=True,
        arrowhead=3,
        arrowsize=1.2,
        arrowwidth=2,
        arrowcolor="black",
        text=""
    )
    # y-axis arrow (↑)
    fig.add_annotation(
        xref="paper", yref="paper",
        x=0.0, y=1.02,
        ax=0.0, ay=1.0,
        showarrow=True,
        arrowhead=3,
        arrowsize=1.2,
        arrowwidth=2,
        arrowcolor="black",
        text=""
    )


SOC_COLORS = {"Light": "#1f77b4", "Moderate": "#ff7f0e", "Heavy": "#2ca02c"}
PWR_COLORS = {"Light": "#9467bd", "Moderate": "#d62728", "Heavy": "#8c564b"}
INP_COLORS = {
    "Brightness": "#17becf",
    "Gate": "#7f7f7f",
    "ScreenTerm": "#e377c2"
}


fig1 = go.Figure()
for s in SCENARIOS:
    fig1.add_trace(go.Scatter(
        x=t, y=100 * data[s]["SOC"],
        mode="lines",
        name=f"SOC - {s}",
        line=dict(width=4, color=SOC_COLORS[s]),
        hovertemplate="Scenario: "+s+"<br>t=%{x:.2f} h<br>SOC=%{y:.2f}%<extra></extra>"
    ))

apply_base_layout(
    fig1,
    title="SOC Trajectory (Continuous-Time ODE)",
    x_title="Time t (hours)",
    y_title="SOC (%)",
    y_range=[0, 105]
)
add_axis_arrows(fig1)

fig1.show()
fig1.write_html("plotly_fig1_SOC.html", include_plotlyjs="cdn")

fig2 = go.Figure()
for s in SCENARIOS:
    fig2.add_trace(go.Scatter(
        x=t, y=data[s]["P"],
        mode="lines",
        name=f"Power - {s}",
        line=dict(width=4, color=PWR_COLORS[s]),
        hovertemplate="Scenario: "+s+"<br>t=%{x:.2f} h<br>P=%{y:.3f} W<extra></extra>"
    ))

fig2.add_trace(go.Scatter(
    x=t, y=[P_BASE_W]*len(t),
    mode="lines",
    name="Baseline P_base",
    line=dict(width=2, color="black", dash="dash"),
    hoverinfo="skip"
))

apply_base_layout(
    fig2,
    title="Instantaneous Power Consumption P(t)",
    x_title="Time t (hours)",
    y_title="Power P(t) (W)"
)
add_axis_arrows(fig2)

fig2.show()
fig2.write_html("plotly_fig2_Power.html", include_plotlyjs="cdn")

m = data["Moderate"]
screen_term = ALPHA_SCREEN_W * m["B"] * m["I"]

fig3 = go.Figure()

fig3.add_trace(go.Scatter(
    x=t, y=m["B"],
    mode="lines",
    name="Brightness B(t) (0–1)",
    line=dict(width=4, color=INP_COLORS["Brightness"]),
    hovertemplate="t=%{x:.2f} h<br>B=%{y:.2f}<extra></extra>"
))

fig3.add_trace(go.Scatter(
    x=t, y=m["I"],
    mode="lines",
    name="Screen-on Gate I_screen(t) (0/1)",
    line=dict(width=4, color=INP_COLORS["Gate"], shape="hv"),
    hovertemplate="t=%{x:.2f} h<br>I=%{y:.0f}<extra></extra>"
))

fig3.add_trace(go.Scatter(
    x=t, y=screen_term,
    mode="lines",
    name="Screen Term α·B·I (W)",
    line=dict(width=4, color=INP_COLORS["ScreenTerm"]),
    hovertemplate="t=%{x:.2f} h<br>αBI=%{y:.3f} W<extra></extra>"
))

apply_base_layout(
    fig3,
    title="Inputs / Gating Functions (Moderate Scenario)",
    x_title="Time t (hours)",
    y_title="Value",
    extra_right_margin=40  
)
add_axis_arrows(fig3)

fig3.show()
fig3.write_html("plotly_fig3_Inputs.html", include_plotlyjs="cdn")

print("Done. Exported HTML: plotly_fig1_SOC.html, plotly_fig2_Power.html, plotly_fig3_Inputs.html")
