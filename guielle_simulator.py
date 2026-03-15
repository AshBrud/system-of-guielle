"""
Simulateur du Système de Guielle
=================================
Auteur : Herolde GUIELLE
Implémentation numérique : RK4 + formules exactes du document

Installation :
    pip install streamlit numpy scipy matplotlib plotly

Lancement :
    streamlit run guielle_simulator.py
"""

import numpy as np
import streamlit as st
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import time

# ─── Configuration de la page ─────────────────────────────────────────────────
st.set_page_config(
    page_title="Système de Guielle — Simulateur",
    page_icon="⚙️",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown("""
<style>
    .metric-box {
        background: #f0f4f8;
        border-radius: 8px;
        padding: 12px 16px;
        border-left: 4px solid #2E75B6;
        margin: 4px 0;
    }
    .metric-label { font-size: 12px; color: #666; margin: 0; }
    .metric-value { font-size: 22px; font-weight: 600; color: #1F4E79; margin: 0; }
    .cond-ok   { background:#e6f4ea; color:#1e7e34; border-radius:6px; padding:4px 10px; font-weight:600; }
    .cond-warn { background:#fff3cd; color:#856404; border-radius:6px; padding:4px 10px; font-weight:600; }
    .cond-fail { background:#fce8e6; color:#c62828; border-radius:6px; padding:4px 10px; font-weight:600; }
    h1 { color: #1F4E79; }
    .stSlider > div > div > div { color: #2E75B6; }
</style>
""", unsafe_allow_html=True)


# ─── Formules exactes du document ─────────────────────────────────────────────

def spring_data(theta: float, i: int, n: int, r: float, l: float):
    """
    Calcule les données géométriques du ressort i.
    Formules fidèles au document de H. GUIELLE.

    Retourne :
      xi_r  : élongation adimensionnelle xi/r
      alpha : angle alpha_i (rad)
      angle : angle (e_rho, e_xi) = -theta + beta_i + alpha_i
      Ax,Ay : coordonnées de A_i (fixes dans le labo)
    """
    beta = i * 2 * np.pi / n                      # angle fixe de A_i
    tau  = beta + theta                            # τ_i = β_i + θ  (doc)
    eps  = 1.0 + l / r                            # ε = 1 + l₀/r
    D    = np.sqrt(1 + eps**2 + 2*eps*np.cos(tau)) # (x_i+l)/r
    xi_r = D - l / r                               # x_i/r  (doc)
    sin_a = np.clip(np.sin(tau) / D, -1, 1)
    alpha = np.arcsin(sin_a)                       # α_i  (doc)
    angle = -theta + beta + alpha                  # (ê_ρ, ê_{xi})  (doc)
    # A_i fixe : O + r·ε·(cos β_i, sin β_i)
    Ax = r * eps * np.cos(beta)
    Ay = r * eps * np.sin(beta)
    return xi_r, alpha, angle, Ax, Ay


def compute_g(theta: float, n: int, r: float, l: float) -> float:
    """Moment adimensionnel g(t) — formule exacte du document."""
    g = 0.0
    for i in range(1, n + 1):
        xi_r, _, angle, _, _ = spring_data(theta, i, n, r, l)
        if xi_r > 0:
            g += xi_r * np.sin(angle)
    return g


def compute_moment(theta: float, n: int, r: float, k: float, l: float) -> float:
    """Moment résultant M (N·m) = −r²·k·g(t)."""
    return -r**2 * k * compute_g(theta, n, r, l)


def compute_energies(theta: float, omega: float, n: int, r: float,
                     k: float, l: float, I: float):
    """Énergies cinétique et élastique totale."""
    Ec = 0.5 * I * omega**2
    Ee = 0.0
    for i in range(1, n + 1):
        xi_r, _, _, _, _ = spring_data(theta, i, n, r, l)
        if xi_r > 0:
            xi = xi_r * r
            Ee += 0.5 * k * xi**2
    return Ec, Ee


# ─── Intégration RK4 ──────────────────────────────────────────────────────────

def deriv(theta: float, omega: float, n: int, r: float, k: float,
          l: float, I: float, b: float, Cs: float):
    M = compute_moment(theta, n, r, k, l)
    friction = b * omega + Cs * np.sign(omega) if omega != 0 else b * omega
    d_omega = (M - friction) / I
    return omega, d_omega


def rk4_step(theta: float, omega: float, dt: float,
             n: int, r: float, k: float, l: float,
             I: float, b: float, Cs: float):
    """Un pas RK4."""
    k1_th, k1_om = deriv(theta, omega, n, r, k, l, I, b, Cs)
    k2_th, k2_om = deriv(theta + k1_th*dt/2, omega + k1_om*dt/2, n, r, k, l, I, b, Cs)
    k3_th, k3_om = deriv(theta + k2_th*dt/2, omega + k2_om*dt/2, n, r, k, l, I, b, Cs)
    k4_th, k4_om = deriv(theta + k3_th*dt,   omega + k3_om*dt,   n, r, k, l, I, b, Cs)
    new_theta = theta + dt * (k1_th + 2*k2_th + 2*k3_th + k4_th) / 6
    new_omega = omega + dt * (k1_om + 2*k2_om + 2*k3_om + k4_om) / 6
    return new_theta, new_omega


def simulate(params: dict, t_end: float = 20.0, dt: float = 0.012):
    """
    Simule le système sur [0, t_end] avec un pas RK4 dt.
    Retourne un dictionnaire de séries temporelles.
    """
    n, r, k, l = params['n'], params['r'], params['k'], params['l']
    m, alpha   = params['m'], params['alpha']
    b, Cs      = params['b'], params['Cs']
    omega0     = params['omega0']
    I = alpha * m * r**2

    steps = int(t_end / dt)
    t_arr     = np.zeros(steps)
    theta_arr = np.zeros(steps)
    omega_arr = np.zeros(steps)
    g_arr     = np.zeros(steps)
    Ec_arr    = np.zeros(steps)
    Ee_arr    = np.zeros(steps)
    M_arr     = np.zeros(steps)

    theta, omega = 0.0, omega0
    for i in range(steps):
        t = i * dt
        g  = compute_g(theta, n, r, l)
        M  = -r**2 * k * g
        Ec, Ee = compute_energies(theta, omega, n, r, k, l, I)
        t_arr[i]     = t
        theta_arr[i] = theta
        omega_arr[i] = omega
        g_arr[i]     = g
        Ec_arr[i]    = Ec
        Ee_arr[i]    = Ee
        M_arr[i]     = M
        theta, omega = rk4_step(theta, omega, dt, n, r, k, l, I, b, Cs)

    return {
        't': t_arr, 'theta': theta_arr, 'omega': omega_arr,
        'g': g_arr, 'Ec': Ec_arr, 'Ee': Ee_arr, 'M': M_arr
    }


# ─── Triple condition ──────────────────────────────────────────────────────────

def check_triple_condition(g_arr: np.ndarray, tail: int = 500):
    """Évalue la triple condition sur la fin de la simulation."""
    gs = g_arr[-tail:] if len(g_arr) > tail else g_arr
    min_abs = np.min(np.abs(gs))
    max_abs = np.max(np.abs(gs))
    rel_var = (max_abs - min_abs) / max_abs if max_abs > 1e-9 else 1.0
    all_pos = np.all(gs >= 0)
    all_neg = np.all(gs <= 0)

    c1 = ("✓ |g| > 0", "ok")   if min_abs > 1e-4 else ("✗ |g| ≈ 0", "fail")
    c2_label = f"✓ g constant (var. {rel_var*100:.1f}%)" if rel_var < 0.15 \
          else (f"~ g quasi-const. ({rel_var*100:.1f}%)" if rel_var < 0.4 \
          else  f"✗ g variable ({rel_var*100:.1f}%)")
    c2_status = "ok" if rel_var < 0.15 else ("warn" if rel_var < 0.4 else "fail")
    c3 = ("✓ Signe fixe", "ok") if (all_pos or all_neg) else ("✗ Signe change", "fail")

    return [c1, (c2_label, c2_status), c3]


# ─── Animation roue ───────────────────────────────────────────────────────────

def draw_wheel_plotly(theta: float, params: dict):
    """Dessine la roue, les ressorts et la force résultante avec Plotly."""
    n, r, k, l = params['n'], params['r'], params['k'], params['l']
    eps = 1.0 + l / r

    fig = go.Figure()
    fig.update_layout(
        width=380, height=380,
        margin=dict(l=10, r=10, t=30, b=10),
        showlegend=False,
        xaxis=dict(range=[-1.8, 1.8], scaleanchor="y", showgrid=False,
                   zeroline=False, showticklabels=False),
        yaxis=dict(range=[-1.8, 1.8], showgrid=False,
                   zeroline=False, showticklabels=False),
        plot_bgcolor="white",
        title=dict(text=f"θ = {theta % (2*np.pi):.2f} rad", font=dict(size=12)),
    )

    # Cercle roue
    ang = np.linspace(0, 2*np.pi, 200)
    fig.add_trace(go.Scatter(
        x=np.cos(ang), y=np.sin(ang),
        mode='lines', line=dict(color='#CCCCCC', width=1),
    ))

    # Point P
    Px, Py = np.cos(theta), np.sin(theta)

    # Rayon OP
    fig.add_trace(go.Scatter(
        x=[0, Px], y=[0, Py],
        mode='lines', line=dict(color='#AAAAAA', width=1, dash='dot'),
    ))

    # Ressorts et points Ai
    Frx, Fry = 0.0, 0.0
    for i in range(1, n + 1):
        xi_r, _, _, Ax, Ay = spring_data(theta, i, n, r, l)
        Ax_n, Ay_n = Ax / r, Ay / r  # normalisé
        tendu = xi_r > 0
        color = '#E8593C' if tendu else '#BBBBBB'
        dash  = 'solid'   if tendu else 'dot'
        width = 2.0       if tendu else 1.0

        fig.add_trace(go.Scatter(
            x=[Px, Ax_n], y=[Py, Ay_n],
            mode='lines', line=dict(color=color, width=width, dash=dash),
        ))
        # Point Ai fixe
        fig.add_trace(go.Scatter(
            x=[Ax_n], y=[Ay_n],
            mode='markers',
            marker=dict(color='#BA7517', size=8, symbol='square'),
        ))

        if tendu:
            dist = np.sqrt((Px - Ax_n)**2 + (Py - Ay_n)**2)
            if dist > 1e-6:
                nx = (Ax_n - Px) / dist
                ny = (Ay_n - Py) / dist
                Frx += k * xi_r * r * nx
                Fry += k * xi_r * r * ny

    # Flèche force résultante
    Fmag = np.sqrt(Frx**2 + Fry**2)
    if Fmag > 1e-6:
        scale = 0.5 / (k * l + 1)
        ex = Px + Frx * scale
        ey = Py + Fry * scale
        fig.add_annotation(
            x=ex, y=ey, ax=Px, ay=Py,
            xref="x", yref="y", axref="x", ayref="y",
            showarrow=True, arrowhead=2, arrowsize=1.5,
            arrowwidth=2.5, arrowcolor='#1D9E75',
        )

    # Centre O
    fig.add_trace(go.Scatter(
        x=[0], y=[0], mode='markers',
        marker=dict(color='#333333', size=8),
    ))
    # Point P
    fig.add_trace(go.Scatter(
        x=[Px], y=[Py], mode='markers',
        marker=dict(color='#3B8BD4', size=10),
    ))

    return fig


# ─── Interface Streamlit ──────────────────────────────────────────────────────

st.title("⚙️ Système de Guielle — Simulateur")
st.caption("Modélisation et simulation numérique (RK4) · Auteur : Herolde GUIELLE")

# ── Sidebar : paramètres ──
with st.sidebar:
    st.header("Paramètres du système")

    st.subheader("Géométrie")
    n      = st.slider("n — Nombre de ressorts", 2, 12, 4, 1)
    r      = st.slider("r — Rayon de la roue (m)", 0.05, 1.0, 0.30, 0.01)
    l0     = st.slider("l₀ — Longueur naturelle (m)", 0.05, 2.0, 0.40, 0.01)
    k      = st.slider("k — Raideur (N/m)", 1, 200, 30, 1)

    st.subheader("Inertie (libre)")
    m_val  = st.slider("m — Masse roue (kg)", 0.1, 30.0, 5.0, 0.1)
    alpha  = st.slider("α — Coefficient inertie I=α·m·r²", 0.10, 1.0, 0.50, 0.01)
    I_val  = alpha * m_val * r**2
    st.caption(f"→ I = {I_val:.4f} kg·m²")

    st.subheader("Moments résistants (cas réel)")
    b_val  = st.slider("b — Frottement visqueux", 0.0, 2.0, 0.0, 0.01)
    Cs_val = st.slider("Cs — Couple résistant (N·m)", 0.0, 10.0, 0.0, 0.05)

    st.subheader("Conditions initiales")
    omega0 = st.slider("ω₀ — Vitesse initiale (rad/s)", 0.1, 15.0, 2.0, 0.1)

    st.subheader("Simulation")
    t_end  = st.slider("Durée (s)", 5, 120, 30, 5)
    dt     = st.select_slider("Pas Δt (s)", [0.004, 0.008, 0.012, 0.02], value=0.012)

    run_btn = st.button("▶ Lancer la simulation", type="primary", use_container_width=True)

params = dict(n=n, r=r, k=k, l=l0, m=m_val, alpha=alpha, b=b_val, Cs=Cs_val, omega0=omega0)

# ── Résultats ──
if run_btn or 'results' not in st.session_state:
    with st.spinner("Simulation en cours (RK4)..."):
        st.session_state.results = simulate(params, t_end=t_end, dt=dt)
        st.session_state.params  = params.copy()

res    = st.session_state.results
p_used = st.session_state.get('params', params)

# ── Métriques instantanées ──
col1, col2, col3, col4, col5, col6 = st.columns(6)
with col1:
    st.metric("ω final (rad/s)", f"{res['omega'][-1]:.4f}")
with col2:
    st.metric("g final", f"{res['g'][-1]:.4f}")
with col3:
    st.metric("Ec final (J)", f"{res['Ec'][-1]:.3f}")
with col4:
    st.metric("Ee final (J)", f"{res['Ee'][-1]:.3f}")
with col5:
    st.metric("M final (N·m)", f"{res['M'][-1]:.4f}")
with col6:
    st.metric("θ final (rad)", f"{res['theta'][-1] % (2*np.pi):.3f}")

st.divider()

# ── Triple condition ──
conditions = check_triple_condition(res['g'])
st.subheader("Triple condition sur g(t) — derniers 500 pas")
cc1, cc2, cc3 = st.columns(3)
for col, (label, status) in zip([cc1, cc2, cc3], conditions):
    css_class = f"cond-{status}"
    col.markdown(f'<div class="{css_class}">{label}</div>', unsafe_allow_html=True)

st.divider()

# ── Graphiques ──
colA, colB = st.columns([1, 2])

with colA:
    st.subheader("Animation — position courante")
    theta_display = float(res['theta'][-1])
    fig_wheel = draw_wheel_plotly(theta_display, p_used)
    st.plotly_chart(fig_wheel, use_container_width=False)

with colB:
    st.subheader("Courbes temporelles")
    tab1, tab2, tab3, tab4 = st.tabs(["ω(t)", "g(t) — Triple condition", "Énergies", "Phase ω(θ)"])

    with tab1:
        fig_w = go.Figure()
        fig_w.add_trace(go.Scatter(x=res['t'], y=res['omega'],
                                   mode='lines', line=dict(color='#3B8BD4', width=1.5),
                                   name='ω(t)'))
        fig_w.add_hline(y=0, line_dash="dash", line_color="#AAAAAA")
        fig_w.update_layout(xaxis_title="t (s)", yaxis_title="ω (rad/s)",
                            height=320, margin=dict(l=50,r=20,t=20,b=40))
        st.plotly_chart(fig_w, use_container_width=True)

    with tab2:
        fig_g = go.Figure()
        fig_g.add_trace(go.Scatter(x=res['t'], y=res['g'],
                                   mode='lines', line=dict(color='#E8593C', width=1.5),
                                   name='g(t)'))
        fig_g.add_hline(y=0, line_dash="dash", line_color="#AAAAAA")
        # Bandes de référence
        g_mean = float(np.mean(res['g'][-500:]))
        if abs(g_mean) > 1e-6:
            fig_g.add_hrect(y0=g_mean*0.85, y1=g_mean*1.15,
                            fillcolor="#E8593C", opacity=0.08,
                            annotation_text="±15% autour de la moyenne")
        fig_g.update_layout(xaxis_title="t (s)", yaxis_title="g(t) adim.",
                            height=320, margin=dict(l=50,r=20,t=20,b=40))
        st.plotly_chart(fig_g, use_container_width=True)

    with tab3:
        fig_e = make_subplots(rows=1, cols=2,
                              subplot_titles=("Énergie cinétique Ec(t)", "Énergie élastique Ee(t)"))
        fig_e.add_trace(go.Scatter(x=res['t'], y=res['Ec'], mode='lines',
                                   line=dict(color='#1D9E75', width=1.5), name='Ec'), row=1, col=1)
        fig_e.add_trace(go.Scatter(x=res['t'], y=res['Ee'], mode='lines',
                                   line=dict(color='#BA7517', width=1.5), name='Ee'), row=1, col=2)
        fig_e.update_xaxes(title_text="t (s)")
        fig_e.update_yaxes(title_text="J", row=1, col=1)
        fig_e.update_layout(height=320, showlegend=False, margin=dict(l=50,r=20,t=40,b=40))
        st.plotly_chart(fig_e, use_container_width=True)

    with tab4:
        theta_mod = res['theta'] % (2 * np.pi)
        fig_ph = go.Figure()
        fig_ph.add_trace(go.Scatter(x=theta_mod, y=res['omega'],
                                    mode='lines',
                                    line=dict(color='#7F77DD', width=1, ),
                                    name='ω(θ mod 2π)'))
        fig_ph.update_layout(xaxis_title="θ mod 2π (rad)", yaxis_title="ω (rad/s)",
                             height=320, margin=dict(l=50,r=20,t=20,b=40))
        st.plotly_chart(fig_ph, use_container_width=True)

st.divider()

# ── Analyse paramétrique : comparaison de n ──
st.subheader("Analyse paramétrique — comparaison pour n = 2 à 8")
if st.button("Lancer l'analyse comparative de n", use_container_width=False):
    with st.spinner("Simulation pour n = 2, 3, 4, 5, 6, 7, 8..."):
        fig_comp = make_subplots(rows=1, cols=2,
                                 subplot_titles=("ω(t) pour différents n",
                                                 "g(t) pour différents n"))
        colors = ['#E8593C','#BA7517','#1D9E75','#3B8BD4',
                  '#7F77DD','#D4537E','#639922']
        for idx, ni in enumerate(range(2, 9)):
            p_n = params.copy(); p_n['n'] = ni
            r_n = simulate(p_n, t_end=min(t_end, 30), dt=dt)
            fig_comp.add_trace(go.Scatter(x=r_n['t'], y=r_n['omega'], mode='lines',
                               line=dict(color=colors[idx], width=1.5),
                               name=f"n={ni}"), row=1, col=1)
            fig_comp.add_trace(go.Scatter(x=r_n['t'], y=r_n['g'], mode='lines',
                               line=dict(color=colors[idx], width=1.5),
                               name=f"n={ni}", showlegend=False), row=1, col=2)
        fig_comp.update_xaxes(title_text="t (s)")
        fig_comp.update_yaxes(title_text="ω (rad/s)", row=1, col=1)
        fig_comp.update_yaxes(title_text="g(t)", row=1, col=2)
        fig_comp.update_layout(height=380, margin=dict(l=50,r=20,t=40,b=40))
        st.plotly_chart(fig_comp, use_container_width=True)

st.divider()

# ── Export CSV ──
st.subheader("Export des données")
import pandas as pd
df = pd.DataFrame({
    't (s)':       res['t'],
    'theta (rad)': res['theta'],
    'omega (r/s)': res['omega'],
    'g(t)':        res['g'],
    'M (N.m)':     res['M'],
    'Ec (J)':      res['Ec'],
    'Ee (J)':      res['Ee'],
})
csv = df.to_csv(index=False).encode('utf-8')
st.download_button(
    label="Télécharger les données (CSV)",
    data=csv,
    file_name="guielle_simulation.csv",
    mime="text/csv",
)

# ── Récapitulatif paramètres ──
with st.expander("Récapitulatif des paramètres utilisés"):
    pcol1, pcol2 = st.columns(2)
    with pcol1:
        st.write(f"**n** = {p_used['n']} ressorts")
        st.write(f"**r** = {p_used['r']} m")
        st.write(f"**l₀** = {p_used['l']} m")
        st.write(f"**k** = {p_used['k']} N/m")
        st.write(f"**ε** = {1 + p_used['l']/p_used['r']:.3f}")
    with pcol2:
        st.write(f"**m** = {p_used['m']} kg")
        st.write(f"**α** = {p_used['alpha']} → I = {p_used['alpha']*p_used['m']*p_used['r']**2:.4f} kg·m²")
        st.write(f"**b** = {p_used['b']}")
        st.write(f"**Cs** = {p_used['Cs']} N·m")
        st.write(f"**ω₀** = {p_used['omega0']} rad/s")

st.caption("Simulateur du Système de Guielle · Formules fidèles au document original · Intégration RK4 · H. GUIELLE 2025")
