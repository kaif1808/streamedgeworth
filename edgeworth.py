import streamlit as st
import numpy as np
import plotly.graph_objects as go
from scipy.optimize import minimize, brentq, minimize_scalar
import re
import warnings

# --- Configuration & Setup ---
st.set_page_config(
    layout="wide", 
    page_title="Edgeworth Box Simulator",
    page_icon="📊"
)
warnings.filterwarnings('ignore')

# Custom CSS for cleaner look
st.markdown("""
<style>
    .block-container {padding-top: 2rem; padding-bottom: 2rem;}
    h1, h2, h3 {font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif;}
</style>
""", unsafe_allow_html=True)

# --- Helper Functions ---
def parse_latex_to_numpy(latex_str):
    if not latex_str: return "0"
    expr = latex_str.lower().replace("^", "**").replace(r"\cdot", "*")
    replacements = {
        r"\\ln": "np.log", r"\\log": "np.log", r"\\exp": "np.exp",
        r"\\sqrt": "np.sqrt", r"\\min": "np.minimum", r"\\max": "np.maximum",
        r"min": "np.minimum", r"max": "np.maximum",
    }
    for tex, py in replacements.items(): expr = re.sub(tex, py, expr)
    expr = expr.replace("{", "(").replace("}", ")")
    expr = re.sub(r'(\d)([xy])', r'\1*\2', expr)
    return expr

def evaluate_custom_utility(x, y, formula):
    try:
        env = {'x': x, 'y': y, 'np': np, 'abs': np.abs, 'log': np.log, 
               'exp': np.exp, 'sqrt': np.sqrt, 'minimum': np.minimum, 'maximum': np.maximum}
        return eval(parse_latex_to_numpy(formula), {"__builtins__": None}, env)
    except (SyntaxError, NameError, TypeError, ZeroDivisionError):
        return np.zeros_like(x) if isinstance(x, np.ndarray) else 0

# Helper for dual inputs (Slider + Number Box + Buttons)
def dual_input(label, min_val, max_val, default_val, step=0.1, key_suffix=""):
    k_slider = f"slider_{key_suffix}"
    k_num = f"num_{key_suffix}"
    
    if k_slider not in st.session_state:
        st.session_state[k_slider] = default_val
    if k_num not in st.session_state:
        st.session_state[k_num] = default_val

    def update_from_slider(): st.session_state[k_num] = st.session_state[k_slider]
    def update_from_num(): st.session_state[k_slider] = st.session_state[k_num]
    def increment():
        new_val = min(max_val, st.session_state[k_num] + step)
        st.session_state[k_num] = new_val
        st.session_state[k_slider] = new_val
    def decrement():
        new_val = max(min_val, st.session_state[k_num] - step)
        st.session_state[k_num] = new_val
        st.session_state[k_slider] = new_val

    st.markdown(f"<small>{label}</small>", unsafe_allow_html=True)
    c1, c2, c3, c4 = st.columns([0.15, 0.45, 0.25, 0.15])
    
    with c1: st.button("➖", key=f"dec_{key_suffix}", on_click=decrement)
    with c2: st.slider("", min_value=float(min_val), max_value=float(max_val), key=k_slider, step=step, on_change=update_from_slider, label_visibility="collapsed")
    with c3: st.number_input("", min_value=float(min_val), max_value=float(max_val), key=k_num, step=step, format="%.2f", on_change=update_from_num, label_visibility="collapsed")
    with c4: st.button("➕", key=f"inc_{key_suffix}", on_click=increment)
        
    return st.session_state[k_num]

# --- Utility Logic ---
def utility_func(x, y, u_type, params):
    x = np.maximum(x, 1e-9)
    y = np.maximum(y, 1e-9)

    if u_type == "Custom (Enter Formula)":
        return evaluate_custom_utility(x, y, params.get('formula', 'x*y'))

    alpha = params.get('alpha', 0.5)
    beta = params.get('beta', 0.5)
    a = params.get('a', 0.0)
    b = params.get('b', 0.0)
    
    if u_type == "Cobb-Douglas":
        return (x ** alpha) * (y ** beta)
    elif u_type == "Perfect Substitutes":
        return alpha * x + beta * y
    elif u_type == "Perfect Complements (Min)": 
        return np.minimum(alpha * x, beta * y)
    elif u_type == "Max Preferences (Convex)": 
        return np.maximum(alpha * x, beta * y)
    elif u_type == "Quasi-Linear (Shifted Product)": 
        return (x + a) * (y + b) 
    elif u_type == "Satiation (Bliss Point)": 
        return -1 * ((x - a)**2 + (y - b)**2)
    elif u_type == "Mixed Cobb-Douglas": 
        return x * (y ** alpha)
    return 0

def check_preference_curvature(u_type, params):
    if u_type in ["Cobb-Douglas", "Perfect Complements (Min)", "Mixed Cobb-Douglas", "Quasi-Linear (Shifted Product)", "Satiation (Bliss Point)"]:
        return "PREFER_AVERAGES"
    if u_type in ["Max Preferences (Convex)"]:
        return "PREFER_EXTREMES"
    if u_type in ["Perfect Substitutes"]:
        return "LINEAR"
        
    # Custom: Check along a budget-like line (1,3) to (3,1)
    try:
        u1 = utility_func(1.0, 3.0, u_type, params)
        u2 = utility_func(3.0, 1.0, u_type, params)
        u_mid = utility_func(2.0, 2.0, u_type, params)
        avg_u = (u1 + u2) / 2.0
        
        if u_mid > avg_u + 1e-6:
             return "PREFER_AVERAGES"
        if u_mid < avg_u - 1e-6:
             return "PREFER_EXTREMES"
        return "LINEAR"
    except:
        return "PREFER_AVERAGES"

def calculate_mrs(x, y, u_type, params):
    h = 1e-5
    u0 = utility_func(x, y, u_type, params)
    ux = (utility_func(x + h, y, u_type, params) - u0) / h
    uy = (utility_func(x, y + h, u_type, params) - u0) / h
    
    if abs(uy) < 1e-9:
        if abs(ux) < 1e-9: return 0 
        return np.inf
    return ux / uy

def verify_pareto_efficiency(xA, yA, total_x, total_y, type_A, params_A, type_B, params_B):
    mrs_A = calculate_mrs(xA, yA, type_A, params_A)
    mrs_B = calculate_mrs(total_x - xA, total_y - yA, type_B, params_B)
    
    # Interior
    if xA > 1e-3 and xA < total_x - 1e-3 and yA > 1e-3 and yA < total_y - 1e-3:
        if np.isinf(mrs_A) or np.isinf(mrs_B): return True
        return abs(mrs_A - mrs_B) < 0.2 # Looser tolerance for numerical stability
        
    # Left Wall (xA ~ 0). A has 0 X. A Buy X?
    # Trade if MRS_A > MRS_B. Efficient if MRS_A <= MRS_B.
    if xA <= 1e-3:
        if np.isinf(mrs_B): return True 
        return mrs_A <= mrs_B + 1e-2

    # Right Wall (xA ~ total). A has All X. A Sell X?
    # Trade if MRS_B > MRS_A. Efficient if MRS_B <= MRS_A.
    if xA >= total_x - 1e-3:
        if np.isinf(mrs_A): return True
        return mrs_A >= mrs_B - 1e-2
        
    # Bottom Wall (yA ~ 0). A has 0 Y. A Buy Y (Sell X).
    # Trade if MRS_B > MRS_A. Efficient if MRS_A >= MRS_B.
    if yA <= 1e-3:
         return mrs_A >= mrs_B - 1e-2
         
    # Top Wall (yA ~ total). A has All Y. A Sell Y (Buy X).
    # Trade if MRS_A > MRS_B. Efficient if MRS_A <= MRS_B.
    if yA >= total_y - 1e-3:
         return mrs_A <= mrs_B + 1e-2
         
    return True

# --- Solver Logic ---
def get_demand(u_type, params, px, py, income, total_x_limit=None, total_y_limit=None):
    """Calculate optimal bundle (x, y) given prices and income."""
    # 1. Analytical Solutions for Standard Types
    alpha = params.get('alpha', 0.5)
    beta = params.get('beta', 0.5)
    
    if u_type in ["Cobb-Douglas", "Mixed Cobb-Douglas"]:
        # CD: x = (alpha/(alpha+beta)) * I / px
        # Mixed CD: U = x * y^alpha -> equivalent to alpha=1, beta=alpha
        if u_type == "Mixed Cobb-Douglas":
            eff_alpha, eff_beta = 1.0, alpha
        else:
            eff_alpha, eff_beta = alpha, beta
            
        x = (eff_alpha / (eff_alpha + eff_beta)) * income / px
        y = (eff_beta / (eff_alpha + eff_beta)) * income / py
        return x, y

    elif u_type == "Perfect Substitutes":
        # MRS = alpha/beta. If px/py < MRS, buy all X. If >, buy all Y.
        mrs = alpha / beta
        price_ratio = px / py
        
        if price_ratio < mrs - 1e-6:
            return income / px, 0.0
        elif price_ratio > mrs + 1e-6:
            return 0.0, income / py
        else:
            return income / px, 0.0 

    elif u_type == "Perfect Complements (Min)":
        x = income / (px + py * (alpha / beta))
        y = (alpha / beta) * x
        return x, y

    elif u_type == "Quasi-Linear (Shifted Product)":
        a = params.get('a', 0.0)
        b = params.get('b', 0.0)
        I_eff = income + px*a + py*b
        
        X = I_eff / (2 * px)
        Y = I_eff / (2 * py)
        
        x = max(0, X - a)
        y = (income - px*x) / py
        return x, y

    # 2. Numerical Solution for Others (Satiation, Custom, Max Prefs)
    if u_type == "Max Preferences (Convex)":
        x1, y1 = income / px, 0
        x2, y2 = 0, income / py
        u1 = utility_func(x1, y1, u_type, params)
        u2 = utility_func(x2, y2, u_type, params)
        return (x1, y1) if u1 >= u2 else (x2, y2)

    # General Numerical Solver
    def obj(v): return -utility_func(v[0], v[1], u_type, params)
    def con_budget(v): return income - (px*v[0] + py*v[1])
    
    b_x = (0, total_x_limit) if total_x_limit else (0, None)
    b_y = (0, total_y_limit) if total_y_limit else (0, None)
    
    x0 = income / (2 * px)
    y0 = income / (2 * py)
    
    res = minimize(obj, [x0, y0], bounds=[b_x, b_y], constraints={'type':'ineq', 'fun':con_budget}, tol=1e-5)
    if res.success:
        return res.x[0], res.x[1]
    
    return x0, y0

def solve_walrasian_equilibrium(total_x, total_y, type_A, params_A, type_B, params_B, endow_A, endow_B):
    # Normalize py = 1. Solve for px.
    py = 1.0
    wAx, wAy = endow_A
    wBx, wBy = endow_B
    
    def excess_demand_x(px):
        if px <= 0: return 1e9 # Penalty for negative price
        
        IA = px * wAx + py * wAy
        IB = px * wBx + py * wBy
        
        xA, yA = get_demand(type_A, params_A, px, py, IA, total_x, total_y)
        xB, yB = get_demand(type_B, params_B, px, py, IB, total_x, total_y)
        
        return (xA + xB) - total_x

    low, high = 0.01, 100.0
    try:
        px_eq = brentq(excess_demand_x, low, high, xtol=1e-4)
    except ValueError:
        res = minimize_scalar(lambda p: abs(excess_demand_x(p)), bounds=(0.01, 100.0), method='bounded')
        px_eq = res.x
    
    IA = px_eq * wAx + py * wAy
    xA, yA = get_demand(type_A, params_A, px_eq, py, IA, total_x, total_y)
    
    return px_eq, (xA, yA)

def solve_contract_curve(total_x, total_y, type_A, params_A, type_B, params_B, uA_w, uB_w, Z_B_min, Z_B_max):
    pareto_x, pareto_y, core_x, core_y = [], [], [], []
    if Z_B_max <= Z_B_min: return pareto_x, pareto_y, core_x, core_y

    curv_A = check_preference_curvature(type_A, params_A)
    curv_B = check_preference_curvature(type_B, params_B)
    
    # Strategy: If either prefers extremes, enable full candidate search
    check_boundaries = (curv_A == "PREFER_EXTREMES" or curv_B == "PREFER_EXTREMES" or 
                        type_A == "Custom (Enter Formula)" or type_B == "Custom (Enter Formula)")

    steps = 50 # Increased precision
    levels_B = np.linspace(Z_B_min, Z_B_max, steps)
    last_x = [total_x / 2, total_y / 2] 

    for ub_val in levels_B:
        candidates = []

        # 1. Interior Optimizer
        def obj(v): return -utility_func(v[0], v[1], type_A, params_A)
        def con(v): return utility_func(total_x - v[0], total_y - v[1], type_B, params_B) - ub_val
        
        bnds = ((0, total_x), (0, total_y))
        res = minimize(obj, last_x, bounds=bnds, constraints={'type':'ineq', 'fun':con}, tol=1e-5)
        
        if res.success:
            candidates.append(res.x)
            last_x = res.x
        else:
            starts = [[0, 0], [total_x, total_y], [0, total_y], [total_x, 0]]
            for s in starts:
                res_retry = minimize(obj, s, bounds=bnds, constraints={'type':'ineq', 'fun':con}, tol=1e-5)
                if res_retry.success:
                     candidates.append(res_retry.x)
                     last_x = res_retry.x
                     break

        # 2. Boundary Candidates (if needed)
        if check_boundaries:
            # Corners
            candidates.append([0, 0])
            candidates.append([total_x, 0])
            candidates.append([0, total_y])
            candidates.append([total_x, total_y])
            
            # Edges
            def f_left(y): return utility_func(total_x, total_y - y, type_B, params_B) - ub_val
            try: 
                y_sol = brentq(f_left, 0, total_y)
                candidates.append([0, y_sol])
            except: pass
            
            def f_right(y): return utility_func(0, total_y - y, type_B, params_B) - ub_val
            try:
                y_sol = brentq(f_right, 0, total_y)
                candidates.append([total_x, y_sol])
            except: pass

            def f_bottom(x): return utility_func(total_x - x, total_y, type_B, params_B) - ub_val
            try:
                x_sol = brentq(f_bottom, 0, total_x)
                candidates.append([x_sol, 0])
            except: pass
            
            def f_top(x): return utility_func(total_x - x, 0, type_B, params_B) - ub_val
            try:
                x_sol = brentq(f_top, 0, total_x)
                candidates.append([x_sol, total_y])
            except: pass

        # 3. Selection
        best_p = None
        best_u = -np.inf
        
        for cand in candidates:
            cx, cy = cand
            cx = np.clip(cx, 0, total_x)
            cy = np.clip(cy, 0, total_y)
            
            ub_real = utility_func(total_x - cx, total_y - cy, type_B, params_B)
            
            if ub_real >= ub_val - 1e-2: 
                ua_val = utility_func(cx, cy, type_A, params_A)
                if ua_val > best_u:
                    best_u = ua_val
                    best_p = [cx, cy]

        if best_p is not None:
            if verify_pareto_efficiency(best_p[0], best_p[1], total_x, total_y, type_A, params_A, type_B, params_B):
                pareto_x.append(best_p[0])
                pareto_y.append(best_p[1])
                if best_u >= uA_w - 1e-3 and ub_real >= uB_w - 1e-3:
                    core_x.append(best_p[0])
                    core_y.append(best_p[1])

    if pareto_x:
        p_points = sorted(zip(pareto_x, pareto_y), key=lambda k: k[0])
        pareto_x, pareto_y = zip(*p_points)
        pareto_x, pareto_y = list(pareto_x), list(pareto_y)

    if core_x:
        c_points = sorted(zip(core_x, core_y), key=lambda k: k[0])
        core_x, core_y = zip(*c_points)
        core_x, core_y = list(core_x), list(core_y)

    return pareto_x, pareto_y, core_x, core_y

# --- Plotting Logic ---
def get_theme_config(theme_name, dark_mode):
    if theme_name == "Modern Professional":
        return {
            "font": "Arial, sans-serif",
            "bg": "#1a1a1a" if dark_mode else "#ffffff",
            "grid": "#333" if dark_mode else "#e5e7eb",
            "text": "#e0e0e0" if dark_mode else "#374151",
            "A": "#d32f2f", # Strong Red
            "B": "#1976d2", # Strong Blue
            "Pareto": "#2e7d32", # Green
            "Core": "#fbc02d", # Amber/Gold
            "Lens": "rgba(46, 125, 50, 0.08)",
            "EndowLineA": "rgba(211, 47, 47, 0.8)",
            "EndowLineB": "rgba(25, 118, 210, 0.8)",
            "origin_A": "#d32f2f",
            "origin_B": "#1976d2"
        }
    else: # Classic Textbook
        return {
            "font": "Times New Roman, serif",
            "bg": "#121212" if dark_mode else "#fcfcfc",
            "grid": "#444" if dark_mode else "#e0e0e0",
            "text": "white" if dark_mode else "black",
            "A": "rgba(180, 0, 0, 0.9)",
            "B": "rgba(0, 0, 180, 0.9)",
            "Pareto": "#388e3c",
            "Core": "#ffa000",
            "Lens": "rgba(100, 100, 100, 0.1)",
            "EndowLineA": "rgba(180, 0, 0, 1)",
            "EndowLineB": "rgba(0, 0, 180, 1)",
            "origin_A": "#800000",
            "origin_B": "#000080"
        }

def plot_edgeworth_box(Z_A, Z_B, x_vec, y_vec, total_x, total_y, 
                       pareto_x, pareto_y, core_x, core_y, 
                       uA_w, uB_w, endow_x, endow_y, 
                       settings, theme_config, we_data=None):
    
    colors = theme_config
    ft_font = theme_config["font"]
    
    fig = go.Figure()

    # 1. Exchange Lens
    if settings.get("show_lens", True):
        lens_mask = np.logical_and(Z_A >= uA_w - 1e-4, Z_B >= uB_w - 1e-4).astype(int)
        fig.add_trace(go.Contour(
            z=lens_mask, x=x_vec, y=y_vec,
            showscale=False,
            contours=dict(start=0.5, end=0.5, coloring='fill'),
            colorscale=[[0, 'rgba(0,0,0,0)'], [1, colors["Lens"]]],
            line=dict(width=0),
            name="Exchange Lens",
            hoverinfo='skip'
        ))

    # 2. Endowment ICs
    if settings.get("show_endow", True):
        fig.add_trace(go.Contour(
            z=Z_A, x=x_vec, y=y_vec, showscale=False,
            contours=dict(type='constraint', value=uA_w, coloring='lines'),
            line=dict(width=3, color=colors["EndowLineA"], dash=settings.get("style_A", "solid")),
            name="UA(ω)", hovertemplate="UA = %{z:.2f}<extra></extra>"
        ))
        fig.add_trace(go.Contour(
            z=Z_B, x=x_vec, y=y_vec, showscale=False,
            contours=dict(type='constraint', value=uB_w, coloring='lines'),
            line=dict(width=3, color=colors["EndowLineB"], dash=settings.get("style_B", "dot")),
            name="UB(ω)", hovertemplate="UB = %{z:.2f}<extra></extra>"
        ))

    # 3. General ICs
    ic_mode = settings.get("ic_mode", "Auto (Density)")
    style_A = settings.get("style_A", "solid")
    style_B = settings.get("style_B", "dot")
    
    if ic_mode == "Auto (Density)":
        n_curves = settings.get("n_curves", 30)
        if settings.get("show_curves_A", True):
            fig.add_trace(go.Contour(
                z=Z_A, x=x_vec, y=y_vec, colorscale='Reds', showscale=False, ncontours=n_curves, 
                contours=dict(coloring='lines', showlabels=False),
                line=dict(width=1, color=colors["A"], dash=style_A), name="UA Map",
                hovertemplate="UA = %{z:.2f}<extra></extra>"
            ))
        if settings.get("show_curves_B", True):
            fig.add_trace(go.Contour(
                z=Z_B, x=x_vec, y=y_vec, colorscale='Blues', showscale=False, ncontours=n_curves,
                contours=dict(coloring='lines', showlabels=False),
                line=dict(width=1, color=colors["B"], dash=style_B), name="UB Map",
                hovertemplate="UB = %{z:.2f}<extra></extra>"
            ))
    else:
        if settings.get("show_curves_A", True):
            count_A = settings.get("n_curves_A", 10)
            levels_A = np.linspace(np.min(Z_A), np.max(Z_A), count_A + 2)[1:-1]
            for val in levels_A:
                fig.add_trace(go.Contour(
                    z=Z_A, x=x_vec, y=y_vec, showscale=False,
                    contours=dict(type='constraint', value=val, coloring='lines'),
                    line=dict(width=1, color=colors["A"], dash=style_A), showlegend=False, hoverinfo='skip'
                ))
        if settings.get("show_curves_B", True):
            count_B = settings.get("n_curves_B", 10)
            levels_B = np.linspace(np.min(Z_B), np.max(Z_B), count_B + 2)[1:-1]
            for val in levels_B:
                fig.add_trace(go.Contour(
                    z=Z_B, x=x_vec, y=y_vec, showscale=False,
                    contours=dict(type='constraint', value=val, coloring='lines'),
                    line=dict(width=1, color=colors["B"], dash=style_B), showlegend=False, hoverinfo='skip'
                ))

    # 4. Pareto Set & Core
    line_mode_enabled = settings.get("line_mode", False)
    def curve_mode(points):
        return 'lines' if line_mode_enabled and len(points) >= 2 else 'markers'

    if settings.get("show_pareto", True) and pareto_x:
        fig.add_trace(go.Scatter(
            x=pareto_x, y=pareto_y, mode=curve_mode(pareto_x), 
            marker=dict(size=6, color=colors["Pareto"], line=dict(width=1, color="white")), 
            line=dict(width=4, color=colors["Pareto"]), name="Pareto Set",
            hovertemplate="Pareto<br>x: %{x:.2f}<br>y: %{y:.2f}<extra></extra>"
        ))

    if settings.get("show_core", True) and core_x:
        fig.add_trace(go.Scatter(
            x=core_x, y=core_y, mode=curve_mode(core_x), 
            marker=dict(size=9, color=colors["Core"], line=dict(width=1, color="white")), 
            line=dict(width=8, color=colors["Core"]), name="The Core",
            hovertemplate="Core<br>x: %{x:.2f}<br>y: %{y:.2f}<extra></extra>"
        ))

    # 5. Endowment Point
    if settings.get("show_endow", True):
        fig.add_trace(go.Scatter(
            x=[endow_x], y=[endow_y], mode='markers', 
            marker=dict(color=colors["text"], size=14, line=dict(width=2, color=colors["bg"])), 
            name="Endowment", hovertemplate="Endowment (ω)<br>x: %{x}<br>y: %{y}<extra></extra>"
        ))
        fig.add_annotation(x=endow_x, y=endow_y, text="ω", font=dict(size=18, color=colors["text"], weight="bold", family=ft_font), ax=15, ay=-15)

    # 6. Walrasian Equilibrium
    if settings.get("show_we", False) and we_data:
        px_eq, (xA_eq, yA_eq) = we_data
        
        x_range = np.linspace(0, total_x, 100)
        y_line = endow_y - px_eq * (x_range - endow_x)
        
        mask = (y_line >= 0) & (y_line <= total_y)
        
        if np.any(mask):
            fig.add_trace(go.Scatter(
                x=x_range[mask], y=y_line[mask], mode='lines',
                line=dict(color=colors["text"], width=2, dash='dashdot'),
                name=f"Budget Line (p={px_eq:.2f})",
                hoverinfo='skip'
            ))
            
        fig.add_trace(go.Scatter(
            x=[xA_eq], y=[yA_eq], mode='markers',
            marker=dict(color='#9c27b0', size=12, symbol='diamond', line=dict(width=1, color='white')), 
            name="Walrasian Eq.",
            hovertemplate="WE Allocation<br>xA: %{x:.2f}<br>yA: %{y:.2f}<extra></extra>"
        ))

    aspect_ratio = total_y / total_x if total_x > 0 else 1
    max_width = 900
    max_height = 700
    
    if aspect_ratio > (max_height / max_width):
        calc_height = max_height
        calc_width = max_height / aspect_ratio
    else:
        calc_width = max_width
        calc_height = max_width * aspect_ratio

    fig.update_layout(
        template="plotly_white" if not dark_mode else "plotly_dark",
        width=calc_width, height=calc_height,
        xaxis=dict(
            title=dict(text="Good X (Agent A)", font=dict(size=16)), 
            range=[0, total_x], constrain='domain', 
            gridcolor=colors["grid"], gridwidth=1, showline=True, linewidth=2, linecolor=colors["text"], mirror=True,
            tickfont=dict(color=colors["text"], family=ft_font), title_font=dict(color=colors["text"], family=ft_font)
        ),
        yaxis=dict(
            title=dict(text="Good Y (Agent A)", font=dict(size=16)), 
            range=[0, total_y], scaleanchor="x", scaleratio=1, 
            gridcolor=colors["grid"], gridwidth=1, showline=True, linewidth=2, linecolor=colors["text"], mirror=True,
            tickfont=dict(color=colors["text"], family=ft_font), title_font=dict(color=colors["text"], family=ft_font)
        ),
        plot_bgcolor=colors["bg"], paper_bgcolor=colors["bg"],
        legend=dict(orientation="h", y=1.02, x=0.5, xanchor="center", bgcolor='rgba(0,0,0,0)', font=dict(color=colors["text"], family=ft_font)),
        margin=dict(l=50, r=50, t=80, b=50),
        hovermode="closest"
    )

    fig.add_annotation(x=0, y=0, text="Origin A", font=dict(size=16, color=colors["origin_A"], weight="bold", family=ft_font), showarrow=False, xanchor="right", yanchor="top", xshift=-5, yshift=-5)
    fig.add_annotation(x=total_x, y=total_y, text="Origin B", font=dict(size=16, color=colors["origin_B"], weight="bold", family=ft_font), showarrow=False, xanchor="left", yanchor="bottom", xshift=5, yshift=5)

    return fig

# --- Main App Logic ---
st.title("Edgeworth Box Simulator")

# Mode Selection
mode = st.sidebar.radio("Mode", ["Continuous Exchange (Edgeworth Box)", "Indivisible Object (Discrete)"], horizontal=True)

if mode == "Indivisible Object (Discrete)":
    st.markdown("### Discrete Allocation: Indivisible Object + Money")
    st.sidebar.header("⚙️ Discrete Setup")
    
    n_agents = st.sidebar.number_input("Number of Agents", 2, 3, 3)
    total_money = st.sidebar.number_input("Total Money", 10, 1000, 100)
    
    agents = []
    cols = st.columns(n_agents)
    for i in range(n_agents):
        name = chr(65+i) # A, B, C
        with cols[i]:
            st.subheader(f"Agent {name}")
            val = st.number_input(f"Valuation {name} (€)", 0, 1000, (i+1)*10, key=f"v_{i}")
            # We don't strictly need initial money endowment for efficiency, but we do for core/fairness.
            # For simple efficiency, just valuation matters.
            # But let's add endowment to check IR constraints if we want (future).
            # For now, user inputs total money, we can distribute it arbitrarily or evenly for display?
            # The problem statement usually gives initial endowments.
            money = st.number_input(f"Initial Money {name}", 0, int(total_money), int(total_money/n_agents), key=f"m_{i}")
            agents.append({"name": name, "val": val, "money": money})
            
    # Calculate Efficient Allocation
    # Efficiency: Agent with highest valuation holds object.
    sorted_agents = sorted(agents, key=lambda x: x['val'], reverse=True)
    highest_val_agent = sorted_agents[0]
    
    # Check for ties in highest valuation
    max_val = highest_val_agent['val']
    efficient_holders = [a for a in agents if a['val'] == max_val]
    
    st.markdown("---")
    st.subheader("🏆 Efficiency Analysis")
    
    if len(efficient_holders) == 1:
        st.success(f"**Pareto Efficient Outcome:** Agent **{efficient_holders[0]['name']}** holds the object.")
        st.write(f"Reason: Agent {efficient_holders[0]['name']} has the strictly highest valuation ({efficient_holders[0]['val']}).")
    else:
        names = ", ".join([a['name'] for a in efficient_holders])
        st.success(f"**Pareto Efficient Outcome:** Any of **{names}** can hold the object.")
        st.write(f"Reason: Agents {names} have the tied highest valuation ({max_val}).")
        
    st.info("In a quasilinear setting with transferrable utility (money), efficiency requires the object goes to the person who values it most, regardless of initial money distribution (provided transfers are possible).")

    # Core Analysis (Simple)
    st.markdown("### 🧩 Core Analysis (Payoff Space)")
    st.write("The core consists of efficient allocations where no coalition can improve upon their initial endowment.")
    
    # This is complex to visualize generally, but we can show the condition.
    # Condition: Agent i gets utility u_i >= u_i(endowment).
    # u_i(endow) = money_i + (val_i if holds else 0)
    
    st.write("#### Individual Rationality Constraints:")
    for a in agents:
        # Does agent hold object initially? 
        # We didn't ask who holds it initially. Let's assume nobody or ask?
        # Problem 4 says "There are 100 euros and one object". 
        # Usually ownership is defined. Let's add a selector for initial owner.
        pass
        
    initial_owner_idx = st.sidebar.selectbox("Initial Object Owner", range(n_agents), format_func=lambda x: agents[x]['name'])
    initial_owner = agents[initial_owner_idx]
    
    st.write(f"**Initial State:** {initial_owner['name']} has the object.")
    
    # Calculate Reservation Utilities
    res_utils = {}
    for i, a in enumerate(agents):
        is_owner = (i == initial_owner_idx)
        u_res = a['money'] + (a['val'] if is_owner else 0)
        res_utils[a['name']] = u_res
        st.write(f"- **Agent {a['name']}**: Min Utility = {u_res} (Money: {a['money']} + Item: {a['val'] if is_owner else 0})")
        
    st.markdown("#### Efficient Allocations in Core:")
    st.write(f"Agent **{efficient_holders[0]['name']}** gets object. Transfers must satisfy IR.")
    
    # Show simple table of condition
    # If A is efficient holder:
    # u_A = m_A' + v_A >= u_res_A
    # u_B = m_B' >= u_res_B
    # ...
    # Sum m' = Total Money
    
    holder = efficient_holders[0]
    st.write(f"If **{holder['name']}** holds object:")
    
    conditions = []
    min_transfer_needed = 0
    
    for a in agents:
        if a['name'] == holder['name']:
            # Holder condition: m' + v >= u_res
            # m' >= u_res - v
            min_m = res_utils[a['name']] - a['val']
            conditions.append(f"Agent {a['name']} Final Money ≥ {min_m}")
        else:
            # Non-holder: m' >= u_res
            min_m = res_utils[a['name']]
            conditions.append(f"Agent {a['name']} Final Money ≥ {min_m}")
            
    for c in conditions:
        st.write(f"- {c}")
        
else:
    # CONTINUOUS MODE (Existing Logic)
    st.sidebar.header("⚙️ Configuration")

    # Presets
    presets = {
        "Custom": {},
        "PS8 Q1 (2025): Shifted CD": {
            "dim": (5, 10), "endow": (3, 3),
            "A": {"type": "Quasi-Linear (Shifted Product)", "b": 3.0},
            "B": {"type": "Quasi-Linear (Shifted Product)", "b": 2.0}
        },
        "PS8 Q1 (Old): Convex vs Linear": {
            "dim": (10, 10), "endow": (7, 6),
            "A": {"type": "Custom (Enter Formula)", "formula": "x**2 + y**2"},
            "B": {"type": "Perfect Substitutes", "alpha": 1.0, "beta": 1.0}
        },
        "PS8 Q2 (2025): Mixed CD vs Perf Subs": {
            "dim": (12, 12), "endow": (3, 3),
            "A": {"type": "Mixed Cobb-Douglas", "alpha": 3.0},
            "B": {"type": "Perfect Substitutes", "alpha": 1.0, "beta": 1.0}
        },
        "PS8 Q3 (2025): Min vs Max": {
            "dim": (6, 6), "endow": (4, 1),
            "A": {"type": "Perfect Complements (Min)", "alpha": 1.0, "beta": 1.0},
            "B": {"type": "Max Preferences (Convex)", "alpha": 1.0, "beta": 1.0}
        },
        "PS8 Q5 (2025): Leontief vs Perf Subs": {
            "dim": (10, 10), "endow": (4, 4),
            "A": {"type": "Perfect Complements (Min)", "alpha": 1.0, "beta": 1.0},
            "B": {"type": "Perfect Substitutes", "alpha": 1.0, "beta": 1.0}
        },
        "PS8 Q6 (2025): Satiation": {
            "dim": (10, 10), "endow": (4, 8),
            "A": {"type": "Satiation (Bliss Point)", "a": 3.0, "b": 3.0},
            "B": {"type": "Cobb-Douglas", "alpha": 1.0, "beta": 1.0}
        },
        "Gen Eq Notes: Core Example": {
            "dim": (10, 10), "endow": (6, 3),
            "A": {"type": "Quasi-Linear (Shifted Product)", "a": 0.0, "b": 2.0}, # x(y+2)
            "B": {"type": "Cobb-Douglas", "alpha": 1.0, "beta": 1.0} # xy
        }
    }

    if st.sidebar.button("Reset to Baseline"):
        # Reset Dimensions
        st.session_state["slider_dim_x"] = 10.0
        st.session_state["num_dim_x"] = 10.0
        st.session_state["slider_dim_y"] = 10.0
        st.session_state["num_dim_y"] = 10.0
        
        # Reset Endowment
        st.session_state["slider_endow_x"] = 5.0
        st.session_state["num_endow_x"] = 5.0
        st.session_state["slider_endow_y"] = 5.0
        st.session_state["num_endow_y"] = 5.0
        
        # Reset Utility Types
        st.session_state["type_A"] = "Cobb-Douglas"
        st.session_state["type_B"] = "Cobb-Douglas"
        
        # Reset Parameters
        st.session_state["slider_al_A"] = 1.0
        st.session_state["num_al_A"] = 1.0
        st.session_state["slider_be_A"] = 1.0
        st.session_state["num_be_A"] = 1.0
        
        st.session_state["slider_al_B"] = 1.0
        st.session_state["num_al_B"] = 1.0
        st.session_state["slider_be_B"] = 1.0
        st.session_state["num_be_B"] = 1.0
        
        st.rerun()

    selected_preset = st.sidebar.selectbox("Load Scenario", list(presets.keys()))
    p_data = presets[selected_preset]

    # Configuration Sections
    with st.sidebar.expander("📦 Dimensions & Endowment", expanded=True):
        tx_def = p_data.get("dim", (10, 10))[0]
        ty_def = p_data.get("dim", (10, 10))[1]
        total_x = dual_input("Total Good X", 1.0, 100.0, float(tx_def), 1.0, "dim_x")
        total_y = dual_input("Total Good Y", 1.0, 100.0, float(ty_def), 1.0, "dim_y")
        
        # Clamping logic for endowments
        if "slider_endow_x" in st.session_state and st.session_state["slider_endow_x"] > total_x:
            st.session_state["slider_endow_x"] = total_x
            st.session_state["num_endow_x"] = total_x
            
        if "slider_endow_y" in st.session_state and st.session_state["slider_endow_y"] > total_y:
            st.session_state["slider_endow_y"] = total_y
            st.session_state["num_endow_y"] = total_y

        st.markdown("---")
        ex_def = p_data.get("endow", (total_x*0.7, total_y*0.6))[0]
        ey_def = p_data.get("endow", (total_x*0.7, total_y*0.6))[1]
        with st.container():
            st.markdown("**Agent A Endowment**")
            endow_x = dual_input("ω_x", 0.0, total_x, float(ex_def), 0.1, "endow_x")
            endow_y = dual_input("ω_y", 0.0, total_y, float(ey_def), 0.1, "endow_y")
        endow_B_x, endow_B_y = total_x - endow_x, total_y - endow_y

    def config_agent(name, prefix, data, expanded=False):
        with st.sidebar.expander(f"👤 Agent {name} Preferences", expanded=expanded):
            default_type = data.get(name, {}).get("type", "Cobb-Douglas")
            opts = ["Cobb-Douglas", "Perfect Substitutes", "Perfect Complements (Min)", 
                    "Max Preferences (Convex)", "Quasi-Linear (Shifted Product)", 
                    "Mixed Cobb-Douglas", "Satiation (Bliss Point)", "Custom (Enter Formula)"]
            
            idx = opts.index(default_type) if default_type in opts else 0
            u_type = st.selectbox(f"Utility Type", opts, index=idx, key=f"type_{prefix}")
            
            params = {}
            defaults = data.get(name, {})
            if u_type == "Custom (Enter Formula)":
                params['formula'] = st.text_input("Formula (e.g. x^0.5 * y^0.5)", value=defaults.get('formula', "x*y"), key=f"f_{prefix}")
            elif u_type == "Satiation (Bliss Point)":
                params['a'] = dual_input("Bliss Point X", -50, 50, defaults.get('a', 3.0), 0.5, f"a_{prefix}")
                params['b'] = dual_input("Bliss Point Y", -50, 50, defaults.get('b', 3.0), 0.5, f"b_{prefix}")
            elif u_type == "Quasi-Linear (Shifted Product)":
                params['a'] = dual_input("Shift Parameter X", -50, 50, defaults.get('a', 0.0), 0.5, f"qa_{prefix}")
                params['b'] = dual_input("Shift Parameter Y", -50, 50, defaults.get('b', 0.0), 0.5, f"qb_{prefix}")
            else:
                def_a = defaults.get('alpha', 1.0)
                def_b = defaults.get('beta', 1.0)
                params['alpha'] = dual_input("Alpha (α)", 0.1, 10.0, float(def_a), 0.1, f"al_{prefix}")
                if u_type != "Mixed Cobb-Douglas":
                    params['beta'] = dual_input("Beta (β)", 0.1, 10.0, float(def_b), 0.1, f"be_{prefix}")
            return u_type, params

    type_A, params_A = config_agent("A", "A", p_data, expanded=True)
    type_B, params_B = config_agent("B", "B", p_data, expanded=False)

    with st.sidebar.expander("🎨 Visual Settings", expanded=False):
        theme_name = st.radio("Theme", ["Modern Professional", "Classic Textbook"], horizontal=True)
        dark_mode = st.checkbox("Dark Mode", value=False)
        st.markdown("---")
        vis_settings = {}
        vis_settings["show_endow"] = st.checkbox("Show Endowment", value=True)
        vis_settings["show_core"] = st.checkbox("Show Core", value=True)
        vis_settings["show_pareto"] = st.checkbox("Show Pareto Set", value=True)
        vis_settings["show_lens"] = st.checkbox("Shade Exchange Lens", value=True)
        vis_settings["show_curves_A"] = st.checkbox("Show Curves (Agent A)", value=True)
        vis_settings["show_curves_B"] = st.checkbox("Show Curves (Agent B)", value=True)
        vis_settings["line_mode"] = st.checkbox("Connect Pareto Points", value=False)
        vis_settings["show_we"] = st.checkbox("Show Walrasian Equilibrium", value=False)
        
        st.markdown("**Line Styles**")
        c1, c2 = st.columns(2)
        style_map = {"Solid": "solid", "Dotted": "dot", "Dashed": "dash"}
        vis_settings["style_A"] = style_map[c1.selectbox("Agent A", ["Solid", "Dotted", "Dashed"], index=0, key="sA")]
        vis_settings["style_B"] = style_map[c2.selectbox("Agent B", ["Solid", "Dotted", "Dashed"], index=1, key="sB")]
        
        st.markdown("**Curve Density**")
        ic_mode = st.radio("Mode", ["Auto", "Manual"], horizontal=True, label_visibility="collapsed")
        vis_settings["ic_mode"] = "Auto (Density)" if ic_mode == "Auto" else "Manual"
        if ic_mode == "Auto":
            vis_settings["n_curves"] = st.slider("", 10, 100, 30)
        else:
            c1, c2 = st.columns(2)
            vis_settings["n_curves_A"] = c1.number_input("N (A)", 1, 50, 10)
            vis_settings["n_curves_B"] = c2.number_input("N (B)", 1, 50, 10)

    # --- Calculation ---
    N = 200 # High res
    x_vec = np.linspace(0, total_x, N)
    y_vec = np.linspace(0, total_y, N)
    X, Y = np.meshgrid(x_vec, y_vec)

    try:
        Z_A = utility_func(X, Y, type_A, params_A)
        if isinstance(Z_A, (float, int)): Z_A = np.full_like(X, Z_A)
    except: Z_A = np.zeros_like(X)

    try:
        Z_B = utility_func(total_x - X, total_y - Y, type_B, params_B)
        if isinstance(Z_B, (float, int)): Z_B = np.full_like(X, Z_B)
    except: Z_B = np.zeros_like(X)

    uA_w = utility_func(endow_x, endow_y, type_A, params_A)
    uB_w = utility_func(endow_B_x, endow_B_y, type_B, params_B)
    mrs_A = calculate_mrs(endow_x, endow_y, type_A, params_A)
    mrs_B = calculate_mrs(endow_B_x, endow_B_y, type_B, params_B)

    pareto_x, pareto_y, core_x, core_y = solve_contract_curve(
        total_x, total_y, type_A, params_A, type_B, params_B, uA_w, uB_w, np.min(Z_B), np.max(Z_B)
    )

    # Calculate Walrasian Equilibrium
    we_data = solve_walrasian_equilibrium(
        total_x, total_y, type_A, params_A, type_B, params_B, (endow_x, endow_y), (endow_B_x, endow_B_y)
    )

    # --- Output ---
    theme_config = get_theme_config(theme_name, dark_mode)
    fig = plot_edgeworth_box(
        Z_A, Z_B, x_vec, y_vec, total_x, total_y,
        pareto_x, pareto_y, core_x, core_y,
        uA_w, uB_w, endow_x, endow_y,
        vis_settings, theme_config, we_data
    )

    c_main = st.container()
    with c_main:
        st.plotly_chart(fig, use_container_width=True)

    # Analytics
    st.markdown("---")
    st.subheader("📊 Economic Analysis")

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Utility A (at ω)", f"{uA_w:.2f}")
    c2.metric("Utility B (at ω)", f"{uB_w:.2f}")
    c3.metric("MRS A", f"{mrs_A:.2f}" if not np.isinf(mrs_A) else "∞")
    c4.metric("MRS B", f"{mrs_B:.2f}" if not np.isinf(mrs_B) else "∞")

    st.markdown("")
    col1, col2 = st.columns([0.6, 0.4])

    with col1:
        st.info(f"""
        **Pareto Efficiency Status**
        
        Difference in MRS: **{abs(mrs_A - mrs_B):.4f}**
        
        { "✅ **Efficient Allocation**" if abs(mrs_A - mrs_B) < 0.05 or np.isinf(mrs_A) or np.isinf(mrs_B) else "⚠️ **Inefficient Allocation** - Trade Opportunities Exist" }
        """)

    with col2:
        if abs(mrs_A - mrs_B) >= 0.05 and not (np.isinf(mrs_A) or np.isinf(mrs_B)):
            if mrs_A > mrs_B:
                st.success("💡 **Trade Idea:** Agent A should **buy X** and **sell Y**.")
            else:
                st.success("💡 **Trade Idea:** Agent B should **buy X** and **sell Y**.")
        else:
            st.write("No mutually beneficial trade possible from this endowment.")

    if we_data:
        st.markdown("---")
        st.subheader("⚖️ Walrasian Equilibrium Analysis")
        px_eq, (xA_eq, yA_eq) = we_data
        
        cw1, cw2, cw3 = st.columns(3)
        cw1.metric("Eq. Price Ratio (Px/Py)", f"{px_eq:.2f}")
        cw2.metric("Agent A Allocation", f"({xA_eq:.2f}, {yA_eq:.2f})")
        
        net_x = xA_eq - endow_x
        net_y = yA_eq - endow_y
        trade_dir = "Buys" if net_x > 0 else "Sells"
        cw3.metric(f"Agent A Trade", f"{trade_dir} {abs(net_x):.2f} X")
        
        st.caption(f"Agent A also {'Sells' if net_y < 0 else 'Buys'} {abs(net_y):.2f} Y. Net Excess Demand for X ≈ 0.00")
