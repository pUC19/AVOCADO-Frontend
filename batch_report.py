# batch_report.py — AVOCADO Batch Report (finale, robuste Version 5)
from __future__ import annotations

import re
import uuid
import threading
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, Optional, Tuple, List
import random

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from flask import Blueprint, request, jsonify, abort, send_file

# ---- optionale SciPy-Abhängigkeiten ----
try:
    from scipy.signal import savgol_filter, find_peaks
except Exception:
    savgol_filter = None
    find_peaks = None
try:
    from scipy.stats import gaussian_kde, mannwhitneyu
except Exception:
    gaussian_kde = None
    mannwhitneyu = None # Platzhalter, falls SciPy nicht installiert ist

report_bp = Blueprint("batch_report", __name__, url_prefix="/report")
JOBS: Dict[str, Dict[str, Any]] = {}

# Farben/Style
COLOR = {"unfolding": "#f39c12", "refolding": "#1f77b4"}
GRID_ALPHA = 0.12
plt.rcParams.update({
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.titleweight": "bold",
    "axes.labelsize": 10,
    "axes.titlesize": 11,
    "legend.frameon": False,
    "figure.dpi": 300
})

# --- KDE-Glättung ---
SMOOTH_BW_FRAC = 0.08
SMOOTH_2ND_PASS = 1.5

# =========================
#         IO / READ
# =========================
def _read_table(path: Path) -> pd.DataFrame:
    """CSV/XLSX robust lesen."""
    if path.suffix.lower() == ".xlsx":
        for head in (0, None):
            try:
                df = pd.read_excel(path, header=head)
                if isinstance(df, pd.DataFrame) and not df.empty:
                    if head is None:
                        df.columns = [f"C{i}" for i in range(df.shape[1])]
                    return df
            except Exception:
                pass
        return pd.DataFrame()
    skip = 0
    try:
        with open(path, "r", encoding="utf-8") as fh:
            if fh.readline().startswith(">"):
                skip = 1
    except Exception:
        pass
    for kw in ({"sep": ",", "skiprows": skip}, {"sep": ";", "decimal": ",", "skiprows": skip}):
        try:
            df = pd.read_csv(path, **kw)
            if isinstance(df, pd.DataFrame) and not df.empty:
                return df
        except Exception:
            pass
    for kw in ({"sep": ",", "header": None, "skiprows": skip}, {"sep": ";", "decimal": ",", "header": None, "skiprows": skip}):
        try:
            df = pd.read_csv(path, **kw)
            if isinstance(df, pd.DataFrame) and not df.empty:
                df.columns = [f"C{i}" for i in range(df.shape[1])]
                return df
        except Exception:
            pass
    return pd.DataFrame()

def _is_smooth_filename(p: Path) -> bool:
    n = p.name.lower()
    return ("smooth" in n) or n.endswith(("_smooth.csv", "-smooth.csv"))

def _detect_direction_from_name(name: str) -> str:
    s = (name or "").lower()
    if re.search(r"(?:^|[_\-\s])(fw|unfold|forward|[0-9]*\s*-?f)(?:[^a-z]|$)", s):
        return "unfolding"
    if re.search(r"(?:^|[_\-\s])(rv|refold|reverse|[0-9]*\s*-?r)(?:[^a-z]|$)", s):
        return "refolding"
    return "unknown"

def _curve_direction_from_df_or_name(df: pd.DataFrame, p: Path) -> str:
    for c in df.columns:
        cl = str(c).lower()
        if "direction" in cl or "orientation" in cl:
            d = str(df[c].iloc[0]).lower()
            if any(k in d for k in ("unfold", "forward", "fw")): return "unfolding"
            if any(k in d for k in ("refold", "reverse", "rv")): return "refolding"
    return _detect_direction_from_name(p.name)

# =========================
#     FD CURVE PARSING
# =========================
def _is_force_like(arr: np.ndarray) -> bool:
    arr = arr[np.isfinite(arr)]
    if arr.size < 10: return False
    p1, p99 = np.percentile(arr, [1, 99])
    rng = p99 - p1
    return (-20 <= p1 <= 50) and (5 <= rng <= 200) and (p99 <= 200)

def _is_dist_nm_like(arr: np.ndarray) -> bool:
    arr = arr[np.isfinite(arr)]
    if arr.size < 10: return False
    p1, p99 = np.percentile(arr, [1, 99])
    return (p99 - p1) >= 150 or p99 >= 300

def _is_dist_um_like(arr: np.ndarray) -> bool:
    arr = arr[np.isfinite(arr)]
    if arr.size < 10: return False
    p1, p99 = np.percentile(arr, [1, 99])
    return (0.1 <= p1 <= 5.0) and (0.2 <= (p99 - p1) <= 20.0) and (p99 <= 50)

def _read_fd_xy(df: pd.DataFrame) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """X=Distanz (nm), Y=Force (pN) erkennen."""
    exact_x = [c for c in df.columns if str(c).strip().lower() == "relative distance, nm"]
    exact_y = [c for c in df.columns if str(c).strip().lower() == "force, pn"]
    if exact_x and exact_y:
        x = pd.to_numeric(df[exact_x[0]], errors="coerce").values
        y = pd.to_numeric(df[exact_y[0]], errors="coerce").values
    elif df.shape[1] == 2:
        a = pd.to_numeric(df.iloc[:, 0], errors="coerce").values
        b = pd.to_numeric(df.iloc[:, 1], errors="coerce").values
        if _is_force_like(a) and (_is_dist_nm_like(b) or _is_dist_um_like(b)):
            x, y = (b * 1000.0 if _is_dist_um_like(b) else b, a)
        elif _is_force_like(b) and (_is_dist_nm_like(a) or _is_dist_um_like(a)):
            x, y = (a * 1000.0 if _is_dist_um_like(a) else a, b)
        else: # Fallback
            rng = lambda u: np.ptp(u[np.isfinite(u)]) if np.any(np.isfinite(u)) else 0
            x, y = (a, b) if rng(a) > rng(b) else (b, a)
            if _is_dist_um_like(x) and not _is_dist_nm_like(x): x *= 1000.0
    else:
        num_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c].dtype)]
        if len(num_cols) < 2: return None
        x, y = df[num_cols[0]].values, df[num_cols[1]].values
    
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 10: return None
    x, y = x[m], y[m]
    if np.nanmedian(np.diff(x)) < 0: x, y = x[::-1], y[::-1]
    return x, y

# =========================
#   CURVE CLEAN & SMOOTH
# =========================
def _prepare_curve(x: np.ndarray, y: np.ndarray, target_n: int = 1400) -> Tuple[np.ndarray, np.ndarray]:
    """Sortieren, Duplikate, clippen, glätten, resamplen."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 10: return np.array([]), np.array([])
    x, y = x[m], y[m]
    
    order = np.argsort(x); x, y = x[order], y[order]
    ux, idx = np.unique(x, return_index=True); x, y = ux, y[idx]
    if x.size < 10: return x, y

    y = np.clip(y, *np.nanquantile(y, (0.005, 0.995)))
    if savgol_filter is not None and x.size >= 41:
        try:
            win = max(41, int(x.size * 0.04) // 2 * 2 + 1)
            y = savgol_filter(y, window_length=win, polyorder=3, mode="interp")
        except Exception: pass
    
    n = min(target_n, x.size)
    xg = np.linspace(x[0], x[-1], n)
    yg = np.interp(xg, x, y)
    return xg, yg

def _collect_fd_curves(analysis_folder: str | Path, max_curves: int = 2000) -> List[Tuple[np.ndarray, np.ndarray, str]]:
    folder = Path(analysis_folder)
    curves: List[Tuple[np.ndarray, np.ndarray, str]] = []
    files = [p for p in folder.iterdir() if p.is_file() and p.suffix.lower() in (".csv", ".xlsx")]
    use = [p for p in files if _is_smooth_filename(p)] or files
    use.sort(key=lambda q: (("fd curve" not in q.name.lower()), q.name.lower()))
    
    for p in use:
        if len(curves) >= max_curves: break
        df = _read_table(p)
        if df.empty: continue
        xy = _read_fd_xy(df)
        if xy is None: continue
        x, y = _prepare_curve(*xy)
        if x.size < 10: continue
        d = _curve_direction_from_df_or_name(df, p)
        curves.append((x, y, d))
    return curves

# =========================
#      STABLE DENSITIES
# =========================
def _smooth_density_1d(values: np.ndarray, gridsize: int = 512):
    """Sehr glatte 1D-Dichte über Histogramm + Gauß-Glättung."""
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 2: return None
    
    x_min, x_max = np.min(x), np.max(x)
    r = x_max - x_min if x_max > x_min else 1.0
    pad = 0.15 * r
    edges = np.linspace(x_min - pad, x_max + pad, gridsize + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    hist, _ = np.histogram(x, bins=edges, density=True)
    
    std = np.std(x, ddof=1)
    iqr = np.subtract(*np.percentile(x, [75, 25]))
    sigma_silver = 0.9 * min(std, iqr / 1.34) * (x.size ** (-1/5))
    sigma_min = SMOOTH_BW_FRAC * r
    sigma_x = max(sigma_silver, sigma_min)
    
    dx = edges[1] - edges[0]
    sigma_bins = max(sigma_x / dx, 0.5)
    half_w = int(4 * sigma_bins)
    kernel = np.exp(-0.5 * (np.arange(-half_w, half_w + 1) / sigma_bins) ** 2)
    kernel /= kernel.sum()
    dens = np.convolve(hist, kernel, mode="same")
    
    if SMOOTH_2ND_PASS > 0:
        sigma_bins2 = sigma_bins * float(SMOOTH_2ND_PASS)
        half_w2 = int(3 * sigma_bins2)
        k2 = np.exp(-0.5 * (np.arange(-half_w2, half_w2 + 1) / sigma_bins2) ** 2)
        dens = np.convolve(dens, k2 / k2.sum(), mode="same")
        
    return centers, dens

def kdeplot_filled(ax, data, label: str, color: str, lw: float = 1.6, alpha: float = 0.35):
    """
    Glatte Dichte plotten, am Rand auf y=0 andocken und die Dichtekurve zurückgeben.
    """
    arr = pd.to_numeric(pd.Series(data), errors="coerce").dropna().values
    if arr.size < 2: return None
    
    res = _smooth_density_1d(arr, gridsize=600)
    if res is None: return None
    xs, ys = res
    
    nonzero_indices = np.where(ys > (ys.max() * 0.001))[0]
    if nonzero_indices.size == 0: return xs, ys

    start_idx = max(0, nonzero_indices[0] - 5)
    end_idx = min(len(ys) - 1, nonzero_indices[-1] + 5)
    xs_t = xs[start_idx : end_idx + 1]
    ys_t = ys[start_idx : end_idx + 1]
    if xs_t.size < 2: return xs, ys

    xs_p = np.concatenate(([xs_t[0]], xs_t, [xs_t[-1]]))
    ys_p = np.concatenate(([0.0], ys_t, [0.0]))
    
    min_len = min(len(xs_p), len(ys_p))
    if min_len < 2: return xs, ys
    xs_p, ys_p = xs_p[:min_len], ys_p[:min_len]

    ax.plot(xs_p, ys_p, lw=lw, label=label, color=color)
    ax.fill_between(xs_p, ys_p, 0, alpha=alpha, color=color, linewidth=0)
    y_scatter = np.interp(arr, xs_p, ys_p)
    ax.scatter(arr, y_scatter, color=color, s=15, alpha=0.8, zorder=5, edgecolor='black', linewidth=0.3)
    
    return xs, ys # Gibt die *originale*, ungepaddete Dichtekurve zurück

# =========================
#   SUMMARY & INTERPRET
# =========================
def _aggregate_key_table(df: pd.DataFrame) -> pd.DataFrame:
    """Kompakte Kennzahlen-Tabelle."""
    out = []
    cols_to_agg = [
        ("mean force [pN]", "Force (pN)"),
        ("Delta Lc abs (nm)", "|ΔLc| (nm)"),
        ("Work [pN*nm]", "Work (pN*nm)")
    ]
    for name, sel in [("unfolding", df[df["Direction"] == "unfolding"]),
                      ("refolding", df[df["Direction"] == "refolding"]),
                      ("all", df)]:
        row = {"Set": name, "n": int(len(sel))}
        for col, key in cols_to_agg:
            if col in sel.columns and not sel[col].dropna().empty:
                v = pd.to_numeric(sel[col], errors="coerce").dropna()
                row[f"{key} mean"] = np.mean(v)
                row[f"{key} sd"] = np.std(v, ddof=1) if len(v) > 1 else 0.0
                row[f"{key} median"] = np.median(v)
        out.append(row)
    t = pd.DataFrame(out).dropna(axis=1, how='all')
    for c in t.columns:
        if c not in ("Set", "n"):
            t[c] = t[c].map(lambda x: f"{x:.2f}" if pd.notna(x) else np.nan)
    return t

def _tabulate(df: pd.DataFrame) -> str:
    """Monospace-Tabellentext für PDF-Seite."""
    if df.empty: return "No data."
    cols = list(df.columns)
    widths = [max(len(str(c)), df[c].astype(str).str.len().max()) + 2 for c in cols]
    header = "".join(str(c).ljust(w) for c, w in zip(cols, widths))
    line = "-" * len(header)
    rows = [header, line]
    for _, r in df.iterrows():
        vals = [str(r[c]) if not pd.isna(r[c]) else "—" for c in cols]
        rows.append("".join(v.ljust(w) for v, w in zip(vals, widths)))
    return "\n".join(rows)

def _format_p_value(p):
    if p is None: return "nicht verfügbar"
    return f"p < 0.001" if p < 0.001 else f"p = {p:.3f}"

def _generate_interpretation(df: pd.DataFrame, p_values: dict) -> str:
    """Erzeugt eine detaillierte biophysikalische Interpretation."""
    unf = df[df["Direction"] == "unfolding"].copy()
    ref = df[df["Direction"] == "refolding"].copy()

    unf["mean force [pN]"] = pd.to_numeric(unf["mean force [pN]"], errors='coerce')
    ref["mean force [pN]"] = pd.to_numeric(ref["mean force [pN]"], errors='coerce')
    unf["Delta Lc abs (nm)"] = pd.to_numeric(unf["Delta Lc abs (nm)"], errors='coerce')
    ref["Delta Lc abs (nm)"] = pd.to_numeric(ref["Delta Lc abs (nm)"], errors='coerce')
    unf.dropna(subset=["mean force [pN]", "Delta Lc abs (nm)"], inplace=True)
    ref.dropna(subset=["mean force [pN]", "Delta Lc abs (nm)"], inplace=True)

    num_files = df['Filename'].nunique() if 'Filename' in df.columns else 'N/A'
    
    out = [
        "## Biophysikalische Interpretation der Ergebnisse", "",
        f"Diese Analyse fasst **{len(df)}** Events aus **{num_files}** Experimenten zusammen. "
        "Die folgenden Beobachtungen basieren auf den aggregierten Verteilungen der Entfaltungs- (unfolding) und "
        "Wiedereinfaltungs- (refolding) Zyklen."
    ]

    if not unf.empty and not ref.empty:
        median_force_unf, std_force_unf = unf["mean force [pN]"].median(), unf["mean force [pN]"].std()
        median_force_ref, std_force_ref = ref["mean force [pN]"].median(), ref["mean force [pN]"].std()
        out.append("\n### 1. Analyse der Faltungskräfte")
        out.append(
            f"Die **Entfaltung** erfordert eine mediane Kraft von **{median_force_unf:.1f} ± {std_force_unf:.1f} pN**. "
            f"Die **Wiedereinfaltung** geschieht bei **{median_force_ref:.1f} ± {std_force_ref:.1f} pN**."
        )
        p_force_str = _format_p_value(p_values.get('force'))
        out.append(f"Ein Mann-Whitney-U-Test bestätigt, dass dieser Unterschied statistisch signifikant ist (**{p_force_str}**).")

        if median_force_unf > median_force_ref + std_force_ref:
            hysteresis = median_force_unf - median_force_ref
            out.append(
                f"Die höhere Kraft bei der Entfaltung deutet auf eine **mechanische Hysterese** von ca. **{hysteresis:.1f} pN** hin. "
                "Dies ist typisch für Experimente außerhalb des thermodynamischen Gleichgewichts und spiegelt die Energie wider, "
                "die zur Überwindung der Faltungsenergiebarriere bei einer endlichen Ziehgeschwindigkeit benötigt wird."
            )
        else:
            out.append(
                "Die Kräfte sind sehr ähnlich, was auf einen Prozess nahe am **thermodynamischen Gleichgewicht** hindeutet."
            )
        out.append(
            "Die **Streuung der Kräfte** (siehe Breite der Peaks im PDF-Diagramm) gibt Aufschluss über die Homogenität des Prozesses. "
            "Enge Verteilungen deuten auf einen wohldefinierten Zwei-Zustands-Übergang hin."
        )

        median_dlc_unf, std_dlc_unf = unf["Delta Lc abs (nm)"].median(), unf["Delta Lc abs (nm)"].std()
        median_dlc_ref, std_dlc_ref = ref["Delta Lc abs (nm)"].median(), ref["Delta Lc abs (nm)"].std()
        out.append("\n### 2. Analyse der Konturlängenänderung (|ΔLc|)")
        out.append(
            f"Die mediane Längenänderung beträgt **{median_dlc_unf:.1f} ± {std_dlc_unf:.1f} nm** (Entfaltung) und **{median_dlc_ref:.1f} ± {std_dlc_ref:.1f} nm** (Wiedereinfaltung)."
        )
        
        if abs(median_dlc_unf - median_dlc_ref) < 0.1 * (median_dlc_unf + median_dlc_ref):
             out.append(
                "Die Übereinstimmung der |ΔLc|-Werte bestätigt, dass der Prozess **strukturell reversibel** ist. "
                "Das Molekül kehrt in seinen ursprünglichen Zustand zurück."
            )
        else:
            out.append(
                "Eine signifikante Abweichung der |ΔLc|-Werte könnte auf eine **inkomplette Wiedereinfaltung** oder fehlgefaltete Zustände hindeuten."
            )
        out.append(
            "Das |ΔLc|-Diagramm mit **automatischer Peak-Erkennung** kann auf das Vorhandensein verschiedener, stabiler Entfaltungs-Events hinweisen (z.B. Entfalten mehrerer Domänen)."
        )

    out.extend([
        "\n### 3. Fazit und Ausblick",
        "Zusammenfassend deuten die Daten auf einen **mechanisch [hysteretischen/reversiblen]** Prozess mit einer **konsistenten strukturellen Änderung** hin. "
        "Der Scatter-Plot (Force vs. |ΔLc|) kann Korrelationen zwischen Stabilität (Kraft) und Größe (|ΔLc|) aufzeigen.", "",
        "**Hinweis:** Diese Analyse basiert auf aggregierten Daten. Die Untersuchung einzelner FD-Kurven ist entscheidend, um heterogenes Verhalten zu identifizieren."
    ])

    final_text = "\n".join(out)
    if not unf.empty and not ref.empty:
        placeholder = "[hysteretischen/reversiblen]"
        replacement = "hysteretischen" if unf["mean force [pN]"].median() > ref["mean force [pN]"].median() else "weitgehend reversiblen"
        final_text = final_text.replace(placeholder, replacement)

    return final_text

# =========================
#  STATISTICAL ANALYSIS
# =========================
def _perform_statistical_tests(df: pd.DataFrame) -> dict:
    if mannwhitneyu is None:
        return {}
    
    p_values = {}
    unf = df[df["Direction"] == "unfolding"]
    ref = df[df["Direction"] == "refolding"]

    for key, col in [("force", "mean force [pN]"), ("dlc", "Delta Lc abs (nm)"), ("work", "Work [pN*nm]")]:
        if col in df.columns:
            d1 = unf[col].dropna()
            d2 = ref[col].dropna()
            if len(d1) >= 5 and len(d2) >= 5:
                try:
                    _, p = mannwhitneyu(d1, d2, alternative='two-sided')
                    p_values[key] = p
                except Exception:
                    p_values[key] = None
    return p_values

# =========================
#        REPORT-PDF
# =========================
def _plot_example_curves(pdf: PdfPages, curves: list, num_examples: int = 3):
    """Plottet einige zufällige Beispielkurven auf einer neuen PDF-Seite."""
    unf_curves = [c for c in curves if c[2] == "unfolding"]
    ref_curves = [c for c in curves if c[2] == "refolding"]
    if not unf_curves and not ref_curves:
        return

    fig, axs = plt.subplots(2, num_examples, figsize=(7, 4.5), constrained_layout=True, sharex=True, sharey=True)
    fig.suptitle("Beispielhafte FD-Kurven", fontsize=12, weight='bold')

    for i in range(num_examples):
        # Unfolding examples
        ax = axs[0, i]
        if i < len(unf_curves):
            x, y, _ = random.choice(unf_curves)
            ax.plot(x, y, color=COLOR["unfolding"])
        ax.set_title(f"Unfolding Bsp. #{i+1}")
        ax.grid(alpha=GRID_ALPHA)
        if i == 0: ax.set_ylabel("Force (pN)")

        # Refolding examples
        ax = axs[1, i]
        if i < len(ref_curves):
            x, y, _ = random.choice(ref_curves)
            ax.plot(x, y, color=COLOR["refolding"])
        ax.set_title(f"Refolding Bsp. #{i+1}")
        ax.grid(alpha=GRID_ALPHA)
        if i == 0: ax.set_ylabel("Force (pN)")
        ax.set_xlabel("Distance (nm)")

    pdf.savefig(fig); plt.close(fig)

def generate_batch_pdf(analysis_folder: str | Path, pdf_name: Optional[str] = None) -> str:
    folder = Path(analysis_folder)
    df = _collect_batch_results(folder)
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    if pdf_name is None:
        pdf_name = f"Batch_Report_{datetime.now().strftime('%Y%m%d-%H%M%S')}.pdf"
    out_path = str(folder / pdf_name)
    
    p_values = _perform_statistical_tests(df)
    has_work_data = "Work [pN*nm]" in df.columns and df["Work [pN*nm]"].notna().sum() > 10

    with PdfPages(out_path) as pdf:
        # Seite 1: Dichte-Verteilungen
        num_plots = 3 if has_work_data else 2
        fig, axs = plt.subplots(1, num_plots, figsize=(3.2 * num_plots, 3.6), constrained_layout=True)
        axs = np.atleast_1d(axs)

        # Plot 1: Force
        ax = axs[0]; ax.grid(alpha=GRID_ALPHA); ax.set_title("PDF: Force")
        ref_force = df.loc[df["Direction"] == "refolding", "mean force [pN]"].dropna()
        if len(ref_force) >= 2:
            kdeplot_filled(ax, ref_force, "refolding", COLOR["refolding"])
        unf_force = df.loc[df["Direction"] == "unfolding", "mean force [pN]"].dropna()
        if len(unf_force) >= 2:
            kdeplot_filled(ax, unf_force, "unfolding", COLOR["unfolding"])
        if not df.empty and "mean force [pN]" in df.columns and not df["mean force [pN]"].dropna().empty:
            ax.set_xlim(0, max(30, df["mean force [pN]"].quantile(0.99, interpolation='higher')))
        ax.set_ylim(bottom=0); ax.set_xlabel("Force (pN)"); ax.set_ylabel("Probability"); ax.legend()
        
        # Plot 2: |ΔLc| mit Peak-Erkennung
        ax = axs[1]; ax.grid(alpha=GRID_ALPHA); ax.set_title("PDF: |ΔLc| mit Peaks")
        for direction in ["refolding", "unfolding"]:
            data = df.loc[df["Direction"] == direction, "Delta Lc abs (nm)"].dropna()
            if len(data) < 2: continue
            
            # ### KORRIGIERTE LOGIK FÜR PEAK-ERKENNUNG ###
            density_data = kdeplot_filled(ax, data, direction, COLOR[direction])
            
            if find_peaks and density_data:
                xs, ys = density_data
                peaks, _ = find_peaks(ys, height=ys.max()*0.1, distance=10)
                for p_idx in peaks:
                    if p_idx < len(xs): # Sicherheitsprüfung
                        ax.axvline(xs[p_idx], color=COLOR[direction], ls='--', alpha=0.8, lw=1)
                        ax.text(xs[p_idx], ys[p_idx], f' {xs[p_idx]:.1f}', color=COLOR[direction], va='bottom', ha='center', fontsize=8)

        if not df.empty and "Delta Lc abs (nm)" in df.columns and not df["Delta Lc abs (nm)"].dropna().empty:
            ax.set_xlim(0, 60) # Festes Limit wie gewünscht
        ax.set_ylim(bottom=0); ax.set_xlabel("|ΔLc| (nm)")
        
        # Plot 3: Work (optional)
        if has_work_data:
            ax = axs[2]; ax.grid(alpha=GRID_ALPHA); ax.set_title("PDF: Faltungsarbeit")
            ref_work = df.loc[df["Direction"] == "refolding", "Work [pN*nm]"].dropna()
            if len(ref_work) >= 2:
                kdeplot_filled(ax, ref_work, "refolding", COLOR["refolding"])
            unf_work = df.loc[df["Direction"] == "unfolding", "Work [pN*nm]"].dropna()
            if len(unf_work) >= 2:
                kdeplot_filled(ax, unf_work, "unfolding", COLOR["unfolding"])
            if not df["Work [pN*nm]"].dropna().empty:
                ax.set_xlim(0, df["Work [pN*nm]"].quantile(0.99, interpolation='higher'))
            ax.set_ylim(bottom=0); ax.set_xlabel("Work (pN*nm)")
        pdf.savefig(fig); plt.close(fig)

        # Seiten 2 & 3: Kurven-Plots
        curves = _collect_fd_curves(folder, max_curves=2000)
        if curves:
            # Overlay
            fig = plt.figure(figsize=(6.5, 5.2), constrained_layout=True)
            ax = plt.gca(); ax.grid(alpha=GRID_ALPHA)
            n_ref = sum(1 for _, _, d in curves if d == "refolding")
            n_unf = sum(1 for _, _, d in curves if d == "unfolding")
            for x, y, d in curves:
                if d in COLOR:
                    ax.plot(np.r_[x[0], x, x[-1]], np.r_[0.0, y, 0.0], lw=0.8, alpha=0.4, color=COLOR[d])
            ax.plot([], [], color=COLOR["refolding"], label=f"refolding (n={n_ref})")
            ax.plot([], [], color=COLOR["unfolding"], label=f"unfolding (n={n_unf})")
            ax.set_xlabel("Relative distance (nm)"); ax.set_ylabel("Force (pN)")
            ax.set_title("Overlay aller geglätteten FD-Kurven"); ax.legend(loc="upper left")
            pdf.savefig(fig); plt.close(fig)
            # Beispiele
            _plot_example_curves(pdf, curves)

        # Seite 4: Scatter-Plot
        plot_df = df.dropna(subset=["mean force [pN]", "Delta Lc abs (nm)", "Direction"])
        if not plot_df.empty:
            fig, ax = plt.subplots(figsize=(6.5, 5.2), constrained_layout=True)
            ax.grid(alpha=GRID_ALPHA)
            for direction, group in plot_df.groupby("Direction"):
                if direction in COLOR:
                    ax.scatter(group["Delta Lc abs (nm)"], group["mean force [pN]"],
                               color=COLOR[direction], label=direction, alpha=0.6, s=12)
            ax.set_xlabel("|ΔLc| (nm)"); ax.set_ylabel("Force (pN)")
            ax.set_title("Scatter-Plot: Force vs. |ΔLc|"); ax.legend(title="Direction")
            pdf.savefig(fig); plt.close(fig)

        # Seite 5: Summary-Tabelle
        if not df.empty:
            keytab = _aggregate_key_table(df)
            fig = plt.figure(figsize=(7.5, 5.0))
            ax = plt.gca(); ax.axis("off")
            header = f"Batch: {folder.name}\nGenerated: {ts}\nFiles: {df['Filename'].nunique() if 'Filename' in df.columns else 0}    Rows: {len(df)}"
            ax.text(0.02, 0.95, header, va="top", fontsize=11)
            ax.text(0.02, 0.82, _tabulate(keytab), family="monospace", fontsize=9, va="top")
            fig.tight_layout(pad=1.5) 
            pdf.savefig(fig); plt.close(fig)
        
        # Seite 6: Interpretation
        if not df.empty:
            interpretation_text = _generate_interpretation(df, p_values)
            fig = plt.figure(figsize=(7, 9.0), constrained_layout=True)
            ax = plt.gca(); ax.axis("off")
            ax.text(0.05, 0.95, interpretation_text, va="top", fontsize=10, linespacing=1.5, wrap=True)
            fig.tight_layout(pad=1.5)
            pdf.savefig(fig); plt.close(fig)

    return out_path

# =========================
#   BATCH AGGREGATION
# =========================
def _collect_batch_results(analysis_folder: Path) -> pd.DataFrame:
    parts = []
    for p in analysis_folder.iterdir():
        if not p.is_file() or p.suffix.lower() not in (".csv", ".xlsx"): continue
        df = _read_table(p)
        if df.empty: continue
        cols = {str(c).lower().strip(): c for c in df.columns}
        
        col_file = cols.get("filename") or cols.get("file")
        col_dir = cols.get("direction") or cols.get("orientation")
        col_force = cols.get("mean force [pn]") or cols.get("force") or cols.get("fc")
        col_sl = cols.get("step length [nm]") or cols.get("step length") or cols.get("steplength")
        col_work = cols.get("work [pn*nm]") or cols.get("work [pn nm]") or cols.get("work")

        if not any([col_file, col_dir, col_force, col_sl, col_work]): continue
        
        out = pd.DataFrame()
        out["Filename"] = df[col_file] if col_file else p.stem
        
        d_raw = df[col_dir].astype(str).str.lower() if col_dir else pd.Series([_detect_direction_from_name(p.name)] * len(df))
        out["Direction"] = d_raw.replace({
            "forward": "unfolding", "fw": "unfolding", "unfold": "unfolding",
            "reverse": "refolding", "rv": "refolding", "refold": "refolding"
        }).fillna("unknown")

        if col_force: out["mean force [pN]"] = pd.to_numeric(df[col_force], errors="coerce")
        if col_work: out["Work [pN*nm]"] = pd.to_numeric(df[col_work], errors="coerce")
        if col_sl:
            out["Step length [nm]"] = pd.to_numeric(df[col_sl], errors="coerce")
            sign = np.where(out["Direction"].str.startswith("refold"), -1.0, 1.0)
            out["Delta Lc (nm)"] = out["Step length [nm]"] * sign
            out["Delta Lc abs (nm)"] = out["Delta Lc (nm)"].abs()
        
        parts.append(out)

    if not parts: return pd.DataFrame()
    return pd.concat(parts, ignore_index=True)

# =========================
#   JOB LAUNCH & ROUTES
# =========================
def _run_batch(job_id: str, folder: str):
    try:
        JOBS[job_id]["status"] = "running"
        pdf_path = generate_batch_pdf(folder)
        JOBS[job_id]["pdf_path"] = pdf_path
        JOBS[job_id]["done"] = True
        JOBS[job_id]["status"] = "done"
    except Exception as e:
        JOBS[job_id]["error"] = str(e)
        JOBS[job_id]["status"] = "error"
        JOBS[job_id]["done"] = True

@report_bp.post("/batch/start")
def start_batch():
    data = request.get_json(silent=True) or {}
    folder = data.get("folder")
    if not folder: return jsonify({"error": "Missing 'folder' in JSON body"}), 400
    if not Path(folder).exists(): return jsonify({"error": f"Folder not found: {folder}"}), 404
    job_id = uuid.uuid4().hex[:12]
    JOBS[job_id] = {
        "id": job_id, "folder": folder, "created": datetime.now().isoformat(timespec="seconds"),
        "status": "queued", "done": False, "error": None, "pdf_path": None,
    }
    threading.Thread(target=_run_batch, args=(job_id, folder), daemon=True).start()
    return jsonify({"job_id": job_id, "status": "queued"})

@report_bp.get("/job/<job_id>")
def job_status(job_id):
    job = JOBS.get(job_id)
    if not job: return jsonify({"error": "job not found"}), 404
    resp = {k: v for k, v in job.items() if k != "pdf_path"}
    resp["report_ready"] = bool(job.get("done") and job.get("pdf_path"))
    return jsonify(resp)

@report_bp.get("/job/<job_id>/report.pdf")
def job_report_pdf(job_id):
    job = JOBS.get(job_id)
    if not job: abort(404)
    if not job.get("done") or not job.get("pdf_path"):
        return jsonify({"error": "report not ready"}), 400
    pdf_path = Path(job["pdf_path"])
    if not pdf_path.exists():
        return jsonify({"error": "report missing on disk"}), 410
    return send_file(pdf_path, as_attachment=True)