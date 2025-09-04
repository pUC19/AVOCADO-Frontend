# batch_report.py — AVOCADO Batch Report (extra-smooth densities, fixed axes, axis-anchored)
from __future__ import annotations

import re
import uuid
import threading
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, Optional, Tuple, List

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from flask import Blueprint, request, jsonify, abort, send_file

# ---- optional SciPy (nur für Savitzky-Golay + KDE-Peaks; Report läuft auch ohne) ----
try:
    from scipy.signal import savgol_filter, find_peaks
except Exception:
    savgol_filter = None
    find_peaks = None
try:
    from scipy.stats import gaussian_kde
except Exception:
    gaussian_kde = None

report_bp = Blueprint("batch_report", __name__, url_prefix="/report")
JOBS: Dict[str, Dict[str, Any]] = {}

# Farben/Style
COLOR = {"unfolding": "#f39c12", "refolding": "#1f77b4"}  # orange / blau
GRID_ALPHA = 0.12
plt.rcParams.update({
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.titleweight": "bold",
    "axes.labelsize": 10,
    "axes.titlesize": 11,
    "legend.frameon": False,
    "figure.dpi": 300  # Erhöhte DPI für schärfere Raster-Elemente
})

# --- KDE-Glättung (Paper-like) ---
SMOOTH_BW_FRAC = 0.08   # 8 % der x-Spannweite als Mindest-Bandbreite (glatt!)
SMOOTH_2ND_PASS = 1.5   # zweiter sanfter Glättungs-Pass (1.5x breiter). 0 = aus.

# =========================
#         IO / READ
# =========================
def _read_table(path: Path) -> pd.DataFrame:
    """CSV/XLSX robust lesen. '>‘-Kommentarzeile ok. Headerlos ok."""
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
            first = fh.readline()
        if first.startswith(">"):
            skip = 1
    except Exception:
        pass

    for kw in ({"sep": ",", "skiprows": skip},
               {"sep": ";", "decimal": ",", "skiprows": skip}):
        try:
            df = pd.read_csv(path, **kw)
            if isinstance(df, pd.DataFrame) and not df.empty:
                return df
        except Exception:
            pass

    for kw in ({"sep": ",", "header": None, "skiprows": skip},
               {"sep": ";", "decimal": ",", "header": None, "skiprows": skip}):
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
    return ("smooth" in n) or n.endswith("_smooth.csv") or n.endswith("-smooth.csv")


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
    """X=Distanz (nm), Y=Force (pN) erkennen (headerlos/smooth inkl. µm→nm)."""
    exact_x = [c for c in df.columns if str(c).strip().lower() == "relative distance, nm"]
    exact_y = [c for c in df.columns if str(c).strip().lower() == "force, pn"]
    if exact_x and exact_y:
        x = pd.to_numeric(df[exact_x[0]], errors="coerce").values
        y = pd.to_numeric(df[exact_y[0]], errors="coerce").values
    elif df.shape[1] == 2:
        a = pd.to_numeric(df.iloc[:, 0], errors="coerce").values
        b = pd.to_numeric(df.iloc[:, 1], errors="coerce").values
        cand: List[Tuple[np.ndarray, np.ndarray]] = []
        if _is_force_like(a) and (_is_dist_nm_like(b) or _is_dist_um_like(b)):
            x1 = b * 1000.0 if _is_dist_um_like(b) else b
            cand.append((x1, a))
        if _is_force_like(b) and (_is_dist_nm_like(a) or _is_dist_um_like(a)):
            x2 = a * 1000.0 if _is_dist_um_like(a) else a
            cand.append((x2, b))
        if cand:
            x, y = cand[0]
        else:
            # Fallback: Spalte mit größerer Range = Distanz
            def rng(u):
                u = u[np.isfinite(u)]
                return np.percentile(u, 99) - np.percentile(u, 1) if u.size else 0.0
            if rng(a) >= rng(b):
                xa, yb = a, b
            else:
                xa, yb = b, a
            xa = xa * 1000.0 if _is_dist_um_like(xa) and not _is_dist_nm_like(xa) else xa
            x, y = xa, yb
    else:
        num = [c for c in df.columns if pd.to_numeric(df[c], errors="coerce").notna().sum() >= 5]
        if len(num) < 2: return None
        x = pd.to_numeric(df[num[0]], errors="coerce").values
        y = pd.to_numeric(df[num[1]], errors="coerce").values

    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 10: return None
    x, y = x[m], y[m]
    if np.nanmedian(np.diff(x)) < 0:
        x, y = x[::-1], y[::-1]
    if np.nanmedian(y) < 0 and (np.nanmax(y) <= 0 or abs(np.nanmedian(y)) > abs(np.nanmax(y))*0.5):
        y = -y
    return x, y


# =========================
#   CURVE CLEAN & SMOOTH
# =========================
def _robust_clip(y: np.ndarray, q=(0.005, 0.995)) -> np.ndarray:
    lo, hi = np.nanquantile(y, q)
    return np.clip(y, lo, hi)

def _prepare_curve(x: np.ndarray, y: np.ndarray, target_n: int = 1400) -> Tuple[np.ndarray, np.ndarray]:
    """Sortieren, Duplikate raus, clippen, glätten, resamplen → schöne Overlay-Linien."""
    x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if x.size < 10: return x, y

    order = np.argsort(x); x, y = x[order], y[order]
    ux, idx = np.unique(x, return_index=True); x, y = ux, y[idx]
    if x.size < 10: return x, y

    y = _robust_clip(y)
    if savgol_filter is not None and x.size >= 31:
        try:
            win = max(41, int(x.size * 0.04) // 2 * 2 + 1)
            y = savgol_filter(y, window_length=win, polyorder=3, mode="interp")
        except Exception:
            pass
    else:
        k = max(7, int(x.size * 0.02) // 2 * 2 + 1)
        if k > 1:
            pad = k // 2
            y_pad = np.pad(y, (pad, pad), mode="edge")
            ker = np.ones(k) / k
            y = np.convolve(y_pad, ker, mode="valid")

    n = min(target_n, max(200, x.size))
    xg = np.linspace(x[0], x[-1], n)
    yg = np.interp(xg, x, y)
    return xg, yg


def _collect_fd_curves(analysis_folder: str | Path, max_curves: int = 2000) -> List[Tuple[np.ndarray, np.ndarray, str]]:
    folder = Path(analysis_folder)
    curves: List[Tuple[np.ndarray, np.ndarray, str]] = []
    files = [p for p in folder.iterdir() if p.is_file() and p.suffix.lower() in (".csv", ".xlsx")]
    smooth = [p for p in files if _is_smooth_filename(p)]
    use = smooth if smooth else files
    use.sort(key=lambda q: (("fd curve" not in q.name.lower()), q.name.lower()))

    for p in use:
        if len(curves) >= max_curves: break
        df = _read_table(p)
        if df is None or df.empty: continue
        xy = _read_fd_xy(df)
        if xy is None: continue
        x, y = _prepare_curve(*xy)
        if x.size < 10: continue
        d = _curve_direction_from_df_or_name(df, p)
        curves.append((x, y, d))
    return curves


# =========================
#      STABLE DENSITIES (EXTRA-SMOOTH)
# =========================
def _smooth_density_1d(values: np.ndarray, gridsize: int = 512):
    """
    Sehr glatte 1D-Dichte über Histogramm + Gauß-Glättung.
    IMMER (N, N) mit N=gridsize. Deutlich erhöhte Mindest-Bandbreite.
    """
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    n = x.size
    if n < 2:
        return None

    x_min, x_max = np.min(x), np.max(x)
    r = x_max - x_min if x_max > x_min else max(abs(x_max), 1.0)
    pad = 0.05 * r
    edges = np.linspace(x_min - pad, x_max + pad, gridsize + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])  # Länge N

    hist, _ = np.histogram(x, bins=edges, density=True)  # Länge N

    std = np.std(x, ddof=1) if n > 1 else 0.0
    iqr = np.subtract(*np.percentile(x, [75, 25])) if n > 1 else 0.0
    sigma_silver = 0.9 * min(std, iqr / 1.34) * (n ** (-1/5)) if (std > 0 or iqr > 0) else 0.0
    sigma_min = max(SMOOTH_BW_FRAC * r, 1e-6)           # breite Mindest-BW -> keine Stacheln
    sigma_x = max(sigma_silver, sigma_min)

    dx = edges[1] - edges[0]
    sigma_bins = max(sigma_x / dx, 0.5)
    half_w = int(4 * sigma_bins)
    g = np.arange(-half_w, half_w + 1, dtype=float)
    kernel = np.exp(-0.5 * (g / sigma_bins) ** 2)
    kernel /= kernel.sum()

    dens = np.convolve(hist, kernel, mode="same")

    if SMOOTH_2ND_PASS and SMOOTH_2ND_PASS > 0:
        sigma_bins2 = sigma_bins * float(SMOOTH_2ND_PASS)
        half_w2 = int(3 * sigma_bins2)
        g2 = np.arange(-half_w2, half_w2 + 1, dtype=float)
        k2 = np.exp(-0.5 * (g2 / sigma_bins2) ** 2)
        k2 /= k2.sum()
        dens = np.convolve(dens, k2, mode="same")

    if dens.size != centers.size:
        x_src = np.linspace(centers[0], centers[-1], dens.size)
        dens = np.interp(centers, x_src, dens)

    return centers, dens


def kdeplot_filled(
    ax,
    data,
    label: str,
    color: str,
    lw: float = 1.6,
    alpha: float = 0.35,
):
    """Glatte Dichte (Histogramm + Gauß). An den Rändern der Kurve auf y=0 andocken."""
    arr = pd.to_numeric(pd.Series(data), errors="coerce").dropna().values
    if arr.size < 2:
        return
    res = _smooth_density_1d(arr, gridsize=600)
    if res is None:
        return
    xs, ys = res
    if xs.size != ys.size:
        n = min(xs.size, ys.size)
        xs, ys = xs[:n], ys[:n]
        
    # --- Kurve am Anfang und Ende auf der x-Achse aufsetzen lassen ---
    valid_indices = np.where(ys > 1e-6)[0]
    if valid_indices.size == 0:
        return
    
    start_idx = valid_indices[0]
    end_idx = valid_indices[-1]

    xs_trimmed = xs[start_idx:end_idx+1]
    ys_trimmed = ys[start_idx:end_idx+1]
    
    xs_padded = np.r_[xs_trimmed[0], xs_trimmed, xs_trimmed[-1]]
    ys_padded = np.r_[0.0, ys_trimmed, 0.0]

    # Bisheriger Code: Zeichnet die Linie und die Füllung
    ax.plot(xs_padded, ys_padded, lw=lw, label=label, color=color)
    ax.fill_between(xs_padded, ys_padded, 0, alpha=alpha, color=color, linewidth=0)

    # --- NEU: Originale Datenpunkte auf die Kurve zeichnen ---
    # Wir interpolieren die y-Position für jeden originalen Datenpunkt aus der geglätteten Kurve.
    # Dies sorgt dafür, dass die Punkte exakt auf der Linie liegen.
    y_scatter = np.interp(arr, xs_padded, ys_padded)
    
    # Scatter-Plot der Punkte hinzufügen
    ax.scatter(arr, y_scatter, color=color, s=15, alpha=0.8, zorder=5, edgecolor='black', linewidth=0.3)


# =========================
#   SUMMARY & INTERPRET
# =========================
def _aggregate_key_table(df: pd.DataFrame) -> pd.DataFrame:
    """Kompakte Kennzahlen-Tabelle pro Richtung + Gesamt."""
    out = []
    for name, sel in [("unfolding", df[df["Direction"] == "unfolding"]),
                      ("refolding", df[df["Direction"] == "refolding"]),
                      ("all", df)]:
        row = {"Set": name, "n": int(len(sel))}
        for col, key in [("mean force [pN]", "Force (pN)"), ("Delta Lc abs (nm)", "|ΔLc| (nm)")]:
            if col in sel.columns and not sel[col].dropna().empty:
                v = pd.to_numeric(sel[col], errors="coerce").dropna()
                row[f"{key} mean"] = np.mean(v)
                row[f"{key} sd"] = np.std(v, ddof=1) if len(v) > 1 else 0.0
                row[f"{key} median"] = np.median(v)
            else:
                row[f"{key} mean"] = np.nan
                row[f"{key} sd"] = np.nan
                row[f"{key} median"] = np.nan
        out.append(row)
    t = pd.DataFrame(out)
    for c in t.columns:
        if c not in ("Set", "n"):
            t[c] = t[c].map(lambda x: np.nan if pd.isna(x) else float(f"{x:.2f}"))
    return t

def _tabulate(df: pd.DataFrame) -> str:
    """Monospace-Tabellentext für PDF-Seite."""
    if df.empty: return "No data."
    cols = ["Set", "n", "Force (pN) mean", "Force (pN) sd", "Force (pN) median",
            "|ΔLc| (nm) mean", "|ΔLc| (nm) sd", "|ΔLc| (nm) median"]
    df = df[cols]
    widths = [10, 6, 14, 12, 16, 16, 12, 18]
    header = " ".join(str(c).ljust(w) for c, w in zip(cols, widths))
    line = "-" * (sum(widths) + len(widths) - 1)
    rows = [header, line]
    for _, r in df.iterrows():
        vals = [str(r[c]) if not pd.isna(r[c]) else "—" for c in cols]
        rows.append(" ".join(v.ljust(w) for v, w in zip(vals, widths)))
    return "\n".join(rows)

def _generate_interpretation(df: pd.DataFrame) -> str:
    unf = df[df["Direction"] == "unfolding"]
    ref = df[df["Direction"] == "refolding"]

    out = [
        "## Interpretation der Stapelanalyse",
        "",
        "Diese Analyse fasst die Ergebnisse mehrerer Einzelmessungen zusammen. Die folgenden Beobachtungen basieren auf den aggregierten Daten."
    ]

    # Force
    if not unf.empty and not ref.empty:
        mean_force_unf = unf["mean force [pN]"].mean()
        mean_force_ref = ref["mean force [pN]"].mean()
        if not np.isnan(mean_force_unf) and not np.isnan(mean_force_ref):
            out.append("### Kräfte (Force)")
            out.append(f"Die mittlere Entfaltungskraft (`unfolding`) liegt bei **{mean_force_unf:.2f} pN**, während die mittlere Wiedereinfaltungskraft (`refolding`) bei **{mean_force_ref:.2f} pN** liegt.")
            if mean_force_unf > mean_force_ref:
                out.append("Das legt nahe, dass zum Entfalten der Proteine eine höhere Kraft benötigt wird als für die Wiedereinfaltung. Dieser Unterschied ist oft ein Anzeichen für eine **hysteretische Schleife**, die durch die molekulare Reibung und die endliche Geschwindigkeit des Experiments verursacht wird.")
            else:
                out.append("Die gemessene Kraft ist bei beiden Prozessen ähnlich, was auf ein weitgehend reversibles System hindeuten könnte.")
    
    # Delta Lc
    if not unf.empty and not ref.empty:
        mean_dlc_unf = unf["Delta Lc abs (nm)"].mean()
        mean_dlc_ref = ref["Delta Lc abs (nm)"].mean()
        if not np.isnan(mean_dlc_unf) and not np.isnan(mean_dlc_ref):
            out.append("### Längenänderung (|ΔLc|)")
            out.append(f"Die durchschnittliche absolute Längenänderung bei der Entfaltung beträgt **{mean_dlc_unf:.2f} nm**, und bei der Wiedereinfaltung **{mean_dlc_ref:.2f} nm**.")
            if abs(mean_dlc_unf - mean_dlc_ref) < 5:  # Schwellenwert für Ähnlichkeit
                out.append("Die Längenänderungen bei Entfaltung und Wiedereinfaltung sind sehr ähnlich, was auf eine konsistente Strukturänderung im Experiment hindeutet. Das Fehlen größerer Abweichungen bestätigt, dass die Proteine in beiden Richtungen den gleichen molekularen Pfad durchlaufen.")
            else:
                out.append("Die unterschiedlichen Längenänderungen könnten auf eine ineffiziente Wiedereinfaltung oder auf eine unterschiedliche Anzahl von Entfaltungsereignissen pro Richtung hinweisen.")

    out.append("---")
    out.append("Dieser Bericht dient der Übersicht. Eine detaillierte Betrachtung der Einzelfälle kann weitere Erkenntnisse liefern.")

    return "\n".join(out)


# =========================
#        REPORT-PDF
# =========================
FORCE_XLIM = (0, 30)   # fix
DLCA_XLIM  = (0, 100)  # fix

def generate_batch_pdf(analysis_folder: str | Path, pdf_name: Optional[str] = None) -> str:
    folder = Path(analysis_folder)
    df = _collect_batch_results(folder)
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    if pdf_name is None:
        pdf_name = f"Batch_Report_{datetime.now().strftime('%Y%m%d-%H%M%S')}.pdf"
    out_path = str(folder / pdf_name)

    if "Delta Lc abs (nm)" in df.columns:
        df["Delta Lc abs (nm)"] = pd.to_numeric(df["Delta Lc abs (nm)"], errors="coerce")
    if "mean force [pN]" in df.columns:
        df["mean force [pN]"] = pd.to_numeric(df["mean force [pN]"], errors="coerce")
    
    # Filtern der Daten für die Plots
    plot_df = df.dropna(subset=["mean force [pN]", "Delta Lc abs (nm)", "Direction"])

    # matplotlib.backends.backend_pdf speichert nativ Vektorgrafiken
    with PdfPages(out_path) as pdf:
        # 1) kompakte Force/ΔLc-Dichten (nebeneinander)
        if not df.empty:
            fig, axs = plt.subplots(1, 2, figsize=(6.5, 3.6), constrained_layout=True)
            # Force
            ax = axs[0]; ax.grid(alpha=GRID_ALPHA)
            refF = df.loc[df["Direction"] == "refolding", "mean force [pN]"]
            unfF = df.loc[df["Direction"] == "unfolding", "mean force [pN]"]
            kdeplot_filled(ax, refF, "refolding", COLOR["refolding"])
            kdeplot_filled(ax, unfF, "unfolding", COLOR["unfolding"])
            ax.set_xlim(*FORCE_XLIM); ax.set_ylim(bottom=0)
            ax.set_xlabel("Force (pN)"); ax.set_ylabel("Wahrscheinlichkeit"); ax.set_title("PDF: Force")
            ax.legend()

            # |ΔLc|
            ax = axs[1]; ax.grid(alpha=GRID_ALPHA)
            refX = df.loc[df["Direction"] == "refolding", "Delta Lc abs (nm)"]
            unfX = df.loc[df["Direction"] == "unfolding", "Delta Lc abs (nm)"]
            kdeplot_filled(ax, refX, "refolding", COLOR["refolding"])
            kdeplot_filled(ax, unfX, "unfolding", COLOR["unfolding"])
            ax.set_xlim(*DLCA_XLIM); ax.set_ylim(bottom=0)
            ax.set_xlabel("|ΔLc| (nm)"); ax.set_title("PDF: |ΔLc|")
            pdf.savefig(fig); plt.close(fig)

        # 2) Overlay FD – kompakte Seite (Start/Ende an x-Achse andocken)
        curves = _collect_fd_curves(folder, max_curves=2000)
        if curves:
            fig = plt.figure(figsize=(6.5, 5.2), constrained_layout=True)
            ax = plt.gca(); ax.grid(alpha=GRID_ALPHA)
            for dir_key in ("refolding", "unfolding"):
                for x, y, d in curves:
                    if d != dir_key: continue
                    # an die x-Achse andocken: y=0 an Anfang/Ende
                    x_plot = np.r_[x[0], x, x[-1]]
                    y_plot = np.r_[0.0,  y, 0.0]
                    ax.plot(x_plot, y_plot, lw=1.0, alpha=0.45, color=COLOR[dir_key])
            n_ref = sum(1 for _, _, d in curves if d == "refolding")
            n_unf = sum(1 for _, _, d in curves if d == "unfolding")
            ax.plot([], [], color=COLOR["refolding"], label=f"refolding (n={n_ref})")
            ax.plot([], [], color=COLOR["unfolding"], label=f"unfolding (n={n_unf})")
            ax.set_xlabel("Relative distance (nm)"); ax.set_ylabel("Force (pN)")
            ax.set_title("Overlay smoothed FD curves"); ax.legend(loc="upper left")
            pdf.savefig(fig); plt.close(fig)

        # 3) Scatter-Plots
        if not plot_df.empty:
            fig, ax = plt.subplots(figsize=(6.5, 5.2), constrained_layout=True)
            ax.grid(alpha=GRID_ALPHA)
            for direction, group in plot_df.groupby("Direction"):
                if direction in COLOR:
                    ax.scatter(group["Delta Lc abs (nm)"], group["mean force [pN]"],
                               color=COLOR[direction], label=direction, alpha=0.6, s=12)
            ax.set_xlabel("|ΔLc| (nm)"); ax.set_ylabel("Force (pN)")
            ax.set_title("Scatter-Plot: Force vs. |ΔLc|")
            ax.legend(title="Direction")
            pdf.savefig(fig); plt.close(fig)

        # 4) Summary-Tabelle (kompakt)
        keytab = _aggregate_key_table(df) if not df.empty else pd.DataFrame()
        fig = plt.figure(figsize=(7, 5.0)) # Etwas breiter gemacht für mehr Platz
        ax = plt.gca(); ax.axis("off")
        header = f"Batch: {folder.name}\nGenerated: {ts}\nFiles: {df['Filename'].nunique() if 'Filename' in df.columns else 0}    Rows: {len(df)}"
        ax.text(0.05, 0.95, header, va="top", fontsize=11) # Linken Rand etwas vergrößert
        ax.text(0.05, 0.86, _tabulate(keytab), family="monospace", fontsize=9, va="top") # Linken Rand etwas vergrößert
        ax.text(0.05, 0.10,
                "Hinweis: Dichten via Histogramm+Gauß (breite Mindest-Bandbreite). KDE & FD an den Achsenenden auf y=0 angedockt.",
                fontsize=8, alpha=0.7)
        
        # NEU: Diese Zeile sorgt für eine automatische Anpassung des Layouts, um Abschneiden zu verhindern
        fig.tight_layout(pad=1.5) 
        
        pdf.savefig(fig); plt.close(fig)
        
        # 5) Interpretation (text-only)
        if not df.empty:
            interpretation_text = _generate_interpretation(df)
            fig = plt.figure(figsize=(7, 5.0), constrained_layout=True)
            ax = plt.gca(); ax.axis("off")
            # Text mit automatischem Umbruch für bessere Lesbarkeit
            ax.text(0.05, 0.95, interpretation_text, va="top", fontsize=10, linespacing=1.5, wrap=True)

            # NEU: Auch hier zur Sicherheit hinzufügen, um das Layout anzupassen
            fig.tight_layout(pad=1.5)

            pdf.savefig(fig); plt.close(fig)

    return out_path


# =========================
#   BATCH AGGREGATION
# =========================
def _collect_batch_results(analysis_folder: Path) -> pd.DataFrame:
    parts = []
    for p in analysis_folder.iterdir():
        if not p.is_file() or p.suffix.lower() not in (".csv", ".xlsx"):
            continue
        df = _read_table(p)
        if df is None or df.empty:
            continue

        cols = {str(c).lower(): c for c in df.columns}
        col_file  = cols.get("filename") or cols.get("file")
        col_dir   = cols.get("direction") or cols.get("orientation")
        col_force = cols.get("mean force [pn]") or cols.get("force") or cols.get("fc")
        col_sl    = cols.get("step length [nm]") or cols.get("step length") \
                    or cols.get("steplength") or cols.get("step_length")

        if not (col_file or col_dir or col_force or col_sl):
            continue

        out = pd.DataFrame()
        out["Filename"] = (df[col_file].astype(str) if (col_file in df) else p.stem)

        if col_dir in df:
            d_raw = df[col_dir].astype(str).str.lower()
        else:
            d_raw = pd.Series([None] * len(out))
        d_norm = d_raw.replace({
            "forward":"unfolding","forwards":"unfolding","fw":"unfolding","unfold":"unfolding",
            "reverse":"refolding","backward":"refolding","rv":"refolding","refold":"refolding"
        }).fillna("unknown")
        d_norm = d_norm.mask(d_norm.eq("unknown"), _detect_direction_from_name(p.name))
        out["Direction"] = d_norm

        if col_force in df:
            out["mean force [pN]"] = pd.to_numeric(df[col_force], errors="coerce")
        if col_sl in df:
            out["Step length [nm]"] = pd.to_numeric(df[col_sl], errors="coerce")

        if "Step length [nm]" in out.columns:
            sign = np.where(out["Direction"].astype(str).str.startswith("refold"), -1.0, 1.0)
            dlc = out["Step length [nm]"].apply(lambda z: float(z) if pd.notna(z) else np.nan)
            out["Delta Lc (nm)"] = dlc * sign
            out["Delta Lc abs (nm)"] = pd.to_numeric(out["Delta Lc (nm)"], errors="coerce").abs()

        parts.append(out)

    if not parts:
        return pd.DataFrame(columns=["Filename","Direction","mean force [pN]","Step length [nm]","Delta Lc (nm)","Delta Lc abs (nm)"])
    df_all = pd.concat(parts, ignore_index=True)
    df_all["Direction"] = df_all["Direction"].fillna("unknown").str.lower().replace({
        "unf":"unfolding","unfold":"unfolding","forward":"unfolding","fw":"unfolding",
        "ref":"refolding","refold":"refolding","reverse":"refolding","rv":"refolding"
    })
    return df_all


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
    if not folder:
        return jsonify({"error": "Missing 'folder' in JSON body"}), 400
    if not Path(folder).exists():
        return jsonify({"error": "Folder not found: {folder}"}), 404

    job_id = uuid.uuid4().hex[:12]
    JOBS[job_id] = {
        "id": job_id, "folder": folder,
        "created": datetime.now().isoformat(timespec="seconds"),
        "status": "queued", "done": False, "error": None, "pdf_path": None,
    }
    threading.Thread(target=_run_batch, args=(job_id, folder), daemon=True).start()
    return jsonify({"job_id": job_id, "status": "queued"})


@report_bp.get("/job/<job_id>")
def job_status(job_id):
    job = JOBS.get(job_id)
    if not job:
        return jsonify({"error": "job not found"}), 404
    resp = {k: v for k, v in job.items() if k != "pdf_path"}
    resp["report_ready"] = bool(job.get("done") and job.get("pdf_path"))
    return jsonify(resp)


@report_bp.get("/job/<job_id>/report.pdf")
def job_report_pdf(job_id):
    job = JOBS.get(job_id)
    if not job:
        abort(404)
    if job.get("status") != "done" or not job.get("pdf_path"):
        return jsonify({"error": "report not ready"}), 400
    pdf_path = job["pdf_path"]
    if not Path(pdf_path).exists():
        return jsonify({"error": "report missing on disk"}), 410
    return send_file(pdf_path, as_attachment=True)
