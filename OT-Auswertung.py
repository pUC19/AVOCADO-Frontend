import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tkinter import Tk, filedialog
from openpyxl import Workbook
from openpyxl.drawing.image import Image as XLImage
from openpyxl.utils.dataframe import dataframe_to_rows

# 1) Datei wählen
root = Tk(); root.withdraw()
xlsx_path = filedialog.askopenfilename(
    title="OT xlsx-Datei auswählen",
    filetypes=[("Excel files", "*.xlsx")]
)
if not xlsx_path:
    print("Abgebrochen."); raise SystemExit(0)
out_path = os.path.splitext(xlsx_path)[0] + "_OT-Auswertung_DeltaLC.xlsx"

# 2) Einlesen
df = pd.read_excel(xlsx_path)
if 'Filename' in df.columns:
    df = df[df['Filename'].notna()].copy()
else:
    df = df.copy()

# 3) Direction robust bestimmen
def detect_direction(fname):
    s = str(fname).lower()
    if any(tag in s for tag in ['1f', 'fw', 'forward', 'unfolding']): return 'unfolding'
    if any(tag in s for tag in ['1r', 'rv', 'reverse', 'refolding']): return 'refolding'
    return 'unknown'

if 'Direction' not in df.columns:
    if 'Filename' in df.columns:
        fname_series = df['Filename'].astype(str)
    else:
        fname_series = pd.Series([None]*len(df), index=df.index, dtype='object')
    df['Direction'] = fname_series.apply(detect_direction)

# 4) ΔLc berechnen
lc_col = 'ss contour Length'
step_col = 'step number' if 'step number' in df.columns else ('step' if 'step' in df.columns else None)

def calc_delta_lc(g):
    gg = g
    if step_col and step_col in gg.columns:
        gg = gg.sort_values(step_col)
    if lc_col not in gg.columns:
        return pd.Series([np.nan]*len(gg), index=gg.index)
    lc = pd.to_numeric(gg[lc_col], errors='coerce')
    delta = lc.diff().abs()
    if len(lc) > 0:
        first = lc.first_valid_index()
        if first is not None:
            delta.loc[first] = lc.loc[first]
    return delta

if lc_col in df.columns:
    df['Delta Lc (nm)'] = (
        df.groupby('Filename', dropna=False, sort=False)
          .apply(calc_delta_lc)
          .reset_index(level=0, drop=True)
    )
else:
    df['Delta Lc (nm)'] = np.nan

# 5) Plots
force_col = "mean force [pN]"
delta_lc_col = "Delta Lc (nm)"
step_len_col = "Step length [nm]"
plotfiles = []

def _save(fig, name):
    fig.tight_layout()
    fig.savefig(name, bbox_inches="tight", dpi=200)
    plt.close(fig)
    plotfiles.append(name)

# 5.1 KDE Force
fig = plt.figure(figsize=(6,4), dpi=120)
for direction, label in zip(['refolding','unfolding'], ['refolding','unfolding']):
    data = pd.to_numeric(df.loc[df['Direction']==direction, force_col], errors='coerce').dropna()
    if not data.empty:
        sns.kdeplot(data=data, fill=True, label=label, lw=2, alpha=0.5)
plt.xlabel("Force (pN)"); plt.ylabel("Wahrscheinlichkeitsdichte"); plt.title("PDF: Force (pN)"); plt.legend()
 _save(fig, "plot_pdf_force.png")

# 5.2 KDE ΔLc
fig = plt.figure(figsize=(6,4), dpi=120)
for direction, label in zip(['refolding','unfolding'], ['refolding','unfolding']):
    data = pd.to_numeric(df.loc[df['Direction']==direction, delta_lc_col], errors='coerce').dropna()
    if not data.empty:
        sns.kdeplot(data=data, fill=True, label=label, lw=2, alpha=0.5)
plt.xlabel("ΔLc (nm)"); plt.ylabel("Wahrscheinlichkeitsdichte"); plt.title("PDF: ΔLc (nm)"); plt.legend()
 _save(fig, "plot_pdf_deltalc.png")

# 5.3 Scatter: Fc ~ ΔLc
fig = plt.figure(figsize=(8,5), dpi=120)
for direction, label in zip(['unfolding','refolding'], ['Unfolding','Refolding']):
    sel = df['Direction']==direction
    x = pd.to_numeric(df.loc[sel, delta_lc_col], errors='coerce')
    y = pd.to_numeric(df.loc[sel, force_col], errors='coerce')
    plt.scatter(x, y, label=label, alpha=0.7, s=20)
plt.xlabel("Delta Lc, nm"); plt.ylabel("Force, pN"); plt.title("Unfolding / Refolding events - Fc ~ ΔLc"); plt.legend()
 _save(fig, "plot_scatter_fc_deltalc.png")

# 5.4 Scatter: Fc ~ Step length (optional)
if step_len_col in df.columns:
    fig = plt.figure(figsize=(8,5), dpi=120)
    for direction, label in zip(['unfolding','refolding'], ['Unfolding','Refolding']):
        sel = df['Direction']==direction
        x = pd.to_numeric(df.loc[sel, step_len_col], errors='coerce')
        y = pd.to_numeric(df.loc[sel, force_col], errors='coerce')
        plt.scatter(x, y, label=label, alpha=0.7, s=20)
    plt.xlabel("Step length, nm"); plt.ylabel("Force, pN"); plt.title("Unfolding / Refolding events - Fc ~ SL"); plt.legend()
     _save(fig, "plot_scatter_fc_sl.png")

# 6) Neues Excel mit Daten + Plots
wb = Workbook()
ws_data = wb.active; ws_data.title = "Data (+ΔLc)"
for r in dataframe_to_rows(df, index=False, header=True):
    ws_data.append(r)

ws_plots = wb.create_sheet("Plots")
row = 1
for pf in plotfiles:
    if os.path.exists(pf):
        img = XLImage(pf); img.width = int(img.width*0.8); img.height = int(img.height*0.8)
        ws_plots.add_image(img, f"A{row}")
        row += 25

wb.save(out_path)
print("Fertig! Gespeichert unter:", out_path)