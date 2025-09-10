# POTATO Frontend (Avocado · Tomato · Carrot)

Flask-Webfrontend zur Analyse von Optical Tweezers Daten mit POTATO. Erstellt von Tobias Fremter/Biochemie3 Uni Regensburg 

## Module

### 🥑 Avocado (Batch)
    (docs/startseite.png)
-   Batch-Analyse von *.h5 / *.csv\* direkt aus Server-Ordnern
-   Multiprocessing-Jobhandling mit Fortschritts-Logs
-   Automatische FD-Plots und Exporte (SMOOTH, PLOT, STEPS, TOTAL, FIT)
-   Download einzelner Ergebnisdateien oder als ZIP

### 🍅 Tomato (Interaktiv)

-   Einzeldateien (Upload oder Serverpfad) laden
-   Interaktive Step-Erkennung und manuelles Setzen
-   DS-Fit und SS-Fits zwischen Steps
-   Berechnung der Arbeit (pN·nm, kT)
-   Export ins standardisierte OT-Schema (CSV + XLSX im ZIP)

### 🥕 Carrot (Excel Postprocessing)

-   Bestehende POTATO-Ergebnis-Excel laden
-   Berechnung von ΔLc via eFJC-Modell (physikalisch)
-   Plots:
    -   **KDE Force**
    -   **KDE \|ΔLc\|**
    -   **Scatter Fc vs \|ΔLc\|** (Unfolding/Refolding, farblich
        getrennt)
    -   **Scatter Fc vs Step length** (falls Spalte vorhanden)
-   Ergebnis: neue Excel mit Datenblatt und eingebetteten Plots

## Beispielplot

![Scatter Beispiel](docs/plot_hexbin_example.png)

## Installation

``` bash
git clone https://github.com/pUC19/AVOCADO.git
cd AVOCADO
pip install -r requirements.txt
python app.py
```

## Anforderungen

-   Python \>= 3.9
-   Flask, NumPy, Pandas, Matplotlib, Seaborn, OpenPyXL

## Start

-   Startseite: http://127.0.0.1:5000/
-   Avocado: http://127.0.0.1:5000/avocado
-   Tomato: http://127.0.0.1:5000/tomato
-   Carrot: http://127.0.0.1:5000/ot
