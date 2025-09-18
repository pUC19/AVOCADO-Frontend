import lumicks.pylake as lk
import tkinter as tk
from tkinter import filedialog

def inspect_lumicks_file(filename):
    """
    Öffnet eine Lumicks H5-Datei und listet die wichtigsten Daten-Gruppen auf.
    (Endgültige korrigierte Version: ändert f.keys() zu list(f) und entfernt f.close())
    """
    print(f"\n--- 🔎 Untersuche Datei: {filename} ---")
    
    try:
        # Öffne die Datei direkt. Pylake kümmert sich selbst um das Schließen.
        f = lk.File(filename)
        
        print("\n## 1. Gespeicherte FD-Kurven (f.fdcurves)")
        # f.fdcurves ist ein Dict-ähnliches Objekt, das .keys() hat
        if f.fdcurves:
            for name in f.fdcurves.keys():
                print(f"  - {name}")
        else:
            print("   (Keine vordefinierten FD-Kurven in dieser Datei)")

        print("\n## 2. Haupt-Daten-Gruppen (Root-Level)")
        # KORREKTUR: Das File-Objekt selbst wird direkt iteriert (wie list(f))
        root_groups = list(f)
        if root_groups:
            for group_name in root_groups:
                print(f"  - {group_name}")
        else:
            print("   (Keine Hauptgruppen gefunden)")

        print("\n## 3. Inhalte der relevanten Gruppen (Hier suchen!)")
        print("   Hier solltest du die Namen suchen, die POTATO nicht finden konnte:")

        if "Distance" in f:
            print("\n   Inhalt von 'Distance':")
            try:
                # Eine Gruppe (wie f['Distance']) kann man wieder direkt iterieren
                for dataset_name in list(f["Distance"]):
                    print(f"   ✅ Gefunden: Distance/{dataset_name}")
            except Exception as e:
                print(f"     Konnte 'Distance' nicht lesen: {e}")

        if "Force HF" in f:
            print("\n   Inhalt von 'Force HF':")
            try:
                for dataset_name in list(f["Force HF"]):
                    print(f"   - Force HF/{dataset_name}")
            except Exception as e:
                print(f"     Konnte 'Force HF' nicht lesen: {e}")

        if "Force LF" in f:
            print("\n   Inhalt von 'Force LF':")
            try:
                for dataset_name in list(f["Force LF"]):
                    print(f"   - Force LF/{dataset_name}")
            except Exception as e:
                print(f"     Konnte 'Force LF' nicht lesen: {e}")

    except Exception as e:
        print(f"\nFEHLER: Konnte Datei nicht verarbeiten. Ist es eine Lumicks H5-Datei?")
        print(f"   Details: {e}")
        # KORREKTUR: Der 'finally'-Block mit f.close() wurde komplett entfernt.

def main():
    """
    Hauptfunktion: Erstellt ein Tkinter-Fenster (unsichtbar) 
    und öffnet den Datei-Auswahl-Dialog. (Diese Funktion ist unverändert)
    """
    root = tk.Tk()
    root.withdraw()

    print("Bitte wähle eine Lumicks .h5 Datei zur Analyse aus...")

    file_path = filedialog.askopenfilename(
        title="Wähle eine Lumicks .h5 Datei",
        filetypes=(("HDF5 Dateien", "*.h5"), ("Alle Dateien", "*.*"))
    )

    if not file_path:
        print("\nKeine Datei ausgewählt. Programm wird beendet.")
    else:
        inspect_lumicks_file(file_path)

    print("\n--- Analyse abgeschlossen. ---")


# --- Script-Start ---
if __name__ == "__main__":
    main()
