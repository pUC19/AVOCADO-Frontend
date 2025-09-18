# -*- coding: utf-8 -*-

# HINWEIS: Bevor Sie dieses Skript ausführen, müssen Sie die ViennaRNA-Bibliothek installieren.
# Öffnen Sie Ihr Terminal oder Ihre Kommandozeile und geben Sie folgenden Befehl ein:
# pip install vienna-rna

import RNA
import re

def calculate_contour_length(sequence, seq_type):
    """
    Berechnet die ungefähre Konturlänge für eine gegebene DNA- oder RNA-Sequenz.

    Argumente:
    sequence (str): Die Nukleinsäuresequenz.
    seq_type (str): Der Typ der Sequenz ('DNA' oder 'RNA').

    Rückgabe:
    float: Die berechnete Konturlänge in Nanometern (nm).
    """
    length_bp = len(sequence)
    
    if seq_type == 'DNA':
        # Annahme: Die Länge eines Basenpaares in der B-DNA-Doppelhelix beträgt ca. 0.34 nm.
        # Die Konturlänge entspricht der Länge der Doppelhelix.
        contour_length = length_bp * 0.34
    elif seq_type == 'RNA':
        # Annahme: Der Abstand zwischen den Basen in einer einzelsträngigen RNA
        # beträgt ca. 0.5 bis 0.7 nm. Wir verwenden hier einen Mittelwert von 0.6 nm.
        contour_length = length_bp * 0.6
    else:
        contour_length = 0.0
        
    return contour_length

def analyze_sequence(sequence):
    """
    Analysiert eine DNA- oder RNA-Sequenz, um Typ, Konturlänge und Sekundärstruktur zu bestimmen.

    Argumente:
    sequence (str): Die zu analysierende Nukleinsäuresequenz.

    Rückgabe:
    dict: Ein Wörterbuch mit den Analyseergebnissen.
    """
    # 1. Eingabe validieren und normalisieren (Großbuchstaben, Leerzeichen entfernen)
    sequence = sequence.upper().strip()
    if not re.match("^[ACGTU]+$", sequence):
        return {"error": "Die Sequenz enthält ungültige Zeichen. Nur A, C, G, T, U sind erlaubt."}

    # 2. Sequenztyp bestimmen (DNA oder RNA)
    seq_type = 'RNA' if 'U' in sequence else 'DNA'

    # Für die Faltung muss DNA in RNA umgewandelt werden (T -> U)
    rna_sequence = sequence.replace('T', 'U')

    # 3. Konturlänge berechnen
    contour_length = calculate_contour_length(sequence, seq_type)

    # 4. Sekundärstruktur vorhersagen mit ViennaRNA
    # RNA.fold() gibt die Struktur in "Dot-Bracket"-Notation und die minimale freie Energie (MFE) zurück.
    (structure, mfe) = RNA.fold(rna_sequence)

    return {
        "sequence": sequence,
        "type": seq_type,
        "length_bp": len(sequence),
        "contour_length_nm": round(contour_length, 2),
        "secondary_structure": structure,
        "mfe_kcal_mol": round(mfe, 2)
    }

def print_results(results):
    """
    Gibt die Analyseergebnisse formatiert in der Konsole aus.
    """
    if "error" in results:
        print(f"\nFehler: {results['error']}")
        return

    print("\n--- Analyseergebnis ---")
    print(f"Sequenztyp:           {results['type']}")
    print(f"Länge:                {results['length_bp']} Basen/Basenpaare")
    print(f"Konturlänge:          ca. {results['contour_length_nm']} nm")
    if results['type'] == 'DNA':
        print("  (Berechnung basiert auf 0.34 nm pro Basenpaar in einer Doppelhelix)")
    else:
        print("  (Berechnung basiert auf 0.6 nm pro Base in einem Einzelstrang)")
    
    print(f"\nMinimale freie Energie (MFE): {results['mfe_kcal_mol']} kcal/mol")
    print("  (Je negativer der Wert, desto stabiler ist die vorhergesagte Struktur)")

    print("\nVorhergesagte Sekundärstruktur (Dot-Bracket-Notation):")
    print(f"Sequenz:    {results['sequence']}")
    print(f"Struktur:   {results['secondary_structure']}")
    
    print("\n--- Erklärung der Notation ---")
    print("  ( ) : Eine öffnende und eine schließende Klammer markieren ein Basenpaar, das einen 'Stamm' (Stem) bildet.")
    print("  .   : Ein Punkt steht für eine ungepaarte Base in einer 'Schleife' (Loop), z.B. Hairpin-, Bulge- oder Internal-Loop.")
    print("---------------------------\n")


# Hauptteil des Skripts
if __name__ == "__main__":
    print("=====================================================")
    print("   DNA/RNA Konturlängen- & Strukturanalyse-Tool")
    print("=====================================================")
    print("Geben Sie eine DNA- oder RNA-Sequenz ein und drücken Sie Enter.")
    print("Beispiel: GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA\n")
    
    try:
        user_sequence = input("Ihre Sequenz: ")
        if user_sequence:
            analysis_results = analyze_sequence(user_sequence)
            print_results(analysis_results)
        else:
            print("\nFehler: Keine Sequenz eingegeben. Bitte starten Sie das Skript erneut.")
    except Exception as e:
        print(f"\nEin unerwarteter Fehler ist aufgetreten: {e}")
        print("Stellen Sie sicher, dass die 'vienna-rna' Bibliothek korrekt installiert ist.")
