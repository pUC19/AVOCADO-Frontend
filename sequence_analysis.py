# -*- coding: utf-8 -*-

# HINWEIS: Bevor Sie dieses Skript ausführen, müssen Sie die ViennaRNA-Bibliothek installieren.
# Öffnen Sie Ihr Terminal oder Ihre Kommandozeile und geben Sie folgenden Befehl ein:
# pip install vienna-rna

import RNA
import re

def calculate_contour_length(sequence, seq_type, structure=None):
    """
    Berechnet die Konturlänge für eine gegebene DNA- oder RNA-Sequenz.
    Für RNA wird die Sekundärstruktur berücksichtigt, um eine genauere Schätzung zu erhalten.

    Argumente:
    sequence (str): Die Nukleinsäuresequenz.
    seq_type (str): Der Typ der Sequenz ('DNA' oder 'RNA').
    structure (str, optional): Die Sekundärstruktur in Dot-Bracket-Notation. Nötig für RNA.

    Rückgabe:
    float: Die berechnete Konturlänge in Nanometern (nm).
    """
    length_nt = len(sequence)
    
    if seq_type == 'DNA':
        # Annahme: Die Länge eines Basenpaares in der B-DNA-Doppelhelix beträgt ca. 0.34 nm.
        # Die Konturlänge entspricht der Länge der Doppelhelix.
        contour_length = length_nt * 0.34
    elif seq_type == 'RNA':
        # Wenn eine Struktur vorhanden ist, wird die Länge genauer berechnet.
        if structure:
            # Zähle die Anzahl der ungepaarten Basen (Loops) und der Basenpaare (Stems)
            num_unpaired_bases = structure.count('.')
            num_base_pairs = structure.count('(')

            # Länge der einzelsträngigen Regionen (Loops): angenommener Abstand 0.59 nm/Base
            length_unpaired = num_unpaired_bases * 0.59
            
            # Länge der doppelsträngigen Regionen (Stems): angenommener Abstand 0.34 nm/Basenpaar
            length_paired = num_base_pairs * 0.34

            # Die Gesamtkonturlänge ist die Summe der Längen der Stems und der Loops
            contour_length = length_paired + length_unpaired
        else:
            # Fallback, falls keine Struktur übergeben wird (einfache Berechnung)
            contour_length = length_nt * 0.59
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
    
    rna_sequence = sequence.replace('T', 'U')

    # 3. Sekundärstruktur *zuerst* vorhersagen, da sie für die Längenberechnung benötigt wird
    (structure, mfe) = RNA.fold(rna_sequence)

    # 4. Konturlänge berechnen (mit der vorhergesagten Struktur für RNA)
    contour_length = calculate_contour_length(sequence, seq_type, structure)

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
        print("  (Berechnung berücksichtigt die Faltung: 0.59 nm/ungepaarte Base + 0.34 nm/Basenpaar)")
    
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

