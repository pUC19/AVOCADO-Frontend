# -*- coding: utf-8 -*-
import dash
import dash_bio as dashbio
from dash import dcc, html
from dash.dependencies import Input, Output, State
import RNA
import re
import pandas as pd

# --- Backend Analysis Logic ---
def calculate_contour_length(sequence, seq_type, structure=None):
    """
    Calculates the contour length for a given DNA or RNA sequence.
    For RNA, the secondary structure is considered for a more accurate estimate.
    """
    length_nt = len(sequence)
    contour_length = 0.0
    
    if seq_type == 'DNA':
        contour_length = length_nt * 0.34
    elif seq_type == 'RNA':
        if structure:
            num_unpaired_bases = structure.count('.')
            num_base_pairs = structure.count('(')
            length_unpaired = num_unpaired_bases * 0.59
            length_paired = num_base_pairs * 0.34
            contour_length = length_paired + length_unpaired
        else:
            contour_length = length_nt * 0.59
            
    return contour_length

def analyze_sequence(sequence):
    """
    Analyzes a DNA or RNA sequence to determine its type, contour length,
    and secondary structure.
    """
    sequence = sequence.upper().strip()
    if not re.match("^[ACGTU]+$", sequence):
        return {"error": "Invalid characters in sequence. Only A, C, G, T, U are allowed."}
    if len(sequence) < 4:
        return {"error": "Sequence is too short. Please enter at least 4 nucleotides."}

    seq_type = 'RNA' if 'U' in sequence else 'DNA'
    rna_sequence = sequence.replace('T', 'U')

    (structure, mfe) = RNA.fold(rna_sequence)

    contour_length = calculate_contour_length(sequence, seq_type, structure)

    return {
        "sequence": sequence,
        "type": seq_type,
        "length_bp": len(sequence),
        "contour_length_nm": round(contour_length, 2),
        "secondary_structure": structure,
        "mfe_kcal_mol": round(mfe, 2),
        "error": None
    }


def add_dash_to_flask(server):
    """
    Initializes the Dash app and mounts it on the provided Flask server.
    """
    app = dash.Dash(
        __name__,
        server=server,
        url_base_pathname='/dna_analyzer/',
        meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}],
        external_stylesheets=[
            {
                "href": "https://fonts.googleapis.com/css2?family=Roboto:wght@300;400;500&display=swap",
                "rel": "stylesheet",
            },
        ]
    )
    app.title = "Sequence Structure Analyzer"

    # Custom CSS and a "Back" button to integrate with the main Flask app's look and feel
    app.index_string = '''
    <!DOCTYPE html>
    <html>
        <head>
            {%metas%}
            <title>{%title%}</title>
            {%favicon%}
            {%css%}
            <style>
                body { font-family: 'Roboto', sans-serif; background-color: #f8f9fa; color: #333; margin: 0; padding: 20px; }
                .container { max-width: 1200px; margin: auto; }
                .header { text-align: center; margin-bottom: 30px; border-bottom: 2px solid #dee2e6; padding-bottom: 20px; }
                h1 { color: #0056b3; }
                .input-section { display: flex; gap: 10px; align-items: center; margin-bottom: 30px; }
                .button { background-color: #007bff; color: white; border: none; padding: 10px 20px; border-radius: 5px; cursor: pointer; font-size: 16px; }
                .button:hover { background-color: #0056b3; }
                .error-message { color: #d9534f; background-color: #f2dede; border: 1px solid #ebccd1; padding: 15px; border-radius: 5px; text-align: center; }
                .results-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(400px, 1fr)); gap: 20px; }
                .result-card { background-color: white; border-radius: 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); padding: 20px; }
                .forna-container svg { width: 100%; height: auto; }
                table { width: 100%; border-collapse: collapse; margin-top: 15px; }
                th, td { text-align: left; padding: 12px; border-bottom: 1px solid #ddd; }
                th { background-color: #f2f2f2; }
                .info-text, .explanation-text { font-size: 0.9em; color: #666; margin-top: 15px; }
                .dot-bracket-display { margin-top: 15px; background-color: #e9ecef; padding: 10px; border-radius: 5px; word-wrap: break-word; }
                .code-label { font-weight: bold; margin-bottom: 5px; }
                code { font-family: 'Courier New', Courier, monospace; font-size: 1.1em; }
            </style>
        </head>
        <body>
            <div class="container">
                <a href="/" style="text-decoration: none; color: #007bff; display: inline-block; margin-bottom: 20px; font-size: 1.1em;">&larr; Zurück zur Hauptseite</a>
            </div>
            {%app_entry%}
            <footer>
                {%config%}
                {%scripts%}
                {%renderer%}
            </footer>
        </body>
    </html>
    '''

    # --- Dash App Layout ---
    app.layout = html.Div(className='container', children=[
        html.Div(className='header', children=[
            html.H1("DNA/RNA Sequence & Structure Analyzer"),
            html.P("Enter a nucleic acid sequence to predict its secondary structure, MFE, and contour length.")
        ]),
        html.Div(className='input-section', children=[
            dcc.Textarea(
                id='sequence-input',
                placeholder='Enter sequence here...\nExample: GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA',
                style={'width': '100%', 'height': 100},
            ),
            html.Button('Analyze Sequence', id='submit-button', n_clicks=0, className='button'),
        ]),
        dcc.Loading(
            id="loading-spinner",
            type="circle",
            children=html.Div(id="output-container")
        )
    ])

    # --- Dash Callback for Interactivity ---
    @app.callback(
        Output('output-container', 'children'),
        Input('submit-button', 'n_clicks'),
        State('sequence-input', 'value'),
        prevent_initial_call=True
    )
    def update_output(n_clicks, sequence):
        if not sequence:
            return html.Div("Please enter a sequence to analyze.", className='error-message')

        results = analyze_sequence(sequence)

        if results['error']:
            return html.Div(results['error'], className='error-message')

        results_df = pd.DataFrame({
            "Parameter": ["Sequence Type", "Length (bases)", "Contour Length (nm)", "MFE (kcal/mol)"],
            "Value": [results['type'], results['length_bp'], f"~ {results['contour_length_nm']}", results['mfe_kcal_mol']]
        })

        forna_sequences = [{
            'sequence': results['sequence'].replace('T', 'U'),
            'structure': results['secondary_structure'],
            'options': {
                'applyForce': True,
                'circularizeExternal': True,
                'labelInterval': 5,
                'initialSize': [600, 500]
            }
        }]
        
        return html.Div(className='results-grid', children=[
            html.Div(className='result-card forna-container', children=[
                html.H3("Predicted Secondary Structure"),
                dashbio.FornaContainer(id='forna-viewer', sequences=forna_sequences)
            ]),
            html.Div(className='result-card', children=[
                html.H3("Analysis Results"),
                dcc.Tabs(children=[
                    dcc.Tab(label='Summary', children=[
                        html.Table([
                            html.Thead(html.Tr([html.Th(col) for col in results_df.columns])),
                            html.Tbody([
                                html.Tr([html.Td(results_df.iloc[i][col]) for col in results_df.columns])
                                for i in range(len(results_df))
                            ])
                        ]),
                        html.P(
                            "Calculation based on 0.59 nm/unpaired base and 0.34 nm/base-pair.",
                            className='info-text'
                        ) if results['type'] == 'RNA' else html.P(
                            "Calculation based on 0.34 nm per base-pair in a double helix.",
                            className='info-text'
                        )
                    ]),
                    dcc.Tab(label='Dot-Bracket', children=[
                        html.Div(className='dot-bracket-display', children=[
                            html.P("Sequence:", className='code-label'),
                            html.Code(results['sequence']),
                            html.P("Structure:", className='code-label'),
                            html.Code(results['secondary_structure'])
                        ])
                    ]),
                     dcc.Tab(label='Notation Guide', children=[
                        html.Div(className='explanation-text', children=[
                            html.P("💡 In the dot-bracket notation:"),
                            html.Ul([
                               html.Li([html.B("( )"), ": A pair of parentheses marks a base pair forming a 'stem'."]),
                               html.Li([html.B("."), ": A dot represents an unpaired base in a 'loop'."])
                            ]),
                            html.P("💡 The MFE (Minimum Free Energy) indicates the stability of the predicted structure. More negative values imply greater stability.")
                        ])
                    ])
                ])
            ])
        ])

    return server