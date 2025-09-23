# -*- coding: utf-8 -*-
# Final, complete version with "applyForce" visualization on click.
import dash
import dash_bio as dashbio
from dash import dcc, html
from dash.dependencies import Input, Output, State
import RNA
import re
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import json

# --- Polymer Physics Models & Helpers ---
KBT_PN_NM = 4.114

def wlc_extension(force, Lp, Lc):
    if force <= 0.01 or Lc <= 0: return 0
    ratio = force * Lp / KBT_PN_NM
    extension = Lc * (1 - 0.5 * np.sqrt(KBT_PN_NM / (force * Lp)) + force / 1000)
    return np.maximum(0, extension)

def fjc_extension(force, Lk, Lc):
    if force <= 0.01 or Lc <= 0: return 0
    x = force * Lk / KBT_PN_NM
    if x > 300: return Lc * (1 - 1/x)
    return Lc * (np.tanh(x)**-1 - 1/x)

def parse_stems_detailed(structure):
    bp_map = [-1] * len(structure)
    stack = []
    for i, char in enumerate(structure):
        if char == '(': stack.append(i)
        elif char == ')' and stack:
            j = stack.pop()
            bp_map[j], bp_map[i] = i, j
    
    stems, visited = [], [False] * len(structure)
    for i in range(len(structure)):
        if bp_map[i] > i and not visited[i]:
            current_stem_pairs = []
            curr_i, curr_j = i, bp_map[i]
            while curr_i < curr_j and bp_map[curr_i] == curr_j and not visited[curr_i]:
                visited[curr_i] = visited[curr_j] = True
                current_stem_pairs.append((curr_i, curr_j))
                curr_i, curr_j = curr_i + 1, curr_j - 1
            
            if len(current_stem_pairs) >= 3: # Minimal stem length
                stems.append({'pairs': current_stem_pairs, 'len': len(current_stem_pairs)})
    return stems

def analyze_and_simulate(sequence, handle_bp, params):
    sequence = sequence.upper().strip()
    if not re.match("^[ACGTU]+$", sequence): return {"error": "Ungültige Zeichen in der Sequenz."}
    if len(sequence) < 10: return {"error": "Sequenz zu kurz."}
    
    try:
        rna_sequence = sequence.replace('T', 'U')
        (structure, mfe) = RNA.fold(rna_sequence)
        
        RNA.pf_fold(rna_sequence)
        n = len(rna_sequence)
        prob_matrix = np.zeros((n, n))
        for i in range(1, n + 1):
            for j in range(i + 1, n + 1):
                prob = RNA.get_pr(i, j)
                if prob > 1e-4: prob_matrix[i-1, j-1] = prob_matrix[j-1, i-1] = prob

        analysis = {
            "sequence": sequence, "length_bp": n,
            "secondary_structure": structure, "mfe_kcal_mol": round(mfe, 2),
            "base_pairing_probabilities": prob_matrix.tolist(),
            "error": None, "stems_info": []
        }

        stems = parse_stems_detailed(structure)
        for s in stems:
            # Construct substring and substructure for energy calculation
            min_idx = min(p[0] for p in s['pairs'])
            max_idx = max(p[1] for p in s['pairs'])
            sub_seq = rna_sequence[min_idx : max_idx + 1]
            sub_struct_list = list('.' * len(sub_seq))
            for p1, p2 in s['pairs']:
                sub_struct_list[p1 - min_idx] = '('
                sub_struct_list[p2 - min_idx] = ')'
            sub_struct = "".join(sub_struct_list)
            
            s['energy_kcal'] = RNA.energy_of_struct(sub_seq, sub_struct)
            
            lc_unfolded_nt = max_idx - min_idx + 1
            s['lc_unfolded'] = lc_unfolded_nt * params['len_nt_ssrna']
            s['lc_folded'] = s['len'] * params['rise_bp_dsrna']
            s['delta_lc'] = s['lc_unfolded'] - s['lc_folded']
            
            energy_pn_nm = abs(s['energy_kcal']) * params['salt_correction_factor']
            s['f_unfold'] = (energy_pn_nm / s['delta_lc']) if s['delta_lc'] > 0.1 else 1000
        
        stems.sort(key=lambda x: x['f_unfold'])
        analysis['stems_info'] = [{'pairs': s['pairs'], 'f_unfold': s['f_unfold'], 'lc_unfolded': s['lc_unfolded'], 'lc_folded': s['lc_folded']} for s in stems]

        force_vec = np.linspace(0.1, 60, 300)
        lc_handles = handle_bp * params['rise_bp_dsdna']
        
        unfolding_x = []
        for f in force_vec:
            unfolded_stems_indices = {i for i, s in enumerate(stems) if f >= s['f_unfold']}
            lc_ssrna_total = sum(stems[i]['lc_unfolded'] for i in unfolded_stems_indices)
            lc_dsrna_total = sum(stems[i]['lc_folded'] for i in range(len(stems)) if i not in unfolded_stems_indices)
            
            ext_handles = wlc_extension(f, params['lp_dsdna'], lc_handles)
            ext_dsrna = wlc_extension(f, params['lp_dsrna'], lc_dsrna_total)
            ext_ssrna = fjc_extension(f, params['lk_ssrna'], lc_ssrna_total)
            unfolding_x.append(ext_handles + ext_dsrna + ext_ssrna)
        
        analysis['fd_unfolding'] = {'x': unfolding_x, 'y': force_vec.tolist()}
        
    except Exception as e:
        return {"error": f"Analysefehler: {e}"}
        
    return analysis

def add_dash_to_flask(server):
    app = dash.Dash(__name__, server=server, url_base_pathname='/dna_analyzer/',
                    meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}],
                    external_stylesheets=["https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css"])
    app.title = "RNA/DNA Analyzer"
    
    app.index_string = '''
    <!DOCTYPE html>
    <html lang="de" data-bs-theme="light">
        <head>
            {%metas%}
            <title>{%title%}</title>
            {%favicon%}
            {%css%}
            <style>
                :root{
                    --radius: 16px;
                    --shadow: 0 10px 30px rgba(0,0,0,.06), 0 2px 6px rgba(0,0,0,.05);
                    --card-bg: color-mix(in oklab, var(--bs-body-bg), transparent 10%);
                }
                [data-bs-theme="dark"]{
                    --shadow: 0 10px 30px rgba(0,0,0,.35), 0 2px 6px rgba(0,0,0,.25);
                    --card-bg: color-mix(in oklab, var(--bs-body-bg), white 4%);
                }
                body {
                    background:
                        radial-gradient(1000px 600px at 10% -10%, rgba(99,102,241,.18), transparent 40%),
                        radial-gradient(800px 500px at 90% 10%, rgba(14,165,233,.18), transparent 40%),
                        linear-gradient(to bottom, var(--bs-body-bg), var(--bs-body-bg));
                    background-attachment: fixed;
                }
                .card-glass{
                    background: var(--card-bg);
                    border: 1px solid color-mix(in oklab, var(--bs-body-color), transparent 90%);
                    border-radius: var(--radius);
                    box-shadow: var(--shadow);
                }
                .section-title{ letter-spacing:.3px }
                .dot-bracket-display { 
                    background-color: color-mix(in oklab, var(--bs-body-bg), transparent 20%);
                    padding: 1rem; border-radius: 8px; word-wrap: break-word;
                    border: 1px solid color-mix(in oklab, var(--bs-body-color), transparent 90%);
                }
                code { font-family: 'Courier New', Courier, monospace; font-size: 1.1em; }
                .nav-link-custom { text-decoration: none; display: inline-block; margin-bottom: 20px; font-size: 1.1em; }
                .nav-link-custom:hover { text-decoration: underline; }
                .table { --bs-table-bg: transparent; }
            </style>
        </head>
        <body>
            {%app_entry%}
            <footer>
                {%config%}
                {%scripts%}
                {%renderer%}
            </footer>
        </body>
    </html>
    '''

    app.layout = html.Div(className='container-xl my-4', children=[
        dcc.Store(id='analysis-results-store'),
        html.A("← Zurück zur Startseite", href="/", className="nav-link-custom"),
        html.H1("Interaktiver DNA/RNA Structure Analyzer", className="section-title h3 mb-4"),
        
        html.Div(className='row g-4', children=[
            html.Div(className='col-lg-4', children=[
                html.Div(className='card card-glass sticky-lg-top', style={'top': '20px'}, children=[
                    html.Div(className='card-body', children=[
                        dcc.Textarea(id='sequence-input', placeholder='RNA/DNA Sequenz hier eingeben...',
                                     className='form-control mb-3', style={'height': 100},
                                     value='GGGCUAUUACGAGCGGUGGAGUUUUUCCUGCGCCGGUAGGCGCCAGGAAGCUUU'),
                        html.H3("In Silico Mutagenese (Optional)", className="h6 mt-4 mb-3"),
                        html.Div(className='row g-3 align-items-center', children=[
                            html.Div(className='col-md-6', children=[
                                dcc.Input(id='mutation-position-input', type='number', placeholder='Position (1-basiert)', className='form-control')
                            ]),
                            html.Div(className='col-md-6', children=[
                                dcc.Dropdown(id='mutation-base-input',
                                             options=[{'label': b, 'value': b} for b in ['A', 'C', 'G', 'T', 'U']],
                                             placeholder='Neue Base', clearable=True)
                            ])
                        ]),
                        html.H3("Simulationsparameter", className="h6 mt-4 mb-3"),
                        html.Div(className='row g-3', children=[
                            html.Div(className='col-md-6', children=[
                                html.Label("Salz-Korrektur", htmlFor="salt-factor-input", className="form-label small"),
                                dcc.Input(id="salt-factor-input", type="number", value=2.9, step=0.1, className="form-control")
                            ]),
                             html.Div(className='col-md-6', children=[
                                html.Label("dsDNA Handles (bp)", htmlFor="handles-input", className="form-label small"),
                                dcc.Input(id="handles-input", type="number", value=1000, step=100, className="form-control")
                            ]),
                        ]),
                        html.Button('Analysieren & Simulieren', id='submit-button', n_clicks=0, className='btn btn-primary w-100 mt-4'),
                    ])
                ])
            ]),
            html.Div(className='col-lg-8', children=[
                 dcc.Loading(id="loading-spinner", type="default", children=html.Div(id="output-container"))
            ])
        ])
    ])

    @app.callback(
        Output('analysis-results-store', 'data'),
        Input('submit-button', 'n_clicks'),
        [State('sequence-input', 'value'),
         State('salt-factor-input', 'value'),
         State('handles-input', 'value'),
         State('mutation-position-input', 'value'),
         State('mutation-base-input', 'value')],
        prevent_initial_call=True
    )
    def run_analysis(n_clicks, sequence, salt_factor, handle_bp, mutation_pos, mutation_base):
        if not sequence: return dash.no_update
        params = {
            'salt_correction_factor': float(salt_factor if salt_factor is not None else 2.9),
            'lp_dsdna': 50.0, 'rise_bp_dsdna': 0.34, 'lp_dsrna': 63.0, 
            'rise_bp_dsrna': 0.29, 'lk_ssrna': 1.5, 'len_nt_ssrna': 0.59
        }
        results_wt = analyze_and_simulate(sequence, int(handle_bp if handle_bp is not None else 1000), params)
        results_mut = None
        if mutation_pos is not None and mutation_base:
            try:
                pos_0_based = int(mutation_pos) - 1
                if 0 <= pos_0_based < len(sequence):
                    s_list = list(sequence)
                    s_list[pos_0_based] = mutation_base
                    mutated_sequence = "".join(s_list)
                    results_mut = analyze_and_simulate(mutated_sequence, int(handle_bp), params)
            except Exception: pass
        return json.dumps({'wt': results_wt, 'mut': results_mut})

    @app.callback(
        Output('output-container', 'children'),
        Input('analysis-results-store', 'data')
    )
    def update_display(stored_data):
        if not stored_data:
            return html.Div("Bitte geben Sie eine Sequenz ein und starten Sie die Analyse.", className="alert alert-info")
        
        data = json.loads(stored_data)
        results_wt, results_mut = data.get('wt'), data.get('mut')

        if not results_wt or results_wt.get('error'):
            return html.Div(results_wt.get('error', "Unbekannter Fehler."), className="alert alert-danger")

        fd_fig = go.Figure()
        fd_fig.add_trace(go.Scatter(x=results_wt['fd_unfolding']['x'], y=results_wt['fd_unfolding']['y'], mode='lines', name='WT Unfolding', line=dict(color='royalblue')))
        if results_mut and not results_mut.get('error'):
             fd_fig.add_trace(go.Scatter(x=results_mut['fd_unfolding']['x'], y=results_mut['fd_unfolding']['y'], mode='lines', name='Mutant Unfolding', line=dict(color='skyblue', dash='dash')))
        fd_fig.update_layout(title='Simulierte FD-Kurve', xaxis_title='Extension (nm)', yaxis_title='Force (pN)', template="plotly_white", plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)', legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01))
        
        forna_wt = [{'sequence': results_wt['sequence'].replace('T', 'U'), 'structure': results_wt['secondary_structure']}]
        
        dot_plot_target = results_mut if results_mut and not results_mut.get('error') else results_wt
        prob_matrix = np.array(dot_plot_target['base_pairing_probabilities'])
        dot_plot_fig = go.Figure(data=go.Heatmap(z=np.sqrt(prob_matrix), colorscale='Blues', showscale=False))
        dot_plot_fig.update_layout(title='Dot-Plot der Basenpaarungswahrscheinlichkeit', yaxis_autorange='reversed', template="plotly_white", plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)')

        structure_comparison_content = []
        if results_mut and not results_mut.get('error'):
            forna_mut = [{'sequence': results_mut['sequence'].replace('T', 'U'), 'structure': results_mut['secondary_structure']}]
            structure_comparison_content = html.Div(className='row mt-3', children=[
                html.Div(className='col-lg-6', children=[html.H4("Wildtyp Struktur", className="h6 text-center mb-2"), dashbio.FornaContainer(sequences=forna_wt, height=500)]),
                html.Div(className='col-lg-6', children=[html.H4("Mutanten Struktur", className="h6 text-center mb-2"), dashbio.FornaContainer(sequences=forna_mut, height=500)])
            ])
        else:
            structure_comparison_content = html.Div(className='mt-3', children=[dashbio.FornaContainer(sequences=forna_wt, height=500)])

        summary_content = [html.H4("Wildtyp Analyse", className="h6"), html.P(f"MFE: {results_wt['mfe_kcal_mol']} kcal/mol")]
        if results_mut and not results_mut.get('error'):
            delta_mfe = results_mut['mfe_kcal_mol'] - results_wt['mfe_kcal_mol']
            summary_content.extend([html.Hr(), html.H4("Mutanten Analyse", className="h6"), html.P(f"MFE: {results_mut['mfe_kcal_mol']} kcal/mol"), html.P(f"ΔMFE: {delta_mfe:.2f} kcal/mol", className="fw-bold")])
            
        return html.Div(className='card card-glass mb-4', children=[
            html.Div(className='card-body', children=[
                dcc.Tabs(id="output-tabs", children=[
                    dcc.Tab(label='Interaktive Analyse', children=html.Div(className='row g-3 pt-3', children=[
                        html.Div(className='col-12', children=[dcc.Graph(id='fd-curve-graph', figure=fd_fig)]),
                        html.Div(className='col-12', children=[
                            html.P("Tippen Sie auf die Kurve, um die entfaltete Struktur zu sehen.", className="text-center text-muted small"),
                            dashbio.FornaContainer(id='forna-container-interactive', sequences=forna_wt)
                        ])
                    ])),
                    dcc.Tab(label='2D-Struktur Vergleich', children=structure_comparison_content),
                    dcc.Tab(label='Weitere Analysen', children=html.Div(className='p-3', children=[
                        dcc.Tabs([
                            dcc.Tab(label='Zusammenfassung', children=html.Div(className='pt-3', children=summary_content)),
                            dcc.Tab(label='Dot-Plot', children=dcc.Graph(figure=dot_plot_fig)),
                            dcc.Tab(label='Dot-Bracket', children=html.Div(className='dot-bracket-display mt-3', children=[
                                html.P("Sequenz:", className='fw-bold mb-1'), html.Code(dot_plot_target['sequence']),
                                html.P("Struktur:", className='fw-bold mb-1 mt-3'), html.Code(dot_plot_target['secondary_structure'])
                            ]))
                        ])
                    ]))
                ])
            ])
        ])
        
    @app.callback(
        Output('forna-container-interactive', 'sequences'),
        Input('fd-curve-graph', 'clickData'),
        State('analysis-results-store', 'data'),
        prevent_initial_call=True
    )
    def update_forna_on_click(clickData, stored_data):
        if not clickData or not stored_data: return dash.no_update
        
        data = json.loads(stored_data)
        point = clickData['points'][0]
        force_clicked = point['y']
        curve_index = point['curveNumber']
        
        target_analysis = data.get('mut') if curve_index == 1 and data.get('mut') and not data['mut'].get('error') else data.get('wt')
        
        if not target_analysis: return dash.no_update
        
        sequence = target_analysis['sequence'].replace('T', 'U')
        original_structure = target_analysis['secondary_structure']
        stems = target_analysis.get('stems_info', [])
        
        # Determine which base pairs are broken at the clicked force
        structure_list = list(original_structure)
        unfolded_stems_pairs = [
            pair for s in stems if s['f_unfold'] <= force_clicked for pair in s['pairs']
        ]

        for i, j in unfolded_stems_pairs:
            structure_list[i] = '.'
            structure_list[j] = '.'

        new_structure = "".join(structure_list)
        
        new_sequences = [{'sequence': sequence, 'structure': new_structure, 'options': {
            'applyForce': True,        # This is the key for the "pulled" look
            'circularizeExternal': False, # This makes it linear
            'labelInterval': 10,
        }}]
        return new_sequences

    return server