# -*- coding: utf-8 -*-
# Final, complete, and corrected version with adjustable parameters and full HTML.
import dash
import dash_bio as dashbio
from dash import dcc, html
from dash.dependencies import Input, Output, State
import RNA
import re
import pandas as pd
import numpy as np
import plotly.graph_objects as go

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

def calculate_contour_length(sequence, seq_type, structure):
    if seq_type == 'DNA': return len(sequence) * 0.34
    num_unpaired = structure.count('.')
    num_pairs = structure.count('(')
    return (num_unpaired * 0.59) + (num_pairs * 0.34)

def parse_structure_final(structure):
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
            stem_pairs = []
            curr_i, curr_j = i, bp_map[i]
            while curr_i < curr_j and bp_map[curr_i] == curr_j and not visited[curr_i]:
                visited[curr_i] = visited[curr_j] = True
                stem_pairs.append((curr_i, curr_j))
                curr_i, curr_j = curr_i + 1, curr_j - 1
            if len(stem_pairs) >= 3:
                start_bp, end_bp = stem_pairs[0][0], stem_pairs[0][1]
                stems.append({'start': start_bp, 'end': end_bp, 'len': len(stem_pairs)})
    return stems

def analyze_and_simulate(sequence, handle_bp, params):
    sequence = sequence.upper().strip()
    if not re.match("^[ACGTU]+$", sequence): return {"error": "Invalid characters."}
    if len(sequence) < 10: return {"error": "Sequence too short."}
    
    try:
        rna_sequence = sequence.replace('T', 'U')
        (structure, mfe) = RNA.fold(rna_sequence)
        contour_len = calculate_contour_length(sequence, 'RNA' if 'U' in sequence else 'DNA', structure)

        RNA.pf_fold(rna_sequence)
        n = len(rna_sequence)
        prob_matrix = np.zeros((n, n))
        for i in range(1, n + 1):
            for j in range(i + 1, n + 1):
                prob = RNA.get_pr(i, j)
                if prob > 1e-4: prob_matrix[i-1, j-1] = prob_matrix[j-1, i-1] = prob

        analysis = {
            "sequence": sequence, "length_bp": n, "contour_length_nm": round(contour_len, 2),
            "secondary_structure": structure, "mfe_kcal_mol": round(mfe, 2),
            "base_pairing_probabilities": prob_matrix.tolist(), "error": None
        }

        stems = parse_structure_final(structure)
        for s in stems:
            sub_seq = rna_sequence[s['start'] : s['end'] + 1]
            sub_struct = structure[s['start'] : s['end'] + 1]
            s['energy_kcal'] = RNA.energy_of_struct(sub_seq, sub_struct)
            
            lc_unfolded_nt = s['end'] - s['start'] + 1
            s['lc_unfolded'] = lc_unfolded_nt * params['len_nt_ssrna']
            s['lc_folded'] = s['len'] * params['rise_bp_dsrna']
            s['delta_lc'] = s['lc_unfolded'] - s['lc_folded']
            
            energy_pn_nm = abs(s['energy_kcal']) * params['salt_correction_factor']
            s['f_unfold'] = (energy_pn_nm / s['delta_lc']) if s['delta_lc'] > 0.1 else 1000

        stems.sort(key=lambda x: x['f_unfold'])
        
        force_vec = np.linspace(0.1, 60, 300)
        lc_handles = handle_bp * params['rise_bp_dsdna']
        
        unfolding_x = []
        unfolded_stems_indices = set()
        for f in force_vec:
            newly_unfolded = {i for i, s in enumerate(stems) if i not in unfolded_stems_indices and f >= s['f_unfold']}
            unfolded_stems_indices.update(newly_unfolded)
            lc_ssrna_total = sum(stems[i]['lc_unfolded'] for i in unfolded_stems_indices)
            lc_dsrna_total = sum(stems[i]['lc_folded'] for i in range(len(stems)) if i not in unfolded_stems_indices)
            ext_handles = wlc_extension(f, params['lp_dsdna'], lc_handles)
            ext_dsrna = wlc_extension(f, params['lp_dsrna'], lc_dsrna_total)
            ext_ssrna = fjc_extension(f, params['lk_ssrna'], lc_ssrna_total)
            unfolding_x.append(ext_handles + ext_dsrna + ext_ssrna)
        
        refolding_x = []
        unfolded_stems_indices = set(range(len(stems)))
        for f in reversed(force_vec):
            f_refold_threshold = f * 1.2
            refolded = {i for i in unfolded_stems_indices if f_refold_threshold >= stems[i]['f_unfold']}
            unfolded_stems_indices -= refolded
            lc_ssrna_total = sum(stems[i]['lc_unfolded'] for i in unfolded_stems_indices)
            lc_dsrna_total = sum(stems[i]['lc_folded'] for i in range(len(stems)) if i not in unfolded_stems_indices)
            ext_handles = wlc_extension(f, params['lp_dsdna'], lc_handles)
            ext_dsrna = wlc_extension(f, params['lp_dsrna'], lc_dsrna_total)
            ext_ssrna = fjc_extension(f, params['lk_ssrna'], lc_ssrna_total)
            refolding_x.insert(0, ext_handles + ext_dsrna + ext_ssrna)

        analysis['fd_unfolding'] = {'x': unfolding_x, 'y': force_vec.tolist()}
        analysis['fd_refolding'] = {'x': refolding_x, 'y': force_vec.tolist()}
        
    except Exception as e:
        return {"error": f"Analysis Error: {e}"}
        
    return analysis

def add_dash_to_flask(server):
    app = dash.Dash(__name__, server=server, url_base_pathname='/dna_analyzer/',
                    meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}],
                    external_stylesheets=["https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css"])
    app.title = "RNA/DNA Analyzer"
    
    # This is the full, correct HTML index string required by Dash.
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
            <script>
                (function(){
                    const key = "potato-theme";
                    const saved = localStorage.getItem(key);
                    const prefersDark = window.matchMedia('(prefers-color-scheme: dark)').matches;
                    const setTheme = t => document.documentElement.setAttribute('data-bs-theme', t);
                    setTheme(saved || (prefersDark ? 'dark' : 'light'));
                })();
            </script>
        </body>
    </html>
    '''

    app.layout = html.Div(className='container-xl my-4', children=[
        html.A("← Zurück zur Startseite", href="/", className="nav-link-custom"),
        html.H1("DNA/RNA Structure Analyzer & FD Simulator", className="section-title h3 mb-4"),
        
        html.Div(className='card card-glass mb-4', children=[
            html.Div(className='card-body', children=[
                dcc.Textarea(id='sequence-input', placeholder='RNA/DNA sequence here...',
                             className='form-control mb-3', style={'height': 100},
                             value='TTTTACCTGATCAGGTGCTTTTTTTTGCACCTGATCAGGTATTT'),
                html.Button('Analysieren & Simulieren', id='submit-button', n_clicks=0, className='btn btn-primary'),
            ])
        ]),

        html.Div(className='card card-glass mb-4', children=[
            html.Div(className='card-body', children=[
                html.H3("Simulationsparameter", className="h5 mb-3"),
                html.Div(className='row g-3', children=[
                    html.Div(className='col-md-4', children=[
                        html.Label("Salz-Korrekturfaktor (für ΔG)", htmlFor="salt-factor-input", className="form-label"),
                        dcc.Input(id="salt-factor-input", type="number", value=2.9, step=0.1, className="form-control")
                    ]),
                    html.Div(className='col-md-4', children=[
                        html.Label("dsDNA Handles (bp)", htmlFor="handles-input", className="form-label"),
                        dcc.Input(id="handles-input", type="number", value=1000, step=100, className="form-control")
                    ]),
                    html.Div(className='col-md-4', children=[
                        html.Label("dsDNA Lp (nm)", htmlFor="lp-dsdna-input", className="form-label"),
                        dcc.Input(id="lp-dsdna-input", type="number", value=50.0, step=1, className="form-control")
                    ]),
                ])
            ])
        ]),
        
        dcc.Loading(id="loading-spinner", type="circle", children=html.Div(id="output-container"))
    ])

    @app.callback(
        Output('output-container', 'children'),
        Input('submit-button', 'n_clicks'),
        [State('sequence-input', 'value'),
         State('salt-factor-input', 'value'),
         State('handles-input', 'value'),
         State('lp-dsdna-input', 'value')]
    )
    def update_all_outputs(n_clicks, sequence, salt_factor, handle_bp, lp_dsdna):
        if not n_clicks: return ""
        if not sequence: return html.Div("Bitte geben Sie eine Sequenz ein.", className='alert alert-warning mt-3')

        params = {
            'salt_correction_factor': float(salt_factor if salt_factor is not None else 2.9),
            'lp_dsdna': float(lp_dsdna if lp_dsdna is not None else 50.0),
            'rise_bp_dsdna': 0.34,
            'lp_dsrna': 63.0, 'rise_bp_dsrna': 0.29,
            'lk_ssrna': 1.5, 'len_nt_ssrna': 0.59
        }
        
        results = analyze_and_simulate(sequence, handle_bp=int(handle_bp if handle_bp is not None else 1000), params=params)
        
        if results.get('error'): return html.Div(results['error'], className='alert alert-danger mt-3')

        results_df = pd.DataFrame({
            "Parameter": ["Länge (Basen)", "Konturlänge (nm)", "MFE (kcal/mol)"],
            "Value": [results['length_bp'], f"~ {results['contour_length_nm']}", results['mfe_kcal_mol']]
        })
        forna_sequences = [{'sequence': results['sequence'].replace('T', 'U'), 'structure': results['secondary_structure'],
                            'options': {'applyForce': True, 'circularizeExternal': True, 'labelInterval': 5}}]
        prob_matrix = np.array(results['base_pairing_probabilities'])
        dot_plot_fig = go.Figure(data=go.Heatmap(z=np.sqrt(prob_matrix), colorscale='Blues', showscale=False))
        dot_plot_fig.update_layout(title='Dot-Plot (√Probability)', yaxis_autorange='reversed', template="plotly_white", plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)')
        
        fd_fig = go.Figure()
        fd_fig.add_trace(go.Scatter(x=results['fd_unfolding']['x'], y=results['fd_unfolding']['y'], mode='lines', name='Unfolding', line=dict(color='royalblue')))
        fd_fig.add_trace(go.Scatter(x=results['fd_refolding']['x'], y=results['fd_refolding']['y'], mode='lines', name='Refolding', line=dict(color='crimson')))
        fd_fig.update_layout(title='Simulierte FD-Kurve', xaxis_title='Extension (nm)', yaxis_title='Force (pN)', template="plotly_white", plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)')

        return html.Div(className='row g-4 mt-2', children=[
            html.Div(className='col-12 col-lg-7', children=[
                html.Div(className='card card-glass h-100', children=[
                    html.Div(className='card-body', children=[
                        dcc.Tabs(children=[
                            dcc.Tab(label='FD Simulation', children=dcc.Graph(figure=fd_fig)),
                            dcc.Tab(label='2D-Struktur (MFE)', children=dashbio.FornaContainer(sequences=forna_sequences, height=500)),
                            dcc.Tab(label='Dot-Plot', children=dcc.Graph(figure=dot_plot_fig)),
                        ])
                    ])
                ])
            ]),
            html.Div(className='col-12 col-lg-5', children=[
                html.Div(className='card card-glass h-100', children=[
                    html.Div(className='card-body', children=[
                        html.H3("Analyse & Simulation", className="h5"),
                        dcc.Tabs(className="mt-3", children=[
                            dcc.Tab(label='Zusammenfassung', children=html.Div(className="mt-3", children=[
                                html.Table(className="table", children=[
                                    html.Thead(html.Tr([html.Th(col) for col in results_df.columns])),
                                    html.Tbody([html.Tr([html.Td(results_df.iloc[i][col]) for col in results_df.columns]) for i in range(len(results_df))])
                                ])
                            ])),
                            dcc.Tab(label='Dot-Bracket', children=html.Div(className='dot-bracket-display mt-3', children=[
                                html.P("Sequenz:", className='fw-bold mb-1'), html.Code(results['sequence']),
                                html.P("Struktur:", className='fw-bold mb-1 mt-3'), html.Code(results['secondary_structure'])
                            ])),
                            dcc.Tab(label='Parameter', children=html.Div(className="mt-3", children=[
                                html.P("Die Simulation wurde mit folgenden Parametern durchgeführt:"),
                                html.Ul([
                                    html.Li(f"Salz-Korrekturfaktor: {params['salt_correction_factor']}"),
                                    html.Li(f"dsDNA Handles: {handle_bp} bp, Lp = {params['lp_dsdna']} nm"),
                                    html.Li(f"dsRNA: Lp = {params['lp_dsrna']} nm"),
                                    html.Li(f"ssRNA: Kuhn-Länge = {params['lk_ssrna']} nm"),
                                ])
                            ]))
                        ])
                    ])
                ])
            ])
        ])

    return server