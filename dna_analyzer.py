# -*- coding: utf-8 -*-
# Final version restoring FornaContainer and keeping the coupled linear annotation view.
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
            
            if len(current_stem_pairs) >= 3:
                stems.append({'pairs': current_stem_pairs, 'len': len(current_stem_pairs)})
    return stems

def run_mutation_scan(sequence, wt_mfe):
    bases = ['A', 'C', 'G', 'U']
    mutations = []
    rna_sequence = sequence.replace('T', 'U')
    seq_list = list(rna_sequence)
    for i in range(len(seq_list)):
        original_base = seq_list[i]
        for new_base in bases:
            if new_base == original_base:
                mutations.append({'Position': i + 1, 'Mutation': f"{original_base} > {new_base}", 'dMFE': 0.0})
                continue
            seq_list[i] = new_base
            mut_seq = "".join(seq_list)
            (_, mut_mfe) = RNA.fold(mut_seq)
            dMFE = round(mut_mfe - wt_mfe, 2)
            mutations.append({'Position': i + 1, 'Mutation': f"{original_base} > {new_base}", 'dMFE': dMFE})
        seq_list[i] = original_base
    return pd.DataFrame(mutations)

def parse_annotations(annotation_string):
    annotations = []
    if not annotation_string: return annotations
    for line in [line.strip() for line in annotation_string.split('\n') if line.strip()]:
        if '=' in line:
            name, parts = [p.strip() for p in line.split('=', 1)]
            for r_str in [r.strip() for r in parts.split(',')]:
                if '-' in r_str:
                    try:
                        start, end = map(int, r_str.split('-'))
                        if start > 0 and end >= start: annotations.append({'name': name, 'start': start, 'end': end})
                    except ValueError: continue
    return annotations

def analyze_and_simulate(sequence, handle_bp, params, annotation_string):
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
        mutation_df = run_mutation_scan(rna_sequence, mfe)
        annotations = parse_annotations(annotation_string)
        analysis = {"sequence": sequence, "length_bp": n, "secondary_structure": structure, "mfe_kcal_mol": round(mfe, 2), "base_pairing_probabilities": prob_matrix.tolist(), "mutation_scan_results": mutation_df.to_dict('records'), "annotations": annotations, "error": None, "stems_info": []}
        stems = parse_stems_detailed(structure)
        for s in stems:
            min_idx = min(p[0] for p in s['pairs']); max_idx = max(p[1] for p in s['pairs'])
            sub_seq = rna_sequence[min_idx : max_idx + 1]
            sub_struct_list = list('.' * len(sub_seq));
            for p1, p2 in s['pairs']: sub_struct_list[p1 - min_idx] = '('; sub_struct_list[p2 - min_idx] = ')'
            sub_struct = "".join(sub_struct_list)
            s['energy_kcal'] = RNA.energy_of_struct(sub_seq, sub_struct)
            lc_unfolded_nt = max_idx - min_idx + 1
            s['lc_unfolded'] = lc_unfolded_nt * params['len_nt_ssrna']; s['lc_folded'] = s['len'] * params['rise_bp_dsrna']
            s['delta_lc'] = s['lc_unfolded'] - s['lc_folded']
            energy_pn_nm = abs(s['energy_kcal']) * params['salt_correction_factor']
            s['f_unfold'] = (energy_pn_nm / s['delta_lc']) if s['delta_lc'] > 0.1 else 1000
        stems.sort(key=lambda x: x['f_unfold'])
        analysis['stems_info'] = [{'pairs': s['pairs'], 'f_unfold': s['f_unfold']} for s in stems]
        force_vec = np.linspace(0.1, 60, 300)
        lc_handles = handle_bp * params['rise_bp_dsdna']
        unfolding_x = []
        for f in force_vec:
            unfolded_stems_indices = {i for i, s in enumerate(stems) if f >= s['f_unfold']}
            lc_ssrna_total = sum(stems[i]['lc_unfolded'] for i in unfolded_stems_indices) if stems and 'lc_unfolded' in stems[0] else 0
            lc_dsrna_total = sum(stems[i]['lc_folded'] for i,s in enumerate(stems) if i not in unfolded_stems_indices) if stems and 'lc_folded' in stems[0] else 0
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
            {%metas%}<title>{%title%}</title>{%favicon%}{%css%}
            <style>
                :root{ --radius: 16px; --shadow: 0 10px 30px rgba(0,0,0,.06); --card-bg: color-mix(in oklab, var(--bs-body-bg), transparent 10%); }
                [data-bs-theme="dark"]{ --shadow: 0 10px 30px rgba(0,0,0,.35); --card-bg: color-mix(in oklab, var(--bs-body-bg), white 4%); }
                body { background: radial-gradient(1000px 600px at 10% -10%, rgba(99,102,241,.18), transparent 40%), radial-gradient(800px 500px at 90% 10%, rgba(14,165,233,.18), transparent 40%), linear-gradient(to bottom, var(--bs-body-bg), var(--bs-body-bg)); background-attachment: fixed; }
                .card-glass{ background: var(--card-bg); border: 1px solid color-mix(in oklab, var(--bs-body-color), transparent 90%); border-radius: var(--radius); box-shadow: var(--shadow); }
                .section-title{ letter-spacing:.3px }
                .dot-bracket-display { background-color: color-mix(in oklab, var(--bs-body-bg), transparent 20%); padding: 1rem; border-radius: 8px; word-wrap: break-word; border: 1px solid color-mix(in oklab, var(--bs-body-color), transparent 90%); }
                .linear-sequence-view { font-family: 'Courier New', Courier, monospace; font-size: 1.1em; padding: 1rem; border-radius: 8px; line-height: 2; }
                .nav-link-custom { text-decoration: none; display: inline-block; margin-bottom: 20px; font-size: 1.1em; }
                .nav-link-custom:hover { text-decoration: underline; }
            </style>
        </head>
        <body> {%app_entry%} <footer> {%config%} {%scripts%} {%renderer%} </footer> </body>
    </html>
    '''
    app.layout = html.Div(className='container-xl my-4', children=[
        dcc.Store(id='analysis-results-store'),
        html.A("← Zurück zur Startseite", href="/", className="nav-link-custom"),
        html.H1("RNA/DNA Analyzer", className="h3 mb-4 section-title"),
        html.Div(className='row g-4', children=[
            html.Div(className='col-lg-4', children=[
                html.Div(className='card card-glass sticky-lg-top', style={'top': '20px'}, children=[
                    html.Div(className='card-body', children=[
                        dcc.Textarea(id='sequence-input', placeholder='RNA/DNA Sequenz hier eingeben...', className='form-control mb-3', style={'height': 100}, value='GGGCUAUUACGAGCGGUGGAGUUUUUCCUGCGCCGGUAGGCGCCAGGAAGCUUU'),
                        html.H3("Struktur-Annotationen (Optional)", className="h6 mt-4 mb-3"),
                        dcc.Textarea(id='annotation-input', placeholder='Beispiel:\nSL1 = 10-25\nBindestelle = 30-40', className='form-control mb-3', style={'height': 80}),
                        html.H3("In Silico Mutagenese (Optional)", className="h6 mt-4 mb-3"),
                        dcc.Input(id='mutation-position-input', type='number', placeholder='Position (1-based)', className='form-control mb-2'),
                        dcc.Dropdown(id='mutation-base-input', options=[{'label': b, 'value': b} for b in ['A', 'C', 'G', 'T', 'U']], placeholder='New Base', clearable=True),
                        html.H3("Simulationsparameter", className="h6 mt-4 mb-3"),
                        dcc.Input(id="salt-factor-input", type="number", value=2.9, step=0.1, className="form-control mb-2", placeholder="Salz-Korrektur"),
                        dcc.Input(id="handles-input", type="number", value=1000, step=100, className="form-control", placeholder="dsDNA Handles (bp)"),
                        html.Button('Analysieren & Simulieren', id='submit-button', n_clicks=0, className='btn btn-primary w-100 mt-4'),
                    ])])]),
            html.Div(className='col-lg-8', children=[dcc.Loading(id="loading-spinner", type="default", children=html.Div(id="output-container"))])])])
    @app.callback(Output('analysis-results-store', 'data'), Input('submit-button', 'n_clicks'), [State('sequence-input', 'value'), State('salt-factor-input', 'value'), State('handles-input', 'value'), State('mutation-position-input', 'value'), State('mutation-base-input', 'value'), State('annotation-input', 'value')], prevent_initial_call=True)
    def run_master_analysis(n_clicks, sequence, salt_factor, handle_bp, mutation_pos, mutation_base, annotation_string):
        if not sequence: return dash.no_update
        params = {'salt_correction_factor': float(salt_factor or 2.9), 'lp_dsdna': 50.0, 'rise_bp_dsdna': 0.34, 'lp_dsrna': 63.0, 'rise_bp_dsrna': 0.29, 'lk_ssrna': 1.5, 'len_nt_ssrna': 0.59}
        results_wt = analyze_and_simulate(sequence, int(handle_bp or 1000), params, annotation_string)
        results_mut = None
        if mutation_pos is not None and mutation_base:
            try:
                pos_0_based = int(mutation_pos) - 1
                if 0 <= pos_0_based < len(sequence):
                    mut_seq = list(sequence); mut_seq[pos_0_based] = mutation_base; mutated_sequence = "".join(mut_seq)
                    results_mut = analyze_and_simulate(mutated_sequence, int(handle_bp), params, annotation_string)
            except Exception: pass
        return json.dumps({'wt': results_wt, 'mut': results_mut})
    @app.callback(Output('output-container', 'children'), Input('analysis-results-store', 'data'))
    def update_main_output(stored_data):
        if not stored_data: return html.Div("Bitte Sequenz eingeben und Analyse starten.", className="alert alert-info")
        data = json.loads(stored_data)
        results_wt = data.get('wt')
        if not results_wt or results_wt.get('error'): return html.Div(results_wt.get('error', "Unbekannter Fehler."), className="alert alert-danger")
        fd_fig, heatmap_fig, dot_plot_fig = build_fd_figure(data), build_heatmap_figure(results_wt), build_dot_plot_figure(data)
        return html.Div(className='card card-glass', children=[html.Div(className='card-body', children=[dcc.Tabs(id="output-tabs", children=[build_interactive_analysis_tab(fd_fig, results_wt), dcc.Tab(label='Mutations-Heatmap', children=build_description_card("Was sehe ich hier?", "Diese Heatmap zeigt die Auswirkung jeder möglichen einzelnen Mutation. Rote Bereiche sind stark destabilisierend (schlecht für die Faltung), blaue Bereiche stabilisierend. So finden Sie strukturelle 'Hotspots'.") + [dcc.Graph(figure=heatmap_fig)]), build_comparison_tab(data), build_other_analyses_tab(data, dot_plot_fig)])])])
    @app.callback(Output('forna-container-interactive', 'sequences'), Input('fd-curve-graph', 'clickData'), State('analysis-results-store', 'data'), prevent_initial_call=True)
    def update_forna_on_click(clickData, stored_data):
        if not clickData or not stored_data: return dash.no_update
        data = json.loads(stored_data); point = clickData['points'][0]; force_clicked, curve_index = point['y'], point['curveNumber']
        target = data.get('mut') if curve_index == 1 and data.get('mut') and not data['mut'].get('error') else data.get('wt')
        if not target: return dash.no_update
        seq, struct, stems = target['sequence'].replace('T', 'U'), target['secondary_structure'], target.get('stems_info', [])
        struct_list = list(struct); unfolded_pairs = [pair for s in stems if s['f_unfold'] <= force_clicked for pair in s['pairs']]
        for i, j in unfolded_pairs: struct_list[i], struct_list[j] = '.', '.'
        color_string, _, _ = get_annotation_colors(target.get('annotations', []))
        return [{'sequence': seq, 'structure': "".join(struct_list), 'options': {'applyForce': True, 'circularizeExternal': False, 'labelInterval': 10, 'customColors': {'domain': [0, 1], 'range': ['#ffffff', '#ffffff'], 'colorValues': {'MyColoring': color_string}, 'colorAccessor': 'MyColoring'}}}]
    def get_annotation_colors(annotations):
        colors = ['rgba(99, 110, 250, 0.4)', 'rgba(239, 85, 59, 0.4)', 'rgba(0, 204, 150, 0.4)', 'rgba(171, 99, 250, 0.4)', 'rgba(255, 161, 90, 0.4)', 'rgba(25, 211, 243, 0.4)']
        color_string_parts, legend_items = [], []
        unique_names = sorted(list(set(a['name'] for a in annotations)))
        name_to_color = {name: colors[i % len(colors)] for i, name in enumerate(unique_names)}
        for name, color in name_to_color.items():
            ranges = [f"{a['start']}-{a['end']}" for a in annotations if a['name'] == name]
            if ranges: color_string_parts.append(f"{color}: " + ", ".join(ranges))
            legend_items.append(html.Span([html.Span(style={'backgroundColor': color.replace('0.4', '1'), 'width': '15px', 'height': '15px', 'display': 'inline-block', 'marginRight': '5px', 'borderRadius': '50%'}), name], className='me-3'))
        return ", ".join(color_string_parts), html.Div(legend_items, className="mt-2 d-flex flex-wrap") if legend_items else None, name_to_color
    def build_annotated_linear_view(sequence, annotations, name_to_color):
        if not annotations: return html.Div(sequence, className="linear-sequence-view")
        
        spans = []
        last_idx = 0
        sorted_annos = sorted(annotations, key=lambda x: x['start'])
        
        for anno in sorted_annos:
            start, end, name = anno['start'] - 1, anno['end'], anno['name']
            if start > last_idx: spans.append(sequence[last_idx:start])
            color = name_to_color.get(name)
            spans.append(html.Span(sequence[start:end], style={'backgroundColor': color, 'borderRadius': '3px', 'padding': '2px 0'}))
            last_idx = end
        if last_idx < len(sequence): spans.append(sequence[last_idx:])
        
        return html.Div(spans, className="linear-sequence-view")
    def build_description_card(title, text):
        return [html.Details([html.Summary(title, className="h6"),html.P(text, className="text-muted mt-2 small")], className="mb-3")]
    def build_fd_figure(data):
        wt, mut = data.get('wt'), data.get('mut'); fig = go.Figure()
        fig.add_trace(go.Scatter(x=wt['fd_unfolding']['x'], y=wt['fd_unfolding']['y'], mode='lines', name='WT Unfolding', line=dict(color='royalblue')))
        if mut and not mut.get('error'): fig.add_trace(go.Scatter(x=mut['fd_unfolding']['x'], y=mut['fd_unfolding']['y'], mode='lines', name='Mutant', line=dict(color='skyblue', dash='dash')))
        fig.update_layout(title_text='Simulierte FD-Kurve', template="plotly_white", plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)')
        return fig
    def build_heatmap_figure(wt):
        df = pd.DataFrame(wt['mutation_scan_results']); bases = ['A', 'C', 'G', 'U']
        fig = go.Figure(data=go.Heatmap(z=df.pivot(index='Mutation', columns='Position', values='dMFE').loc[[f"{b_orig} > {b_new}" for b_orig in df['Mutation'].apply(lambda x: x.split(' > ')[0]).unique() for b_new in bases]].values, x=df['Position'].unique(), y=bases, colorscale='RdBu_r', zmid=0))
        fig.update_layout(title_text='Mutations-Sensitivitäts-Heatmap (ΔMFE)', template="plotly_white", plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)')
        return fig
    def build_dot_plot_figure(data):
        wt, mut = data.get('wt'), data.get('mut'); target = mut if mut and not mut.get('error') else wt
        fig = go.Figure(data=go.Heatmap(z=np.sqrt(np.array(target['base_pairing_probabilities'])), colorscale='Blues', showscale=False))
        fig.update_layout(title_text='Dot-Plot', yaxis_autorange='reversed', template="plotly_white", plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)')
        return fig
    def build_interactive_analysis_tab(fd_fig, wt):
        color_string, legend, name_to_color = get_annotation_colors(wt.get('annotations', []))
        sequences = [{'sequence': wt['sequence'].replace('T', 'U'), 'structure': wt['secondary_structure'], 'options': {'applyForce': False, 'circularizeExternal': True, 'labelInterval': 10, 'customColors': {'domain': [0, 1], 'range': ['#ffffff', '#ffffff'], 'colorValues': {'MyColoring': color_string}, 'colorAccessor': 'MyColoring'}}}]
        linear_view = build_annotated_linear_view(wt['sequence'], wt.get('annotations', []), name_to_color)
        return dcc.Tab(label='Interaktive Analyse', children=html.Div(className='p-3', children=
            build_description_card("Was sehe ich hier?", "Diese Ansicht koppelt die mechanische Stabilität (oben) mit der 2D-Faltung. Ein Tipp auf die FD-Kurve aktualisiert die 2D-Struktur und zeigt den entfalteten Zustand. Darunter sehen Sie Ihre Annotationen in der 2D- und der linearen Sequenz-Ansicht.") +
            [dcc.Graph(id='fd-curve-graph', figure=fd_fig), dashbio.FornaContainer(id='forna-container-interactive', sequences=sequences), html.Hr(), linear_view, legend or '']))
    def build_comparison_tab(data):
        wt, mut = data.get('wt'), data.get('mut')
        color_string_wt, legend_wt, _ = get_annotation_colors(wt.get('annotations', []))
        sequences_wt = [{'sequence': wt['sequence'].replace('T', 'U'), 'structure': wt['secondary_structure'], 'options': {'applyForce': False, 'circularizeExternal': True, 'labelInterval': 10, 'customColors': {'domain': [0, 1], 'range': ['#ffffff', '#ffffff'], 'colorValues': {'MyColoring': color_string_wt}, 'colorAccessor': 'MyColoring'}}}]
        content = [html.H4("Wildtyp Struktur", className="h6 text-center mb-2"), dashbio.FornaContainer(sequences=sequences_wt, height=500), legend_wt or '']
        if mut and not mut.get('error'):
            color_string_mut, legend_mut, _ = get_annotation_colors(mut.get('annotations', []))
            sequences_mut = [{'sequence': mut['sequence'].replace('T', 'U'), 'structure': mut['secondary_structure'], 'options': {'applyForce': False, 'circularizeExternal': True, 'labelInterval': 10, 'customColors': {'domain': [0, 1], 'range': ['#ffffff', '#ffffff'], 'colorValues': {'MyColoring': color_string_mut}, 'colorAccessor': 'MyColoring'}}}]
            content = html.Div(className='row', children=[html.Div(className='col-lg-6', children=[html.H4("Wildtyp", className="h6 text-center"), dashbio.FornaContainer(sequences=sequences_wt, height=500), legend_wt or '']), html.Div(className='col-lg-6', children=[html.H4("Mutante", className="h6 text-center"), dashbio.FornaContainer(sequences=sequences_mut, height=500), legend_mut or ''])])
        return dcc.Tab(label='2D-Struktur Vergleich', children=html.Div(className='p-3', children=content))
    def build_other_analyses_tab(data, dot_plot_fig):
        wt, mut = data.get('wt'), data.get('mut'); target = mut if mut and not mut.get('error') else wt
        summary = [html.H4("Zentrale Kennzahlen", className="h6"), html.P(f"Minimale freie Energie (MFE): {wt['mfe_kcal_mol']} kcal/mol. Dies ist ein Maß für die Faltungsstabilität.")]
        if mut and not mut.get('error'): summary.extend([html.Hr(), html.H4("Mutante", className="h6"), html.P(f"MFE: {mut['mfe_kcal_mol']} kcal/mol"), html.P(f"ΔMFE: {mut['mfe_kcal_mol'] - wt['mfe_kcal_mol']:.2f} kcal/mol", className="fw-bold")])
        return dcc.Tab(label='Weitere Analysen', children=dcc.Tabs(className='p-3', children=[
            dcc.Tab(label='Zusammenfassung', children=html.Div(className='pt-3', children=summary)),
            dcc.Tab(label='Dot-Plot', children=html.Div(className='pt-3', children=
                build_description_card("Was sehe ich hier?", "Der Dot-Plot visualisiert die Wahrscheinlichkeit, mit der jede Base (x-Achse) mit jeder anderen Base (y-Achse) paart. Dunkle Quadrate entlang der Diagonale repräsentieren hochwahrscheinliche und stabile Helices (Stems).") +
                [dcc.Graph(figure=dot_plot_fig)])),
            dcc.Tab(label='Dot-Bracket', children=html.Div(className='dot-bracket-display mt-3', children=[html.Code(f"{target['sequence']}\n{target['secondary_structure']}")]))
            ]))
    return server