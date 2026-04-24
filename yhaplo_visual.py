import pandas as pd
import numpy as np
import re
import plotly.graph_objects as go
import plotly.colors
import math
from pathlib import Path

# === 1. ГЛОБАЛЬНАЯ ПАЛИТРА ===
PALETTE = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
    '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5',
    '#c49c94', '#f7b6d2', '#c7c7c7', '#dbdb8d', '#9edae5'
]
COLOR_MAP = {}
COLOR_IDX = 0

def get_color(label):
    global COLOR_MAP, COLOR_IDX
    key = label.split('-')[0] if '-' in str(label) else str(label)
    if key not in COLOR_MAP:
        COLOR_MAP[key] = PALETTE[COLOR_IDX % len(PALETTE)]
        COLOR_IDX += 1
    return COLOR_MAP[key]

# === 2. ФУНКЦИИ ===
def parse_hg_hierarchy(hg: str, max_depth: int = 4):
    if pd.isna(hg) or hg == 'NA' or not isinstance(hg, str):
        return [None] * max_depth
    hg_clean = re.sub(r'\*\(x[^)]*\)', '*', hg).strip()
    hg_main = hg_clean.split('*')[0]
    parts = []
    tokens = re.split(r'[-_]', hg_main)
    for i in range(len(tokens)):
        parts.append('-'.join(tokens[:i+1]))
    while len(parts) < max_depth:
        parts.append(None)
    return parts[:max_depth]

def get_sunburst_trace(df_subset, pop_name):
    MAX_DEPTH = 4
    rows = []
    valid_males_count = len(df_subset)
    
    root_id = f"{pop_name}_root"
    rows.append({
        'id': root_id, 'parent': '', 'value': valid_males_count,
        'label': pop_name, 'is_root': True
    })
    
    for _, row in df_subset.iterrows():
        levels = parse_hg_hierarchy(row['Hg'], MAX_DEPTH)
        path = []
        for lvl in levels:
            if pd.isna(lvl): break
            path.append(lvl)
            node_id = '/'.join(path)
            parent = '/'.join(path[:-1]) if len(path) > 1 else root_id
            rows.append({
                'id': f"{pop_name}::{node_id}",
                'parent': f"{pop_name}::{parent}" if parent != root_id else root_id,
                'value': 1, 'label': lvl, 'is_root': False
            })
    
    agg = pd.DataFrame(rows).groupby(['id', 'parent', 'label', 'is_root'], as_index=False)['value'].sum()
    
    colors = []
    for idx, row in agg.iterrows():
        if row['is_root']:
            colors.append('#FFFFFF')
        else:
            colors.append(get_color(row['label']))
    
    return go.Sunburst(
        ids=agg['id'],
        parents=agg['parent'],
        labels=agg['label'],
        values=agg['value'],
        branchvalues='total',
        maxdepth=4,
        marker=dict(colors=colors),
        textinfo='label',
        insidetextfont=dict(size=13, color='white'),
        hovertemplate='<b>%{label}</b><br>Samples: %{value}<extra></extra>',
    )

# === 3. ЗАГРУЗКА ДАННЫХ ===
HG_FILE = 'hg_prediction.hg'
META_FILE = 'groups.csv'

df_hg = pd.read_csv(HG_FILE, sep='\t')
df_meta = pd.read_csv(META_FILE, sep=';')

df_merged = pd.merge(df_hg, df_meta, left_on='Sample_name', right_on='ID', how='inner')
df_males = df_merged[(df_merged['Hg'] != 'NA') & (df_merged['Hg'].notna())].copy()

populations = sorted(df_males['pop'].unique())
n_groups = len(populations)
print(f"Популяций: {n_groups}")
print(f"Образцов мужчин: {len(df_males)}")

# === 4. ПОСТРОЕНИЕ СЕТКИ ===
cols = math.ceil(math.sqrt(n_groups))
rows = math.ceil(n_groups / cols)

fig = go.Figure()
x_step = 1.0 / cols
y_step = 1.0 / rows

for i, pop in enumerate(populations):
    r = i // cols
    c = i % cols
    x0, x1 = c * x_step, (c + 1) * x_step
    y0, y1 = 1.0 - (r + 1) * y_step, 1.0 - r * y_step
    
    margin = 0.03  # ВАШ ОТСТУП
    domain = {'x': [x0 + margin, x1 - margin], 'y': [y0 + margin, y1 - margin]}
    
    df_pop_males = df_males[df_males['pop'] == pop]
    total_males = len(df_pop_males)
    
    trace = get_sunburst_trace(df_pop_males, pop)
    trace.domain = domain
    fig.add_trace(trace)
    
    # === ЦЕНТРАЛЬНАЯ ПОДПИСЬ ===
    fig.add_annotation(
        text=f"<b>{pop}</b><br><span style='font-size:12px'>M: {total_males}</span>",
        x=(x0 + x1) / 2,
        y=(y0 + y1) / 2,
        xref="paper",
        yref="paper",
        xanchor="center",
        yanchor="middle",
        showarrow=False,
        font=dict(size=14, color="black"),
        align="center",
        bgcolor="rgba(255,255,255,0.9)",
        bordercolor="gray",
        borderwidth=1,
        borderpad=5
    )

# === 5. МАКЕТ ===
fig.update_layout(
    title_text="Y-Chromosome Haplogroups by Population",
    height=400 * rows,
    width=400 * cols,
    showlegend=False,
    margin=dict(t=50, b=10, l=10, r=10)
)

# === 6. СОХРАНЕНИЕ ===
OUT_HTML = 'populations_grid_final.html'

fig.write_html(OUT_HTML, include_plotlyjs='cdn')
print(f"Файл: {OUT_HTML}")