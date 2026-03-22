import numpy as np
import pandas as pd
import csv
import math
from matplotlib.backends.backend_pdf import PdfPages

def count_motif(sequence: str, motif: str) -> int:
    sequence = sequence.upper()
    motif = motif.upper()

    if not motif or len(motif) > len(sequence):
        return 0

    m = len(motif)
    return sum(1 for i in range(len(sequence) - m + 1) if sequence[i:i+m] == motif)

def get_positions_single(sequence: str, motif: str):
    sequence = sequence.upper()
    motif = motif.upper()

    if not motif or len(motif) > len(sequence):
        return np.array([], dtype=int)

    positions = []
    start = 0
    while True:
        pos = sequence.find(motif, start)
        if pos == -1:
            break
        positions.append(pos)
        start = pos + 1

    return np.array(positions, dtype=int)
def build_counts_dataframe(sequences, motifs):
    motifs = [m.upper() for m in motifs]

    data = {}

    for seq in sequences:
        name = seq["name"]
        sequence = seq["sequence"]

        data[name] = [
            count_motif(sequence, m)
            for m in motifs
        ]

    df = pd.DataFrame(data, index=motifs)
    df.index.name = "Motywy"

    return df

def count_motifs_for_sequence(sequence, motifs):
    sequence = sequence.upper()
    return [count_motif(sequence, m) for m in motifs]

def get_motif_positions(sequences, motifs):
    motifs = [m.upper() for m in motifs]

    data = {}

    for seq in sequences:
        name = seq["name"]
        sequence = seq["sequence"].upper()

        data[name] = {}

        for motif in motifs:
            data[name][motif] = get_positions_single(sequence, motif).tolist()

    return data

def build_export_data(app, motifs):
    motifs = [m.upper() for m in motifs]

    csv_data = []
    headers = ["Motywy"] + [seq["name"] for seq in app.sequences]
    csv_data.append(headers)

    for motif in motifs:
        row = [motif]
        for seq in app.sequences:
            row.append(app.analysis_cache[seq["name"]][motif]["count"])
        csv_data.append(row)

    csv_data.append([])
    csv_data.append([])
    csv_data.append(["POZYCJE MOTYWÓW"])

    for seq in app.sequences:
        seq_name = seq["name"]
        csv_data.append([f"=== {seq_name} ==="])

        for motif in motifs:
            positions = app.analysis_cache[seq_name][motif]["positions"]
            csv_data.append([f"{motif}: {positions}"])

        csv_data.append([])

    return csv_data, app.export_figures

def save_csv(path, data):
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f, delimiter=";")
        writer.writerows(data)

def save_pdf(path, figures):
    with PdfPages(path) as pdf:
        for fig in figures:
            pdf.savefig(fig)

def build_heatmap_matrix_from_positions(sequence, motifs, motif_positions, segment_size=None):
    sequence = sequence.upper()
    motifs = [m.upper() for m in motifs]

    seq_len = len(sequence)

    if segment_size is None:
        target_segments = 25
        segment_size = max(10, seq_len // target_segments)

    n_segments = int(np.ceil(seq_len / segment_size))
    if n_segments > 50:
        segment_size = max(1, int(np.ceil(seq_len / 50)))
        n_segments = int(np.ceil(seq_len / segment_size))

    heatmap = np.zeros((len(motifs), n_segments), dtype=float)

    for m_idx, motif in enumerate(motifs):
        for pos in motif_positions[motif]:
            seg_idx = min(pos // segment_size, n_segments - 1)
            heatmap[m_idx, seg_idx] += 1

    for seg_idx in range(n_segments):
        start = seg_idx * segment_size
        end = min(start + segment_size, seq_len)
        seg_len = max(end - start, 1)
        heatmap[:, seg_idx] /= seg_len

    return heatmap