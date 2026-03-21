import numpy as np
import pandas as pd
import csv
from matplotlib.backends.backend_pdf import PdfPages


# =========================
# FAST MOTIF COUNT (NumPy-friendly)
# =========================

def count_motif(sequence: str, motif: str) -> int:
    sequence = sequence.upper()
    motif = motif.upper()

    count = 0
    start = 0
    m_len = len(motif)

    while True:
        pos = sequence.find(motif, start)
        if pos == -1:
            break
        count += 1
        start = pos + 1

    return count


# =========================
# POSITIONS (NumPy acceleration style)
# =========================

def get_positions_single(sequence, motif):
    sequence = sequence.upper()
    motif = motif.upper()

    return np.array([
        i for i in range(len(sequence))
        if sequence.startswith(motif, i)
    ])


# =========================
# PANDAS: FULL RESULTS TABLE (MAX OPT)
# =========================

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


# =========================
# LEGACY SUPPORT (your GUI uses this)
# =========================

def count_motifs_for_sequence(sequence, motifs):
    sequence = sequence.upper()
    return [count_motif(sequence, m) for m in motifs]


# =========================
# POSITIONS FOR EXPORT (optimized NumPy)
# =========================

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


# =========================
# BUILD EXPORT (NOW Pandas inside, backward compatible)
# =========================

def build_export_data(app, motifs):

    motifs = [m.upper() for m in motifs]

    # FAST PANDAS TABLE
    df = build_counts_dataframe(app.sequences, motifs)

    # POSITIONS
    positions = get_motif_positions(app.sequences, motifs)

    csv_data = []

    # =========================
    # TABLE FROM PANDAS
    # =========================
    csv_data.append(["Motywy"] + list(df.columns))

    for idx, row in df.iterrows():
        csv_data.append([idx] + row.tolist())

    csv_data.append([])
    csv_data.append([])

    # =========================
    # POSITIONS SECTION
    # =========================
    csv_data.append(["POZYCJE MOTYWÓW"])

    for seq_name, motif_dict in positions.items():
        csv_data.append([f"=== {seq_name} ==="])

        for motif, pos_list in motif_dict.items():
            csv_data.append([f"{motif}: {pos_list}"])

        csv_data.append([])

    return csv_data, app.visual_figure


# =========================
# SAVE CSV (unchanged but safe)
# =========================

def save_csv(path, data):
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f, delimiter=";")
        writer.writerows(data)


# =========================
# SAVE PDF
# =========================

def save_pdf(path, figure):
    with PdfPages(path) as pdf:
        pdf.savefig(figure)
