import csv
from matplotlib.backends.backend_pdf import PdfPages

def count_motif(sequence: str, motif: str) -> int:
    sequence = sequence.upper()
    motif = motif.upper()

    count = 0
    start = 0

    while True:
        pos = sequence.find(motif, start)
        if pos == -1:
            break
        count += 1
        start = pos + 1  # pozwala na nakładanie się motywów

    return count
def find_positions(sequences, motifs) -> str:
    result = []

    for seq in sequences:
        name = seq["name"]
        sequence = seq["sequence"].upper()

        result.append(f"=== {name} ===")

        for motif in motifs:
            motif = motif.upper()
            positions = []

            start = 0
            while True:
                pos = sequence.find(motif, start)
                if pos == -1:
                    break
                positions.append(pos)
                start = pos + 1

            result.append(f"{motif}: {positions}")

        result.append("")

    return "\n".join(result)

def get_motif_positions(sequences, motifs):
    data = {}

    for seq in sequences:
        name = seq["name"]
        sequence = seq["sequence"].upper()

        data[name] = {}

        for motif in motifs:
            motif = motif.upper()
            positions = []

            start = 0
            while True:
                pos = sequence.find(motif, start)
                if pos == -1:
                    break
                positions.append(pos)
                start = pos + 1

            data[name][motif] = positions

    return data

def build_export_data(app, motifs):
    """
    Eksport:
    - tabela wyników (motywy vs sekwencje)
    - tabela pozycji motywów
    - figura matplotlib
    """

    csv_data = []

    # ===== 1. TABELA LICZB WYSTĄPIEŃ =====
    headers = ["Motywy"]
    for seq in app.sequences:
        headers.append(seq["name"])

    csv_data.append(headers)

    for row in range(app.results_table.rowCount()):
        row_data = []
        for col in range(app.results_table.columnCount()):
            item = app.results_table.item(row, col)
            row_data.append(item.text() if item else "")
        csv_data.append(row_data)

    # ===== odstęp =====
    csv_data.append([])
    csv_data.append([])

    # ===== 2. POZYCJE MOTYWÓW =====

    positions = get_motif_positions(app.sequences, motifs)

    csv_data.append(["POZYCJE MOTYWÓW"])

    for seq_name, motif_dict in positions.items():
        csv_data.append([f"=== {seq_name} ==="])

        for motif, pos_list in motif_dict.items():
            csv_data.append([f"{motif}: {pos_list}"])

        csv_data.append([])

    # ===== figura =====
    fig = app.visual_figure

    return csv_data, fig

def save_csv(path, data):
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f, delimiter=";")

        for row in data:
            writer.writerow(row)

def save_pdf(path, figure):
    try:
        with PdfPages(path) as pdf:
            pdf.savefig(figure)
    except Exception as e:
        print(f"Błąd zapisu PDF: {e}")
        raise

