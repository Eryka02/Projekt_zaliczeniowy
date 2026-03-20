import sys
import ssl
import csv

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QTextEdit, QListWidget, QTabWidget, QFileDialog,
    QLabel, QInputDialog, QMessageBox, QTableWidget, QTableWidgetItem,
    QListWidgetItem, QGraphicsView, QGraphicsScene, QGraphicsTextItem
)

from PyQt6.QtGui import QAction, QBrush, QColor, QPen
from PyQt6.QtCore import Qt, QRectF

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from core.analysis import (
    build_export_data,
    save_csv,
    save_pdf,
    count_motif,
    find_positions,
    get_motif_positions
)

from core.fasta import parse_fasta_file
from core.ncbi import fetch_sequence_from_ncbi

from utils.json_utils import load_json, save_json


class DNAAnalyzer(QMainWindow):

    # 1
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Analiza motywów DNA")
        self.resize(1000, 650)

        self.sequences = []
        self.motifs = []

        self.recent_files = load_json("recent.json", [])
        self.recent_ncbi = load_json("recent_ncbi.json", [])
        self.max_recent = 5

        self.init_ui()

        self.refresh_recent_list()
        self.refresh_recent_menu()


    #2 ncbi
    def fetch_from_ncbi(self):

        query, ok = QInputDialog.getText(
            self,
            "Pobierz z NCBI",
            "Wpisz nazwę genu lub ID:"
        )

        if not ok or not query:
            return


        self.log(f"Szukam sekwencji dla: {query}...")

        entry = fetch_sequence_from_ncbi(query)

        if not entry:
            self.log("Nie znaleziono sekwencji")
            QMessageBox.warning(self, "Brak wyników", "Nie znaleziono")
            return

        self.sequences.append({
            "name": entry["query"],
            "sequence": entry["sequence"]
        })

        self.refresh_sequence_view()

        self.update_recent_ncbi(entry)

        self.update_recent_files(
            f"NCBI: {query} ({entry['id']})"
        )

    def show_recent_ncbi(self):

        if not self.recent_ncbi:
            QMessageBox.information(
                self,
                "Ostatnie pobrania NCBI",
                "Brak historii"
            )
            return

        items = [f"{e['query']} (ID: {e['id']})" for e in self.recent_ncbi]

        item, ok = QInputDialog.getItem(
            self,
            "Ostatnie pobrania NCBI",
            "Wybierz:",
            items,
            0,
            False
        )

        if ok and item:
            index = items.index(item)
            entry = self.recent_ncbi[index]

            self.sequences.append({
                "name": entry["query"],
                "sequence": entry["sequence"]
            })

            self.refresh_sequence_view()
            self.log(f"Dodano z historii NCBI: {entry['query']}")

    def update_recent_ncbi(self, entry):

        self.recent_ncbi = [
            e for e in self.recent_ncbi if e["id"] != entry["id"]
        ]

        self.recent_ncbi.insert(0, entry)

        if len(self.recent_ncbi) > self.max_recent:
            self.recent_ncbi = self.recent_ncbi[:self.max_recent]

        save_json("recent_ncbi.json", self.recent_ncbi)


    # 3. GUI
    def init_ui(self):
        # MENU
        menu = self.menuBar()

        self.file_menu = menu.addMenu("Plik")
        self.motif_menu = menu.addMenu("Motywy")
        self.ncbi_menu = menu.addMenu("NCBI")
        self.export_menu = menu.addMenu("Eksport")
        menu.addMenu("Pomoc")

        #  Plik
        open_action = QAction("Otwórz FASTA", self)
        open_action.triggered.connect(self.load_file)
        self.file_menu.addAction(open_action)

        self.file_menu.addSeparator()
        self.recent_menu = self.file_menu.addMenu("Ostatnie pliki")
        self.file_menu.addSeparator()

        exit_action = QAction("Wyjście", self)
        exit_action.triggered.connect(self.close)
        self.file_menu.addAction(exit_action)

        # Motywy
        add_motif_action = QAction("Dodaj motyw", self)
        add_motif_action.triggered.connect(self.add_motif)

        remove_motif_action = QAction("Usuń zaznaczony motyw", self)
        remove_motif_action.triggered.connect(self.remove_selected_motif)

        self.motif_menu.addAction(add_motif_action)
        self.motif_menu.addAction(remove_motif_action)

        # NCBI
        fetch_action = QAction("Pobierz sekwencję z NCBI", self)
        fetch_action.triggered.connect(self.fetch_from_ncbi)
        self.ncbi_menu.addAction(fetch_action)

        recent_ncbi_action = QAction("Ostatnie pobrania NCBI", self)
        recent_ncbi_action.triggered.connect(self.show_recent_ncbi)
        self.ncbi_menu.addAction(recent_ncbi_action)

        # Eksport
        export_action = QAction("Eksportuj wyniki", self)
        export_action.triggered.connect(self.export_results)
        self.export_menu.addAction(export_action)

        # CENTRAL WIDGET

        main_widget = QWidget()
        main_layout = QVBoxLayout()
        content_layout = QHBoxLayout()

        # SIDEBAR

        sidebar = QVBoxLayout()

        btn_load = QPushButton("Wczytaj plik")
        btn_ncbi = QPushButton("Pobierz z NCBI")
        btn_add_motif = QPushButton("Dodaj motyw")
        btn_run = QPushButton("Uruchom analizę")
        btn_export = QPushButton("Eksportuj CSV/PDF")

        btn_load.clicked.connect(self.load_file)
        btn_ncbi.clicked.connect(self.fetch_from_ncbi)
        btn_add_motif.clicked.connect(self.add_motif)
        btn_run.clicked.connect(self.run_analysis)
        btn_export.clicked.connect(self.export_results)

        sidebar.addWidget(btn_load)
        sidebar.addWidget(btn_ncbi)
        sidebar.addWidget(btn_add_motif)
        sidebar.addWidget(btn_run)
        sidebar.addWidget(btn_export)

        sidebar.addWidget(QLabel("Ostatnio używane pliki"))

        self.recent_list = QListWidget()
        self.recent_list.itemClicked.connect(self.open_recent_file)
        sidebar.addWidget(self.recent_list)

        sidebar.addStretch()

        sidebar_widget = QWidget()
        sidebar_widget.setLayout(sidebar)

        # TABS

        self.tabs = QTabWidget()

        # Podgląd sekwencji
        self.sequence_view = QTextEdit()
        self.sequence_view.setReadOnly(True)
        self.tabs.addTab(self.sequence_view, "Podgląd sekwencji")

        # Motywy
        motif_widget = QWidget()
        motif_layout = QVBoxLayout()

        self.motif_list = QListWidget()

        motif_layout.addWidget(QLabel("Lista motywów"))
        motif_layout.addWidget(self.motif_list)

        motif_widget.setLayout(motif_layout)
        self.tabs.addTab(motif_widget, "Wybór motywów")

        # Wyniki
        results_widget = QWidget()
        results_layout = QVBoxLayout()

        self.results_table = QTableWidget()
        results_layout.addWidget(self.results_table)

        self.visual_scene = QGraphicsScene()
        self.visual_view = QGraphicsView(self.visual_scene)
        self.visual_view.setMinimumHeight(220)
        results_layout.addWidget(self.visual_view)

        self.positions_view = QTextEdit()
        self.positions_view.setReadOnly(True)
        results_layout.addWidget(self.positions_view)

        results_widget.setLayout(results_layout)
        self.tabs.addTab(results_widget, "Wyniki analizy")

        # Wizualizacja
        self.visual_figure = Figure(figsize=(6, 4))
        self.visual_canvas = FigureCanvas(self.visual_figure)

        visual_widget = QWidget()
        visual_layout = QVBoxLayout()
        visual_layout.addWidget(self.visual_canvas)

        visual_widget.setLayout(visual_layout)
        self.tabs.addTab(visual_widget, "Wizualizacja")

        # LOGI

        self.logs = QTextEdit()
        self.logs.setReadOnly(True)
        self.logs.setMaximumHeight(120)

        # SKŁADANIE UI

        content_layout.addWidget(sidebar_widget, 1)
        content_layout.addWidget(self.tabs, 4)

        main_layout.addLayout(content_layout)
        main_layout.addWidget(self.logs)

        main_widget.setLayout(main_layout)
        self.setCentralWidget(main_widget)

    # LOG

    def log(self, text):
        self.logs.append("> " + text)

    # WCZYTYWANIE FASTA

    def load_file(self):

        path, _ = QFileDialog.getOpenFileName(
            self,
            "Wybierz FASTA",
            "",
            "FASTA (*.fasta *.fa *.txt)"
        )

        if path:
            self.open_file(path)

    def open_file(self, path):

        try:

            name, sequence = parse_fasta_file(path)

            self.sequences.append({
                "name": name,
                "sequence": sequence
            })

            self.refresh_sequence_view()
            self.log(f"Plik wczytany: {path}")

            self.update_recent_files(path)

        except Exception as e:
            self.log(f"Błąd wczytywania pliku: {path} | {e}")

    def open_recent_file(self, item):

        text = item.text()

        try:
            # ===== NCBI =====
            if text.startswith("NCBI:"):

                if "(" not in text:
                    return

                seq_id = text.split("(")[-1].replace(")", "").strip()

                for entry in self.recent_ncbi:
                    if entry["id"] == seq_id:
                        self.sequences.append({
                            "name": entry["query"],
                            "sequence": entry["sequence"]
                        })

                        self.refresh_sequence_view()
                        self.log(f"Dodano sekwencję z sidebar NCBI: {entry['query']}")
                        return

            # PLIK
            self.open_file(text)

        except Exception as e:
            self.log(f"Błąd kliknięcia recent: {e}")

    def refresh_sequence_view(self):

        text = ""

        for seq in self.sequences:
            text += f">{seq['name']}\n"
            text += seq["sequence"] + "\n\n"

        self.sequence_view.setText(text)

    # MOTYWY
    def add_motif(self):

        motif, ok = QInputDialog.getText(self, "Dodaj motyw", "Motyw DNA:")

        if ok and motif:
            motif = motif.upper()

            item = QListWidgetItem(motif)

            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)

            item.setCheckState(Qt.CheckState.Checked)

            self.motif_list.addItem(item)

            self.log(f"Motyw dodany: {motif}")

    def remove_selected_motif(self):
        if self.motif_list.count() == 0:
            self.log("Brak motywów do usunięcia")
            return

        motifs = [self.motif_list.item(i).text() for i in range(self.motif_list.count())]

        motif, ok = QInputDialog.getItem(
            self,
            "Usuń motyw",
            "Wybierz motyw do usunięcia:",
            motifs,
            0,
            False
        )

        if ok and motif:
            for i in range(self.motif_list.count()):
                item = self.motif_list.item(i)
                if item.text() == motif:
                    self.motif_list.takeItem(i)
                    self.log(f"Usunięto motyw: {motif}")
                    break

    # EKSPORT
    def export_results(self):

        if self.results_table.rowCount() == 0:
            self.log("Brak danych do eksportu")
            return

        path, _ = QFileDialog.getSaveFileName(
            self,
            "Zapisz CSV",
            "",
            "CSV (*.csv)"
        )

        if not path:
            return

        try:
            from core.analysis import build_export_data, save_csv, save_pdf

            # ✔ MOTYWY zbierane w GUI (bez Qt w core!)
            motifs = [
                self.motif_list.item(i).text()
                for i in range(self.motif_list.count())
                if self.motif_list.item(i).checkState() == Qt.CheckState.Checked
            ]

            csv_data, pdf_fig = build_export_data(self, motifs)

            save_csv(path, csv_data)

            pdf_path = path.replace(".csv", ".pdf")
            save_pdf(pdf_path, pdf_fig)

            self.log(f"Zapisano CSV: {path}")
            self.log(f"Zapisano PDF: {pdf_path}")

        except Exception as e:
            self.log(f"Błąd eksportu: {e}")


    # ANALIZA DNA

    def run_analysis(self):

        if not self.sequences:
            self.log("Brak sekwencji")
            return

        motifs = []

        for i in range(self.motif_list.count()):
            item = self.motif_list.item(i)
            if item.checkState() == Qt.CheckState.Checked:
                motifs.append(item.text())

        if not motifs:
            self.log("Brak motywów")
            return

        columns = ["Motyw"]
        for seq in self.sequences:
            columns.append(seq["name"])

        self.results_table.setColumnCount(len(columns))
        self.results_table.setHorizontalHeaderLabels(columns)
        self.results_table.setRowCount(len(motifs))

        for row, motif in enumerate(motifs):

            self.results_table.setItem(row, 0, QTableWidgetItem(motif))

            for col, seq in enumerate(self.sequences):

                count = count_motif(seq["sequence"], motif)

                item = QTableWidgetItem(str(count))

                if count > 0:
                    item.setBackground(QColor("lightgreen"))

                self.results_table.setItem(row, col + 1, item)

        self.results_table.resizeColumnsToContents()

        self.draw_genome_map(motifs)
        self.draw_motif_bar_chart(motifs)

        pos_text = find_positions(self.sequences, motifs)

        self.positions_view.setText(pos_text)

        self.log("Analiza zakończona")

    def draw_genome_map(self, motifs):

        self.visual_scene.clear()

        y = 40
        scale = 20

        nucleotide_colors = {
            "A": QColor(173, 255, 47),
            "T": QColor(255, 182, 193),
            "G": QColor(255, 255, 153),
            "C": QColor(135, 206, 250)
        }

        motif_colors = [
            QColor(255, 160, 122),
            QColor(216, 191, 216),
            QColor(144, 238, 144),
            QColor(221, 160, 221),
            QColor(255, 228, 181)
        ]

        for seq in self.sequences:

            sequence = seq["sequence"]
            name = seq["name"]

            self.visual_scene.addText(name).setPos(10, y - 25)

            # rysowanie DNA
            for i, nucleotide in enumerate(sequence):
                x = 120 + i * scale
                rect = QRectF(x, y, scale, scale)

                color = nucleotide_colors.get(
                    nucleotide.upper(),
                    QColor(211, 211, 211)
                )

                self.visual_scene.addRect(rect, QPen(Qt.GlobalColor.black), QBrush(color))

                text_item = QGraphicsTextItem(nucleotide)
                text_item.setPos(x + 4, y)
                self.visual_scene.addItem(text_item)

            # nakładanie motywów (dalej UI, ale bez logiki find w core)
            for m_index, motif in enumerate(motifs):

                start = 0

                while True:
                    pos = sequence.find(motif, start)
                    if pos == -1:
                        break

                    motif_color = motif_colors[m_index % len(motif_colors)]

                    for j, nucleotide in enumerate(motif):
                        x = 120 + (pos + j) * scale
                        rect = QRectF(x, y, scale, scale)

                        self.visual_scene.addRect(
                            rect,
                            QPen(Qt.GlobalColor.black),
                            QBrush(motif_color)
                        )

                        text_item = QGraphicsTextItem(nucleotide)
                        text_item.setPos(x + 4, y)
                        self.visual_scene.addItem(text_item)

                    start = pos + 1

            y += 40


    def draw_motif_bar_chart(self, motifs):

        n_seq = len(self.sequences)
        self.visual_figure.clear()
        ax = self.visual_figure.add_subplot(111)

        n_motif = len(motifs)

        bar_width = 0.8 / n_seq
        x = range(n_motif)

        for i, seq in enumerate(self.sequences):
            counts = [seq["sequence"].count(m) for m in motifs]

            positions = [xi + i * bar_width for xi in x]

            ax.bar(positions, counts, width=bar_width, label=seq["name"])

        ax.set_xticks([xi + bar_width * (n_seq / 2) - (bar_width / 2) for xi in x])
        ax.set_xticklabels(motifs)

        ax.set_ylabel("Liczba wystąpień")
        ax.set_title("Wystąpienia motywów DNA w sekwencjach")
        ax.legend()

        self.visual_canvas.draw()
    # RECENT FILES

    def update_recent_files(self, path):

        if path in self.recent_files:
            self.recent_files.remove(path)

        self.recent_files.insert(0, path)

        if len(self.recent_files) > self.max_recent:
            self.recent_files = self.recent_files[:self.max_recent]

        self.refresh_recent_list()
        self.refresh_recent_menu()

        save_json("recent.json", self.recent_files)

    def refresh_recent_list(self):

        self.recent_list.clear()

        for path in self.recent_files:
            self.recent_list.addItem(path)

    def refresh_recent_menu(self):

        self.recent_menu.clear()

        for path in self.recent_files:

            action = QAction(path, self)

            action.triggered.connect(
                lambda checked=False, p=path: self.open_file(p)
            )

            self.recent_menu.addAction(action)






