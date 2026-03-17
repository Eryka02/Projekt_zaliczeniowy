import sys
import json
import os
import ssl

# WYŁĄCZENIE SSL dla NCBI (tylko testy)
ssl._create_default_https_context = ssl._create_unverified_context

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QTextEdit, QListWidget, QTabWidget,
    QFileDialog, QLabel, QInputDialog, QMessageBox, QTableWidgetItem
)
from PyQt6.QtGui import (
    QAction, QBrush, QColor, QPen
)

from PyQt6.QtCore import (Qt, QRectF)
from PyQt6.QtWidgets import (
    QListWidgetItem, QGraphicsView, QGraphicsScene, QTableWidget
)

from Bio import Entrez
Entrez.email = "erykasworczuk993@gmail.com"


RECENT_FILE_STORAGE = "recent_files.json"
RECENT_NCBI_STORAGE = "recent_ncbi.json"

class DNAAnalyzer(QMainWindow):

    def __init__(self):
        super().__init__()

        self.setWindowTitle("Analiza motywów DNA")
        self.resize(1000, 650)

        self.sequences = []
        self.motifs = []

        self.recent_files = []
        self.recent_ncbi = []
        self.max_recent = 5

        self.init_ui()
        self.load_recent_files()
        self.load_recent_ncbi()

    def fetch_from_ncbi(self):
        # zapytanie od użytkownika
        query, ok = QInputDialog.getText(self, "Pobierz z NCBI", "Wpisz nazwę genu lub ID:")

        if not ok or not query:
            return

        try:
            self.log(f"Szukam sekwencji dla: {query}...")

            # wyszukiwanie w bazie nucleotide
            handle = Entrez.esearch(db="nucleotide", term=query, retmax=1)
            record = Entrez.read(handle)
            handle.close()

            if not record["IdList"]:
                self.log(f"Nie znaleziono sekwencji dla '{query}'")
                QMessageBox.warning(self, "Brak wyników", f"Nie znaleziono sekwencji dla '{query}'")
                return

            seq_id = record["IdList"][0]

            # pobranie sekwencji w formacie fasta
            handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
            fasta_data = handle.read()
            handle.close()

            # zapis do zmiennej sequence
            lines = fasta_data.splitlines()

            sequence = "".join([l.strip() for l in lines if not l.startswith(">")])

            self.sequences.append({
                "name": query,
                "sequence": sequence
            })

            self.refresh_sequence_view()

            self.log(f"Pobrano sekwencję '{query}' z NCBI (ID: {seq_id})")

            # zapisanie do historii NCBI
            ncbi_entry = {
                "query": query,
                "id": seq_id,
                "sequence": sequence
            }
            self.update_recent_ncbi(ncbi_entry)

        except Exception as e:
            self.log(f"Błąd pobierania z NCBI: {e}")
            QMessageBox.critical(self, "Błąd", f"Nie udało się pobrać sekwencji z NCBI:\n{str(e)}")

    def show_recent_ncbi(self):
        """Wyświetla okno dialogowe z ostatnimi pobraniami z NCBI"""
        if not self.recent_ncbi:
            QMessageBox.information(self, "Ostatnie pobrania NCBI", "Brak historii pobierań z NCBI")
            return

        items = [f"{entry['query']} (ID: {entry['id']})" for entry in self.recent_ncbi]

        item, ok = QInputDialog.getItem(
            self,
            "Ostatnie pobrania NCBI",
            "Wybierz sekwencję do wczytania:",
            items,
            0,
            False
        )

        if ok and item:
            index = items.index(item)
            entry = self.recent_ncbi[index]

            self.sequence = entry["sequence"]
            self.sequence_view.setText(self.sequence)
            self.log(f"Wczytano sekwencję z historii NCBI: {entry['query']}")

    def update_recent_ncbi(self, entry):
        """Aktualizuje listę ostatnich pobrań z NCBI"""
        self.recent_ncbi = [e for e in self.recent_ncbi if e["id"] != entry["id"]]
        self.recent_ncbi.insert(0, entry)

        if len(self.recent_ncbi) > self.max_recent:
            self.recent_ncbi = self.recent_ncbi[:self.max_recent]

        self.save_recent_ncbi()

    def save_recent_ncbi(self):
        """Zapisuje historię NCBI do pliku JSON"""
        try:
            with open(RECENT_NCBI_STORAGE, "w") as f:
                json.dump(self.recent_ncbi, f, indent=2)
        except Exception as e:
            self.log(f"Błąd zapisu historii NCBI: {e}")

    def load_recent_ncbi(self):
        """Wczytuje historię NCBI z pliku JSON"""
        if os.path.exists(RECENT_NCBI_STORAGE):
            try:
                with open(RECENT_NCBI_STORAGE) as f:
                    self.recent_ncbi = json.load(f)
            except Exception as e:
                self.log(f"Błąd wczytywania historii NCBI: {e}")
                self.recent_ncbi = []
    # =========================
    # GUI
    # =========================

    def init_ui(self):

        # ===== MENU =====

        menu = self.menuBar()

        self.file_menu = menu.addMenu("Plik")
        menu.addMenu("Motywy")
        self.ncbi_menu = menu.addMenu("NCBI")
        menu.addMenu("Eksport")
        menu.addMenu("Pomoc")

        open_action = QAction("Otwórz FASTA", self)
        open_action.triggered.connect(self.load_file)
        self.file_menu.addAction(open_action)

        self.file_menu.addSeparator()

        self.recent_menu = self.file_menu.addMenu("Ostatnie pliki")

        self.file_menu.addSeparator()

        exit_action = QAction("Wyjście", self)
        exit_action.triggered.connect(self.close)
        self.file_menu.addAction(exit_action)

        fetch_action = QAction("Pobierz sekwencję z NCBI", self)
        fetch_action.triggered.connect(self.fetch_from_ncbi)  # podłączamy funkcję
        self.ncbi_menu.addAction(fetch_action)

        recent_ncbi_action = QAction("Ostatnie pobrania NCBI", self)
        recent_ncbi_action.triggered.connect(self.show_recent_ncbi)
        self.ncbi_menu.addAction(recent_ncbi_action)

        # ===== LAYOUT GŁÓWNY =====

        main_widget = QWidget()
        main_layout = QVBoxLayout()

        content_layout = QHBoxLayout()

        # =========================
        # PANEL BOCZNY
        # =========================

        sidebar = QVBoxLayout()

        btn_load = QPushButton("Wczytaj plik")
        btn_ncbi = QPushButton("Pobierz z NCBI")
        btn_ncbi.clicked.connect(self.fetch_from_ncbi)
        btn_add_motif = QPushButton("Dodaj motyw")
        btn_run = QPushButton("Uruchom analizę")
        btn_export = QPushButton("Eksportuj CSV/PDF")

        btn_load.clicked.connect(self.load_file)
        btn_add_motif.clicked.connect(self.add_motif)
        btn_run.clicked.connect(self.run_analysis)

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

        # =========================
        # ZAKŁADKI
        # =========================

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
        # ===== WYNIKI =====

        results_widget = QWidget()
        results_layout = QVBoxLayout()

        # tabela wyników
        self.results_table = QTableWidget()
        results_layout.addWidget(self.results_table)

        # wizualizacja genomu
        self.visual_scene = QGraphicsScene()
        self.visual_view = QGraphicsView(self.visual_scene)
        self.visual_view.setMinimumHeight(220)

        results_layout.addWidget(self.visual_view)

        # lista pozycji motywów
        self.positions_view = QTextEdit()
        self.positions_view.setReadOnly(True)

        results_layout.addWidget(self.positions_view)

        results_widget.setLayout(results_layout)

        self.tabs.addTab(results_widget, "Wyniki analizy")

        # Wizualizacja
        self.visual_view = QTextEdit()
        self.visual_view.setText("Tutaj można dodać wizualizację genomu")
        self.tabs.addTab(self.visual_view, "Wizualizacja")

        # Eksport
        self.export_view = QTextEdit()
        self.export_view.setText("Panel eksportu wyników")
        self.tabs.addTab(self.export_view, "Eksport")

        # =========================
        # LOGI
        # =========================

        self.logs = QTextEdit()
        self.logs.setReadOnly(True)
        self.logs.setMaximumHeight(120)

        # =========================
        # SKŁADANIE GUI
        # =========================

        content_layout.addWidget(sidebar_widget, 1)
        content_layout.addWidget(self.tabs, 4)

        main_layout.addLayout(content_layout)
        main_layout.addWidget(self.logs)

        main_widget.setLayout(main_layout)

        self.setCentralWidget(main_widget)

    # =========================
    # LOG
    # =========================

    def log(self, text):
        self.logs.append("> " + text)

    # =========================
    # WCZYTYWANIE FASTA
    # =========================

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
            with open(path) as f:
                lines = f.readlines()

            sequence = "".join(
                [l.strip() for l in lines if not l.startswith(">")]
            )

            self.sequences.append({
                "name": os.path.basename(path),
                "sequence": sequence
            })

            self.refresh_sequence_view()

            self.log(f"Plik wczytany: {path}")

            self.update_recent_files(path)

        except Exception:
            self.log(f"Błąd wczytywania pliku: {path}")

    def refresh_sequence_view(self):

        text = ""

        for seq in self.sequences:
            text += f">{seq['name']}\n"
            text += seq["sequence"] + "\n\n"

        self.sequence_view.setText(text)
    # =========================
    # MOTYWY
    # =========================

    def add_motif(self):

        motif, ok = QInputDialog.getText(self, "Dodaj motyw", "Motyw DNA:")

        if ok and motif:
            motif = motif.upper()

            item = QListWidgetItem(motif)

            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)

            item.setCheckState(Qt.CheckState.Checked)

            self.motif_list.addItem(item)

            self.log(f"Motyw dodany: {motif}")

    # =========================
    # ANALIZA DNA
    # =========================

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

                count = seq["sequence"].count(motif)

                item = QTableWidgetItem(str(count))

                if count > 0:
                    item.setBackground(QColor("lightgreen"))

                self.results_table.setItem(row, col + 1, item)

        self.results_table.resizeColumnsToContents()

        self.draw_genome_map(motifs)

        pos_text = self.get_motif_positions(motifs)

        self.positions_view.setText(pos_text)

        self.log("Analiza zakończona")

    def draw_genome_map(self, motifs):

        self.visual_scene.clear()

        y = 40
        scale = 0.3

        colors = [
            QColor("yellow"),
            QColor("orange"),
            QColor("cyan"),
            QColor("pink"),
            QColor("green")
        ]

        for seq in self.sequences:

            sequence = seq["sequence"]
            name = seq["name"]

            length = len(sequence)

            self.visual_scene.addText(name).setPos(10, y - 20)

            genome_line = QRectF(120, y, length * scale, 6)

            self.visual_scene.addRect(genome_line, QPen(), QBrush(QColor("lightgray")))

            for m_index, motif in enumerate(motifs):

                start = 0

                while True:

                    pos = sequence.find(motif, start)

                    if pos == -1:
                        break

                    x = 120 + pos * scale
                    w = len(motif) * scale

                    rect = QRectF(x, y - 5, w, 16)

                    color = colors[m_index % len(colors)]

                    self.visual_scene.addRect(rect, QPen(), QBrush(color))

                    start = pos + 1

            y += 60

    def get_motif_positions(self, motifs):

        text = ""

        for seq in self.sequences:

            sequence = seq["sequence"]
            name = seq["name"]

            text += f"{name}\n"

            for motif in motifs:

                positions = []

                start = 0

                while True:

                    pos = sequence.find(motif, start)

                    if pos == -1:
                        break

                    positions.append(pos)

                    start = pos + 1

                if positions:
                    pos_text = ", ".join(map(str, positions))
                else:
                    pos_text = "brak"

                text += f"{motif} : {pos_text}\n"

            text += "\n"

        return text

    # =========================
    # RECENT FILES
    # =========================

    def update_recent_files(self, path):

        if path in self.recent_files:
            self.recent_files.remove(path)

        self.recent_files.insert(0, path)

        if len(self.recent_files) > self.max_recent:
            self.recent_files.pop()

        self.refresh_recent_list()
        self.refresh_recent_menu()
        self.save_recent_files()

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

    def open_recent_file(self, item):

        path = item.text()
        self.open_file(path)

    # =========================
    # JSON STORAGE
    # =========================

    def save_recent_files(self):

        try:
            with open(RECENT_FILE_STORAGE, "w") as f:
                json.dump(self.recent_files, f)

        except Exception:
            pass

    def load_recent_files(self):

        if os.path.exists(RECENT_FILE_STORAGE):

            try:
                with open(RECENT_FILE_STORAGE) as f:
                    self.recent_files = json.load(f)

            except Exception:
                self.recent_files = []

        self.refresh_recent_list()
        self.refresh_recent_menu()


# =========================
# START PROGRAMU
# =========================

app = QApplication(sys.argv)

window = DNAAnalyzer()
window.show()

sys.exit(app.exec())