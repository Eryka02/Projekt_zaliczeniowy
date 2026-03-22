import sys

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QTextEdit, QListWidget, QTabWidget, QFileDialog,
    QLabel, QInputDialog, QMessageBox, QTableWidget, QTableWidgetItem,
    QListWidgetItem, QGraphicsView, QGraphicsScene, QGraphicsTextItem,
    QScrollArea
)
from PyQt6.QtGui import QAction, QBrush, QColor, QPen, QPainter
from PyQt6.QtCore import Qt, QRectF, QTimer

from pathlib import Path

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from core.analysis import (
    build_export_data,
    save_csv,
    save_pdf,
    build_heatmap_matrix_from_positions
)
from core.fasta import parse_fasta_file
from core.ncbi import fetch_sequence_from_ncbi
from utils.json_utils import load_json, save_json
from services.analysis_service import AnalysisService



class DNAAnalyzer(QMainWindow):

    VIEW_HEIGHT = 220
    LEFT_MARGIN = 120
    SCALE = 20
    ROW_HEIGHT = 40

    def __init__(self):
        super().__init__()

        self.setWindowTitle("Analiza motywów DNA")
        self.resize(1000, 650)

        self.sequences = []
        self.motifs = []
        self.analysis_cache = {}

        self.recent_files = load_json("recent.json", [])
        self.recent_ncbi = load_json("recent_ncbi.json", [])
        self.max_recent = 5

        self.init_ui()

        self.refresh_recent_list()
        self.refresh_recent_menu()
        self.export_figures = []

        self.viewport_size = 80
        self.SCALE = 16

    # GUI
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
        self.visual_view = ZoomGraphicsView(self.visual_scene, self)
        self.visual_view.setMinimumHeight(self.VIEW_HEIGHT)
        self.visual_view.setMinimumHeight(self.VIEW_HEIGHT)
        results_layout.addWidget(self.visual_view)
        self.zoom_start = 0
        self.viewport_size = 150

        self.positions_view = QTextEdit()
        self.positions_view.setReadOnly(True)
        results_layout.addWidget(self.positions_view)

        results_widget.setLayout(results_layout)
        self.tabs.addTab(results_widget, "Wyniki analizy")

        # Wizualizacja
        self.scroll_area = QScrollArea()
        self.scroll_area.setWidgetResizable(True)

        self.visual_container = QWidget()
        self.visual_layout = QVBoxLayout()

        self.visual_container.setLayout(self.visual_layout)
        self.scroll_area.setWidget(self.visual_container)

        self.tabs.addTab(self.scroll_area, "Wizualizacja")

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



    # NCBI
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

    # LOGI

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
    def refresh_sequence_view(self):

        text = "\n\n".join(
            f">{seq['name']}\n{seq['sequence']}"
            for seq in self.sequences
        )

        self.sequence_view.setText(text)


    # MOTYWY
    def add_motif(self):
        motif, ok = QInputDialog.getText(self, "Dodaj motyw", "Motyw DNA:")

        if not ok or not motif:
            return

        motif = motif.strip().upper()

        if not motif:
            return

        if any(ch not in "ATGC" for ch in motif):
            QMessageBox.warning(self, "Błąd", "Motyw może zawierać tylko A, T, G, C")
            return

        existing = {self.motif_list.item(i).text() for i in range(self.motif_list.count())}
        if motif in existing:
            self.log(f"Motyw już istnieje: {motif}")
            return

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

    def get_selected_motifs(self):
        return [
            self.motif_list.item(i).text()
            for i in range(self.motif_list.count())
            if self.motif_list.item(i).checkState() == Qt.CheckState.Checked
        ]

    # ANALIZA DNA

    def run_analysis(self):

        self.export_figures = []

        if not self.sequences:
            self.log("Brak sekwencji")
            return

        motifs = self.get_selected_motifs()

        if not motifs:
            self.log("Brak motywów")
            return

        columns = ["Motyw"]
        for seq in self.sequences:
            columns.append(seq["name"])

        self.results_table.setColumnCount(len(columns))
        self.results_table.setHorizontalHeaderLabels(columns)
        self.results_table.setRowCount(len(motifs))

        service = AnalysisService(self.sequences)
        self.analysis_cache = service.analyze_all(motifs)

        for row_index, motif in enumerate(motifs):
            self.results_table.setItem(row_index, 0, QTableWidgetItem(motif))

            for col_index, seq in enumerate(self.sequences):
                seq_name = seq["name"]
                count = self.analysis_cache[seq_name][motif]["count"]

                item = QTableWidgetItem(str(count))

                if count > 0:
                    item.setBackground(QColor("lightgreen"))

                self.results_table.setItem(row_index, col_index + 1, item)

        self.results_table.resizeColumnsToContents()

        self.draw_genome_map(motifs)
        self.draw_combined_visualization(motifs)

        text = ""
        for seq_name, motif_dict in self.analysis_cache.items():
            text += f"=== {seq_name} ===\n"
            for motif, data in motif_dict.items():
                text += f"{motif}: {data['positions']}\n"
            text += "\n"

        self.positions_view.setText(text)

        self.log("Analiza zakończona")

    def draw_genome_map(self, motifs):

        self.visual_scene.clear()

        y_info = 5
        y_dna = 30
        y_motif = 70

        scale = self.SCALE
        viewport = self.viewport_size

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

            full = seq["sequence"]
            length = len(full)

            start = max(0, min(self.zoom_start, length - viewport))
            end = min(start + viewport, length)

            visible = full[start:end]

            self.visual_scene.addText(
                f"{seq['name']} | {start} - {end}"
            ).setPos(10, y_info)


            for i, n in enumerate(visible):
                x = self.LEFT_MARGIN + i * scale

                self.visual_scene.addRect(
                    QRectF(x, y_dna, scale, scale),
                    QPen(Qt.GlobalColor.black),
                    QBrush(nucleotide_colors.get(n, QColor(200, 200, 200)))
                )

                if scale >= 14:
                    text = QGraphicsTextItem(n)

                    rect = text.boundingRect()
                    text.setPos(
                        x + (scale - rect.width()) / 2,
                        y_dna + (scale - rect.height()) / 2 - 1
                    )

                    self.visual_scene.addItem(text)

            # MOTIFS
            for m_index, motif in enumerate(motifs):

                positions = self.analysis_cache[seq["name"]][motif]["positions"]

                for pos in positions:

                    if start <= pos < end:

                        local = pos - start
                        color = motif_colors[m_index % len(motif_colors)]

                        for j in range(len(motif)):
                            x = self.LEFT_MARGIN + (local + j) * scale

                            self.visual_scene.addRect(
                                QRectF(x, y_motif, scale, scale),
                                QPen(Qt.GlobalColor.black),
                                QBrush(color)
                            )

            y_info += self.ROW_HEIGHT + 30
            y_dna += self.ROW_HEIGHT + 30
            y_motif += self.ROW_HEIGHT + 30

    def draw_combined_visualization(self, motifs):

        self.export_figures = []

        # czyści stare wykresy
        while self.visual_layout.count():
            item = self.visual_layout.takeAt(0)
            widget = item.widget()
            if widget:
                widget.deleteLater()

        # BAR CHART
        fig_bar = Figure(figsize=(8, 4))
        canvas_bar = FigureCanvas(fig_bar)
        self.export_figures.append(fig_bar)
        canvas_bar.setMinimumHeight(400)
        ax = fig_bar.add_subplot(111)

        n_seq = len(self.sequences)
        n_motif = len(motifs)

        bar_width = 0.8 / n_seq
        x = range(n_motif)

        for i, seq in enumerate(self.sequences):
            counts = [
                self.analysis_cache[seq["name"]][motif]["count"]
                for motif in motifs
            ]
            positions = [xi + i * bar_width for xi in x]
            ax.bar(positions, counts, width=bar_width, label=seq["name"])

        ax.set_xticks([
            xi + bar_width * (n_seq / 2) - (bar_width / 2)
            for xi in x
        ])
        ax.set_xticklabels(motifs)
        ax.set_title("Liczba wystąpień motywów")
        ax.set_ylabel("Liczba")

        ax.legend()

        self.visual_layout.addWidget(canvas_bar)

        # HEATMAPY

        for seq in self.sequences:

            fig = Figure(figsize=(10, 4))
            canvas = FigureCanvas(fig)
            self.export_figures.append(fig)
            canvas.setMinimumHeight(400)
            ax = fig.add_subplot(111)

            motif_positions = {
                motif: self.analysis_cache[seq["name"]][motif]["positions"]
                for motif in motifs
            }
            matrix = build_heatmap_matrix_from_positions(seq["sequence"], motifs, motif_positions)

            im = ax.imshow(
                matrix,
                aspect="equal",
                interpolation="none",
                cmap="magma",
                rasterized=True
            )

            # Y = motywy
            ax.set_yticks(range(len(motifs)))
            ax.set_yticklabels(motifs)

            # X = segmenty
            segments = [f"S{j + 1}" for j in range(matrix.shape[1])]
            ax.set_xticks(range(len(segments)))

            if len(segments) > 20:
                ax.set_xticklabels(segments, rotation=90, fontsize=6)
            else:
                ax.set_xticklabels(segments, rotation=45)

            ax.set_title(f"Heatmapa: {seq['name']}")

            fig.colorbar(im, ax=ax)

            self.visual_layout.addWidget(canvas)

        self.visual_layout.addStretch()
        total_widgets = self.visual_layout.count()
        self.visual_container.setMinimumHeight(total_widgets * 420)

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

            motifs = self.get_selected_motifs()

            csv_data, figures = build_export_data(self, motifs)

            save_csv(path, csv_data)

            pdf_path = str(Path(path).with_suffix(".pdf"))
            save_pdf(pdf_path, figures)


            self.log(f"Zapisano CSV: {path}")
            self.log(f"Zapisano PDF: {pdf_path}")

        except Exception as e:
            self.log(f"Błąd eksportu: {e}")


class ZoomGraphicsView(QGraphicsView):

    def __init__(self, scene, parent=None):
        super().__init__(scene, parent)

        self.parent = parent

        self.setRenderHint(QPainter.RenderHint.Antialiasing)
        self.setDragMode(QGraphicsView.DragMode.NoDrag)

        self.last_x = None
        self.dragging = False

        # anti crash
        self._blocked = False


    def wheelEvent(self, event):
        pass

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self.dragging = True
            self.last_x = event.position().x()

    def mouseReleaseEvent(self, event):
        self.dragging = False
        self.last_x = None

    def mouseMoveEvent(self, event):

        if not self.dragging:
            return

        if self.last_x is None:
            self.last_x = event.position().x()

        delta = event.position().x() - self.last_x
        self.last_x = event.position().x()

        self.parent.zoom_start -= int(delta / 3)

        if self.parent.zoom_start < 0:
            self.parent.zoom_start = 0

        if self.parent.sequences:
            max_len = max(len(seq["sequence"]) for seq in self.parent.sequences)
            max_start = max(0, max_len - self.parent.viewport_size)

            if self.parent.zoom_start > max_start:
                self.parent.zoom_start = max_start

        self.schedule_redraw()

    def schedule_redraw(self):

        if self._blocked:
            return

        self._blocked = True

        QTimer.singleShot(60, self._do_redraw)

    def _do_redraw(self):

        self._blocked = False

        try:
            self.parent.draw_genome_map(
                self.parent.get_selected_motifs()
            )
        except Exception as e:
            print("REDRAW ERROR:", e)