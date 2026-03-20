import sys
from PyQt6.QtWidgets import QApplication
from ui.main_window import DNAAnalyzer

app = QApplication(sys.argv)

window = DNAAnalyzer()
window.show()

sys.exit(app.exec())