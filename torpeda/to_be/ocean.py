from PyQt5.QtWidgets import *

from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.figure import Figure

class Ocean(QWidget):

    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        # Создаем холст для отрисовки
        self.canvas = FigureCanvas(Figure())
        vertical_layout = QVBoxLayout()
        vertical_layout.addWidget(self.canvas)
        # Добавляем одну ось для построения графиков
        self.canvas.axes = self.canvas.figure.add_subplot(1, 1, 1)
        self.setLayout(vertical_layout)
