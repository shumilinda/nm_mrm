import ctypes
import os
import sys
import warnings

from PyQt5 import uic
from PyQt5.QtWidgets import QApplication, QMainWindow, QTableWidgetItem, QWidget, QLabel, QVBoxLayout

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import numpy as np

warnings.filterwarnings('ignore')

ui_file = 'UI\MainWindow.ui'

def create_plot(parent):
    parent.fig = Figure(figsize = (parent.width() / 80, parent.height() / 80))
    parent.canvas = FigureCanvas(parent.fig)
    parent.plot = parent.fig.add_subplot(111, projection = "3d", position=[0.01, 0.25, 0.75, 0.75])
    return parent.plot

class SecondWindow(QMainWindow):
    def __init__(self, text):
        super().__init__()
        
        self.setWindowTitle("Результаты")
        
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        
        layout = QVBoxLayout(self.central_widget)
        self.label = QLabel(text)
        layout.addWidget(self.label)

class MainWindow(QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        uic.loadUi(ui_file, self)
        
        self.setWindowTitle("Лабораторная работа №2, метод минимальных невязок, Шумилин Денис 3821Б1ПМоп2")

        self.plt = create_plot(self.plot)
        self.plot.canvas.setParent(self.plot)

        self.plot_toolBar = NavigationToolbar(self.plot.canvas, self)
        self.addToolBar(self.plot_toolBar)

        self.b_plot.clicked.connect(self.plotting)
        self.b_delete.clicked.connect(self.clear_plots)

        self.plot.plot.set_xlabel("x")
        self.plot.plot.set_ylabel("y")
        self.plot.plot.set_zlabel("z")

    def clear_plots(self):
        self.clear_plot(self.plot)

    def clear_plot(self, cur_plot_widget):
        cur_plot_widget.plot.cla()
        cur_plot_widget.canvas.draw()
        
    def processing_results(self, filename):
        matrices = []
        
        with open(filename, 'r') as file:
            matrix = []
            lines = file.readlines()
            for line in lines:
                if line != '\n':
                    matrix.append(list(map(float, line.strip().split())))
                else:
                    matrices.append(matrix.copy())
                    matrix.clear()
            
        return matrices
    
    def get_informations(self, filename):
        matrices = self.processing_results(filename)
    
        if self.combo_task.currentText() == "Тестовая":

            res, diff, eps, S, time = map(lambda x: x[0], matrices[0])
            time_dS = round(float(time) / float(S), 6) 
        
            stop_reason = "превышения количества итераций" if res == 0.0 else "выхода на точность"
            text_to_display = (
                f"Счёт остановлен из-за: {stop_reason}\n"
                f"Достигнутая точность метода: {round(diff, 12)}\n"
                f"Тестовая задача решена с погрешностью: {round(eps, 12)}\n"
                f"Времени затрачено: {time} с\n"
                f"Количество итераций: {int(S)}\n"
                f"Время на одну итерацию: {time_dS} с"
            )
        else:
            res, time, S, eps  = map(lambda x: x[0], matrices[0])
            res_, S_, eps_, prec, time_ = map(lambda x: x[0], matrices[1])
            
            sum_time = round(time + time_, 6)
            time_dS = round(float(time) / float(S), 6)
            time_dS_ = round(float(time_) / float(S_), 6)
        
            stop_reason_main = "превышения количества итераций" if res == 0.0 else "выхода на точность"
            stop_reason_half = "превышения количества итераций" if res_ == 0.0 else "выхода на точность"
        
            text_to_display = (
                "Задача на основной сетке:\n"
                f"Счёт остановлен из-за: {stop_reason_main}\n"
                f"Достигнутая точность мотода: {round(eps, 12)}\n"
                f"Времени затрачено: {time} с\n"
                f"Количество итераций: {int(S)}\n"
                f"Время на одну итерацию: {time_dS} с\n"
                "Задача на контрольной сетке:\n"
                f"Счёт остановлен из-за: {stop_reason_half}\n"
                f"Максимальная разность двух приближений: {round(prec, 12)}\n"
                f"Точность итерационного метода: {round(eps_, 12)}\n"
                f"Времени затрачено: {time_} с\n"
                f"Количество итераций: {int(S_)}\n"
                f"Время на одну итерацию: {time_dS_} с\n"
                f"Всего времени затрачено: {sum_time} с"
            )

        self.second_window = SecondWindow(text_to_display)
        self.second_window.show()
        
    def set_header_table(self, table, X, Y, n, m):
        for col in range(n + 1):
            item_x = QTableWidgetItem(str(X[col]))
            table.setHorizontalHeaderItem(col, item_x)
        for row in range(m + 1):
            item_y = QTableWidgetItem(str(Y[row]))
            table.setVerticalHeaderItem(row, item_y)

    def plotting(self):
        n = int(self.le_n.text())
        m = int(self.le_m.text())
        S = int(self.le_S.text())
        eps = float(self.le_eps.text())
        S_ = int(self.le_S_.text())
        eps_ = float(self.le_eps_.text())
        
        lib_dir = os.path.join(os.curdir, "dll", "Release", "lab_lib.dll")
        lib = ctypes.windll.LoadLibrary(lib_dir)
        
        if (self.combo_task.currentText() == "Тестовая"):
            my_func = lib.solve_test
            my_func.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int,  ctypes.c_double]
            my_func.restype = ctypes.c_void_p
            my_func(n, m, S, eps)
        else:
            my_func = lib.solve_main
            my_func.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int,  ctypes.c_double, ctypes.c_int,  ctypes.c_double]
            my_func.restype = ctypes.c_void_p
            my_func(n, m, S, eps, S_, eps_)
        
        tables = [self.table_1, self.table_2, self.table_3]
        
        matrices = self.processing_results("result//numerical.txt")
        X = matrices[0][0].copy()
        Y = matrices[1][0].copy()
        if (self.combo_task.currentText() == "Тестовая"):
            for i in range(len(tables)):
                if (i == 0):
                    tables[i].setColumnCount(n + 1)
                    tables[i].setRowCount(m + 1)
                
                    self.set_header_table(tables[i], X, Y, n, m)
                    
                    for row in range(m + 1):
                        for col in range(n + 1):
                            item = QTableWidgetItem(str(matrices[2][row][col]))
                            tables[i].setItem(row, col, item)
                elif (i == 1):
                    auxiliary_matrix = self.processing_results("result//auxiliary.txt")
                    tables[i].setColumnCount(n + 1)
                    tables[i].setRowCount(m + 1)
                    
                    self.set_header_table(tables[i], X, Y, n, m)
                    
                    for row in range(m + 1):
                        for col in range(n + 1):
                            item = QTableWidgetItem(str(auxiliary_matrix[0][row][col]))
                            tables[i].setItem(row, col, item)
        else:
            X_ = matrices[2][0].copy()
            Y_ = matrices[3][0].copy()
            for i in range(len(tables)):
                if (i == 0):
                    tables[i].setColumnCount(n + 1)
                    tables[i].setRowCount(m + 1)
                    
                    self.set_header_table(tables[i], X, Y, n, m)
                    
                    for row in range(m + 1):
                        for col in range(n + 1):
                            item = QTableWidgetItem(str(matrices[4][row][col]))
                            tables[i].setItem(row, col, item)
                elif (i == 1):
                    tables[i].setColumnCount(2 * n + 1)
                    tables[i].setRowCount(2 * m + 1)
                    
                    self.set_header_table(tables[i], X_, Y_, 2 * n, 2 * m)
                    
                    for row in range(2 * m + 1):
                        for col in range(2 * n + 1):
                            item = QTableWidgetItem(str(matrices[5][row][col]))
                            tables[i].setItem(row, col, item)
        
        x = matrices[0][0]
        y = matrices[1][0]
        if (self.combo_task.currentText() == "Тестовая"):
            z_numeral = np.array(matrices[2])
        else:
            z_numeral = np.array(matrices[4])
        
        matrices = self.processing_results("result//diff_solutions.txt")
        
        tables[i].setColumnCount(n + 1)
        tables[i].setRowCount(m + 1)
        
        self.set_header_table(tables[i], X, Y, n, m)
        
        for row in range(m + 1):
            for col in range(n + 1):
                item = QTableWidgetItem(str(matrices[0][row][col]))
                tables[i].setItem(row, col, item)
        
        z_diff = np.array(matrices[0])

        X, Y = np.meshgrid(x, y)

        self.plt.plot_surface(X, Y, z_numeral, cmap = 'magma', alpha=0.75)
        graph = self.plt.plot_surface(X, Y, z_diff, cmap = 'ocean', alpha=0.95)
        self.plt.set_xlim(auto=True)
        self.plt.set_ylim(auto=True)
        self.plt.set_zlim(auto=True)
        
        self.plot.canvas.draw()
        
        self.get_informations("result//results.txt")
            

if __name__ == '__main__':
    app = QApplication(sys.argv)
    mainwindow = MainWindow()
    mainwindow.show()

    sys.exit(app.exec_())