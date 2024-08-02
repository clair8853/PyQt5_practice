import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QAction, QVBoxLayout, QWidget
from data_processing import process_data
from ui_elements import open_file_dialog, plot_umap

class MyApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.statusBar()

        openFile = QAction('Open File', self)
        openFile.setShortcut('Ctrl+O')
        openFile.setStatusTip('Open New Anndata File')
        openFile.triggered.connect(self.load_data)

        menubar = self.menuBar()
        menubar.setNativeMenuBar(True)
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(openFile)

        self.setWindowTitle('My First PyQt5 App')
        self.setGeometry(100, 100, 800, 600)  # 창 크기를 좀 더 크게 조정
        self.show()
    
    def load_data(self):
        file_path = open_file_dialog()
        if file_path:
            adata = process_data(file_path)
            canvas = plot_umap(adata)
            
            central_widget = QWidget()
            layout = QVBoxLayout(central_widget)
            layout.addWidget(canvas)
            self.setCentralWidget(central_widget)

app = QApplication(sys.argv)
main_window = MyApp()
main_window.show()
sys.exit(app.exec_())
