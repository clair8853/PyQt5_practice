import sys
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QMainWindow, QAction, QFileDialog, QMessageBox, QVBoxLayout, QWidget, QLabel, QGroupBox, QCheckBox, QListWidget, QListWidgetItem, QHBoxLayout
from PyQt5.QtGui import QIcon
import scanpy as sc
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt

class MyApp(QMainWindow):

    def __init__(self):
        super().__init__()
        self.adata = None
        self.initUI()

    def initUI(self):
        openFile = QAction(QIcon('./app_proj/exit.png'), 'Open File', self)
        openFile.setShortcut('Ctrl+O')
        openFile.setStatusTip('Open Anndata File')
        openFile.triggered.connect(self.load_data)

        self.statusBar()

        menubar = self.menuBar()
        menubar.setNativeMenuBar(False)
        filemenu = menubar.addMenu('&File')
        filemenu.addAction(openFile)

        self.main_widget = QWidget()
        self.setCentralWidget(self.main_widget)
        
        self.main_layout = QHBoxLayout(self.main_widget)
        
        self.create_cluster_all_group_box()
        self.create_plot_area()
        self.create_gene_group_box()

        self.setWindowTitle('Analysis Tools')
        self.setGeometry(100, 100, 1200, 800)
        self.show()

    def create_cluster_all_group_box(self):        
        self.cluster_all_group_box = QGroupBox("Cluster All")
        self.cluster_all_layout = QVBoxLayout()

        self.cluster1_group_box = QGroupBox("Cluster 1")
        self.cluster1_layout = QVBoxLayout()
        self.cluster1_list = QListWidget()
        self.cluster1_layout.addWidget(self.cluster1_list)
        self.cluster1_group_box.setLayout(self.cluster1_layout)

        self.cluster2_group_box = QGroupBox("Cluster 2")
        self.cluster2_group_box.setCheckable(True)
        self.cluster2_group_box.setChecked(False)
        self.cluster2_layout = QVBoxLayout()
        self.cluster2_list = QListWidget()
        self.cluster2_layout.addWidget(self.cluster2_list)
        self.cluster2_group_box.setLayout(self.cluster2_layout)

        self.cluster_all_layout.addWidget(self.cluster1_group_box)
        self.cluster_all_layout.addWidget(self.cluster2_group_box)
        self.cluster_all_group_box.setLayout(self.cluster_all_layout)
        
        self.main_layout.addWidget(self.cluster_all_group_box)

        self.cluster1_list.itemChanged.connect(self.sync_lists)
        self.cluster2_list.itemChanged.connect(self.sync_lists)

    def create_gene_group_box(self):
        self.gene_group_box = QGroupBox("Gene")
        self.gene_layout = QVBoxLayout()
        self.gene_list = QListWidget()
        self.gene_layout.addWidget(self.gene_list)
        self.gene_group_box.setLayout(self.gene_layout)
        self.main_layout.addWidget(self.gene_group_box)

    def create_plot_area(self):
        self.plot_widget = QWidget()
        self.plot_layout = QVBoxLayout(self.plot_widget)
        self.plot_canvas = FigureCanvas(plt.Figure())
        self.plot_layout.addWidget(self.plot_canvas)
        self.main_layout.addWidget(self.plot_widget)

    def load_data(self):
        if self.adata is not None:
            reply = QMessageBox.question(self, 'Message',
                                         "A file is already loaded. Do you want to load a new file and discard the current data?",
                                         QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if reply == QMessageBox.No:
                return
        
        options = QFileDialog.Options()
        file_name, _ = QFileDialog.getOpenFileName(self, "Open Anndata File", "", "Anndata Files (*.hdf5);;All Files (*)", options=options)
        if file_name:
            try:
                self.adata = sc.read_h5ad(file_name)
                self.statusBar().showMessage(f'Loaded {file_name}')
                self.update_data_info()
                self.update_lists()
                self.plot_umap()
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Could not load file:\n{str(e)}")

    def update_data_info(self):
        if self.adata is not None:
            info = f"Data loaded: {self.adata.shape[0]} cells, {self.adata.shape[1]} genes"
            self.statusBar().showMessage(info)

    def update_lists(self):
        if self.adata is not None:
            self.update_list(self.cluster1_list)
            self.update_list(self.cluster2_list)
            self.update_list(self.gene_list, is_gene=True)

    def update_list(self, list_widget, is_gene=False):
        list_widget.clear()
        items = self.adata.var_names if is_gene else sorted(self.adata.obs['cell_type_2'].unique())
        for item in items:
            list_item = QListWidgetItem(item)
            list_item.setFlags(list_item.flags() | Qt.ItemIsUserCheckable)
            list_item.setCheckState(Qt.Unchecked)
            list_widget.addItem(list_item)

    def sync_lists(self, item):
        self.cluster1_list.blockSignals(True)
        self.cluster2_list.blockSignals(True)

        if item.listWidget() == self.cluster1_list:
            other_list = self.cluster2_list
        else:
            other_list = self.cluster1_list

        for i in range(other_list.count()):
            other_item = other_list.item(i)
            if other_item.text() == item.text():
                other_item.setFlags(other_item.flags() & ~Qt.ItemIsEnabled if item.checkState() == Qt.Checked else other_item.flags() | Qt.ItemIsEnabled)

        self.cluster1_list.blockSignals(False)
        self.cluster2_list.blockSignals(False)

    def plot_umap(self):
        if self.adata is not None:
            self.plot_canvas.figure.clear()
            ax = self.plot_canvas.figure.add_subplot(111)
            sc.pl.umap(self.adata, color=["cell_type_2"], title='UMAP - Cell Type 2', frameon=False, ax=ax, show=False)
            self.plot_canvas.draw()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = MyApp()
    sys.exit(app.exec_())
