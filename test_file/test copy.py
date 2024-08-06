import sys
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QMainWindow, QAction, QFileDialog, QMessageBox, QVBoxLayout, QWidget, QLabel, QGroupBox, QCheckBox, QListWidget, QListWidgetItem, QHBoxLayout, QRadioButton, QPushButton, QGridLayout, QDialog
from PyQt5.QtGui import QIcon, QFont
import scanpy as sc
import squidpy as sq
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt

class PlotDialog(QDialog):
    def __init__(self, adata, plot_type, genes=None, clusters1=None, clusters2=None, parent=None):
        super().__init__(parent)
        self.adata = adata
        self.plot_type = plot_type
        self.genes = genes
        self.clusters1 = clusters1
        self.clusters2 = clusters2
        self.initUI()

    def initUI(self):
        self.setWindowTitle('Plot')
        self.setGeometry(100, 100, 800, 600)
        self.layout = QVBoxLayout(self)

        self.plot_canvas = FigureCanvas(plt.Figure())
        self.layout.addWidget(self.plot_canvas)

        self.plot()

    def plot(self):
        self.plot_canvas.figure.clear()
        if self.plot_type == 'dotplot':
            sc.pl.dotplot(self.adata, self.genes, groupby="cell_type_2", categories_order=sorted(self.adata.obs['cell_type_2'].unique()))
        elif self.plot_type == 'violin':
            sc.pl.stacked_violin(self.adata, self.genes, groupby="cell_type_2", categories_order=sorted(self.adata.obs['cell_type_2'].unique()))
        elif self.plot_type == 'feature':
            sc.pl.umap(self.adata, color=self.genes, ncols=4)
        elif self.plot_type == 'spatial':
            ax = self.plot_canvas.figure.add_subplot(111)
            sq.pl.spatial_scatter(self.adata[self.adata.obs['cell_type_2'].isin(self.clusters1)], shape=None, color='cell_type_2', legend_loc=None, ax=ax, frameon=False)
        elif self.plot_type == 'deg':
            sc.tl.rank_genes_groups(self.adata, groupby='cell_type_2', groups=self.clusters1, reference=self.clusters2)
            sc.pl.rank_genes_groups(self.adata, groups=self.clusters1, n_genes=20, sharey=False)
        self.plot_canvas.draw()

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
        self.cluster_all_group_box = QGroupBox()
        self.cluster_all_group_box.setFixedWidth(250)
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

        self.spatial_plot_button = QPushButton("Spatial plot")
        self.spatial_plot_button.clicked.connect(self.spatial_plot)

        self.deg_analysis_button = QPushButton("DEG Analysis")
        self.deg_analysis_button.clicked.connect(self.deg_analysis)

        self.cluster_all_layout.addWidget(self.cluster1_group_box)
        self.cluster_all_layout.addWidget(self.cluster2_group_box)
        self.cluster_all_layout.addWidget(self.spatial_plot_button)
        self.cluster_all_layout.addWidget(self.deg_analysis_button)
        self.cluster_all_group_box.setLayout(self.cluster_all_layout)
        
        self.main_layout.addWidget(self.cluster_all_group_box)

        self.cluster1_list.itemChanged.connect(self.sync_lists)
        self.cluster2_list.itemChanged.connect(self.sync_lists)

    def create_gene_group_box(self):
        self.gene_group_box = QGroupBox("Gene")
        self.gene_group_box.setFixedWidth(200)
        self.gene_layout = QVBoxLayout()
        self.gene_list = QListWidget()
        self.gene_layout.addWidget(self.gene_list)
        self.gene_group_box.setLayout(self.gene_layout)

        self.analysis_group_box = QGroupBox("Analysis")
        self.analysis_layout = QVBoxLayout()
        self.dot_plot_radio = QRadioButton("Dot plot")
        self.violin_plot_radio = QRadioButton("Violin plot")
        self.feature_plot_radio = QRadioButton("Feature plot")
        self.start_plotting_button = QPushButton("Start Plotting")
        self.start_plotting_button.clicked.connect(self.start_plotting)
        self.analysis_layout.addWidget(self.dot_plot_radio)
        self.analysis_layout.addWidget(self.violin_plot_radio)
        self.analysis_layout.addWidget(self.feature_plot_radio)
        self.analysis_layout.addWidget(self.start_plotting_button)
        self.analysis_group_box.setLayout(self.analysis_layout)

        self.grid_layout = QGridLayout()
        self.grid_layout.addWidget(self.gene_group_box, 0, 0)
        self.grid_layout.addWidget(self.analysis_group_box, 1, 0)
        self.main_layout.addLayout(self.grid_layout)

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
            sq.pl.spatial_scatter(self.adata, shape=None, color=['cell_type_2'], legend_loc=None, ax=ax, frameon=False)
            self.plot_canvas.draw()

    def start_plotting(self):
        if self.adata is not None:
            selected_genes = [item.text() for item in self.gene_list.findItems("*", Qt.MatchWildcard) if item.checkState() == Qt.Checked]
            
            if not selected_genes:
                QMessageBox.warning(self, "Warning", "Please select at least one gene.")
                return

            plot_type = None
            if self.dot_plot_radio.isChecked():
                plot_type = 'dotplot'
            elif self.violin_plot_radio.isChecked():
                plot_type = 'violin'
            elif self.feature_plot_radio.isChecked():
                plot_type = 'feature'

            if plot_type is not None:
                self.plot(selected_genes, plot_type)

    def plot(self, genes, plot_type):
        if plot_type == 'dotplot':
            sc.pl.dotplot(self.adata, genes, groupby="cell_type_2", categories_order=sorted(self.adata.obs['cell_type_2'].unique()))
        elif plot_type == 'violin':
            sc.pl.stacked_violin(self.adata, genes, groupby="cell_type_2", categories_order=sorted(self.adata.obs['cell_type_2'].unique()))
        elif plot_type == 'feature':
            sc.pl.umap(self.adata, color=genes, ncols=4)
        plt.show()

    def spatial_plot(self):
        selected_clusters1 = [item.text() for item in self.cluster1_list.findItems("*", Qt.MatchWildcard) if item.checkState() == Qt.Checked]
        if not selected_clusters1:
            QMessageBox.warning(self, "Warning", "Please select at least one cluster in Cluster 1.")
            return
        self.plot_dialog = PlotDialog(self.adata, 'spatial', clusters1=selected_clusters1)
        self.plot_dialog.exec_()

    def deg_analysis(self):
        selected_clusters1 = [item.text() for item in self.cluster1_list.findItems("*", Qt.MatchWildcard) if item.checkState() == Qt.Checked]
        selected_clusters2 = [item.text() for item in self.cluster2_list.findItems("*", Qt.MatchWildcard) if item.checkState() == Qt.Checked]
        if not selected_clusters1 or not selected_clusters2:
            QMessageBox.warning(self, "Warning", "Please select at least one cluster in both Cluster 1 and Cluster 2.")
            return
        self.plot_dialog = PlotDialog(self.adata, 'deg', clusters1=selected_clusters1, clusters2=selected_clusters2)
        self.plot_dialog.exec_()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = MyApp()
    sys.exit(app.exec_())
