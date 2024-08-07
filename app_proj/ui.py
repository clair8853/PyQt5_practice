from PyQt5.QtWidgets import QMainWindow, QAction, QVBoxLayout, QWidget, QHBoxLayout, QListWidget, QPushButton, QGroupBox, QRadioButton, QGridLayout, QListWidgetItem
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import Qt
import scanpy as sc
import squidpy as sq
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
from logic import handle_load_data, handle_plot_umap, handle_sync_lists, handle_deg_plot, handle_spatial_plot, handle_start_plotting, handle_save_file

class MyApp(QMainWindow):

    def __init__(self):
        super().__init__()
        self.adata = None
        self.initUI()

    def initUI(self):
        openFile = QAction(QIcon('./app_proj/openfile.jpg'), 'Open File', self)
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

        self.cluster1_group_box = QGroupBox("Clusters")
        self.cluster1_group_box.setStyleSheet("QGroupBox { font-size: 15pt; }")
        self.cluster1_layout = QVBoxLayout()
        self.cluster1_list = QListWidget()
        self.cluster1_list.setStyleSheet("QListWidget { font-size: 12pt; }")
        self.cluster1_layout.addWidget(self.cluster1_list)
        self.cluster1_group_box.setLayout(self.cluster1_layout)

        self.cluster2_group_box = QGroupBox("Other Clusters")
        self.cluster2_group_box.setStyleSheet("QGroupBox { font-size: 15pt; }")
        self.cluster2_group_box.setCheckable(True)
        self.cluster2_group_box.setChecked(False)
        self.cluster2_layout = QVBoxLayout()
        self.cluster2_list = QListWidget()
        self.cluster2_list.setStyleSheet("QListWidget { font-size: 12pt; }")
        self.cluster2_layout.addWidget(self.cluster2_list)
        self.cluster2_group_box.setLayout(self.cluster2_layout)

        self.all_spatial_button = QPushButton("All Spatial plot")
        self.all_spatial_button.setStyleSheet("QPushButton { font-size: 15pt; }")
        self.all_spatial_button.clicked.connect(self.all_spatial)
        self.spatial_plot_button = QPushButton("Spatial plot")
        self.spatial_plot_button.setStyleSheet("QPushButton { font-size: 13pt; }")
        self.spatial_plot_button.clicked.connect(self.spatial_plot)
        self.deg_analysis_button = QPushButton("DEG Analysis")
        self.deg_analysis_button.setStyleSheet("QPushButton { font-size: 13pt; }")
        self.deg_analysis_button.clicked.connect(self.deg_plot)
        self.save_deg_button = QPushButton("Save DEG Result")
        self.save_deg_button.setStyleSheet("QPushButton { font-size: 13pt; }")
        self.save_deg_button.clicked.connect(self.save_file)
        self.save_deg_button.setEnabled(False)

        self.cluster_all_layout.addWidget(self.all_spatial_button)
        self.cluster_all_layout.addWidget(self.cluster1_group_box)
        self.cluster_all_layout.addWidget(self.spatial_plot_button)
        self.cluster_all_layout.addWidget(self.cluster2_group_box)
        self.cluster_all_layout.addWidget(self.deg_analysis_button)
        self.cluster_all_layout.addWidget(self.save_deg_button)
        self.cluster_all_group_box.setLayout(self.cluster_all_layout)
        
        self.main_layout.addWidget(self.cluster_all_group_box)

        self.cluster1_list.itemChanged.connect(self.sync_lists)
        self.cluster2_list.itemChanged.connect(self.sync_lists)

    def create_gene_group_box(self):
        self.gene_group_box = QGroupBox("Gene")
        self.gene_group_box.setStyleSheet("QGroupBox { font-size: 15pt; }")
        self.gene_group_box.setFixedWidth(200)
        self.gene_layout = QVBoxLayout()
        self.gene_list = QListWidget()
        self.gene_list.setStyleSheet("QListWidget { font-size: 12pt; }")
        self.gene_layout.addWidget(self.gene_list)
        self.gene_group_box.setLayout(self.gene_layout)

        self.analysis_group_box = QGroupBox("Plotting")
        self.analysis_group_box.setStyleSheet("QGroupBox { font-size: 15pt; }")
        self.analysis_layout = QVBoxLayout()
        self.dot_plot_radio = QRadioButton("Dot plot")
        self.violin_plot_radio = QRadioButton("Violin plot")
        self.feature_plot_radio = QRadioButton("Feature plot")
        self.dot_plot_radio.setStyleSheet("QRadioButton { font-size: 12pt; }")
        self.violin_plot_radio.setStyleSheet("QRadioButton { font-size: 12pt; }")
        self.feature_plot_radio.setStyleSheet("QRadioButton { font-size: 12pt; }")
        self.start_plotting_button = QPushButton("Start Plotting")
        self.start_plotting_button.setStyleSheet("QPushButton { font-size: 13pt; }")
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

    def load_data(self):
        handle_load_data(self)

    def sync_lists(self, item):
        handle_sync_lists(self, item)

    def all_spatial(self):
        handle_plot_umap(self)

    def plot_umap(self):
        handle_plot_umap(self)

    def start_plotting(self):
        handle_start_plotting(self)

    def plot(self, genes, plot_type):
        if plot_type == 'dotplot':
            sc.pl.dotplot(self.adata, genes, groupby="cell_type_2", categories_order=sorted(self.adata.obs['cell_type_2'].unique()))
        elif plot_type == 'violin':
            sc.pl.stacked_violin(self.adata, genes, groupby="cell_type_2", categories_order=sorted(self.adata.obs['cell_type_2'].unique()))
        elif plot_type == 'feature':
            sc.pl.umap(self.adata, color=genes, ncols=4)
        plt.show()

    def spatial_plot(self):
        handle_spatial_plot(self)

    def deg_plot(self):
        handle_deg_plot(self)

    def save_file(self):
        handle_save_file(self)
