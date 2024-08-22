from PyQt5.QtWidgets import QMainWindow, QAction, QVBoxLayout, QWidget, QHBoxLayout, QListWidget, QPushButton, QGroupBox, QRadioButton, QGridLayout, QListWidgetItem, QLabel, QLineEdit, QMessageBox
from PyQt5.QtGui import QIcon, QFontMetrics
from PyQt5.QtCore import Qt
import scanpy as sc
import squidpy as sq
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
from matplotlib.transforms import Affine2D
import numpy as np
from logic import handle_load_data, handle_plot_umap, handle_sync_lists, handle_deg_plot, handle_spatial_plot, handle_start_plotting, handle_save_file, handle_plotting_clusters
from style import cluster_all_styles, cluster1_styles, cluster2_styles, gene_group_box_styles, analysis_group_box_styles

class MyApp(QMainWindow):

    def __init__(self, config):
        super().__init__()
        self.adata = None
        self.config = config
        self.initUI()

    def initUI(self):
        openFile = QAction(QIcon('openfile.jpg'), 'Open File', self)
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

        self.setWindowTitle('Vizgen Analysis Tool')
        self.setGeometry(100, 100, self.config.get('window_width', 1200), self.config.get('window_height', 800))
        self.show()

    def create_cluster_all_group_box(self):
        self.cluster_all_group_box = QGroupBox()
        self.cluster_all_group_box.setFixedWidth(250)
        self.cluster_all_group_box.setStyleSheet(cluster_all_styles)
        self.cluster_all_layout = QVBoxLayout()

        self.cluster1_group_box = QGroupBox("Clusters")
        self.cluster1_group_box.setStyleSheet(cluster1_styles)
        self.cluster1_layout = QVBoxLayout()
        self.cluster1_list = QListWidget()
        self.cluster1_list.setStyleSheet(cluster1_styles)
        self.cluster1_layout.addWidget(self.cluster1_list)
        self.cluster1_group_box.setLayout(self.cluster1_layout)

        # Adding the Spatial Plot group box under Clusters
        self.spatial_plot_group_box = QGroupBox('Spatial Plot')
        self.spatial_plot_layout = QVBoxLayout()
        self.spatial_plot_group_box.setLayout(self.spatial_plot_layout)

        # First row: All/Select radio buttons
        radio_layout = QHBoxLayout()
        self.radio_all = QRadioButton('All')
        self.radio_select = QRadioButton('Select')
        radio_layout.addWidget(self.radio_all)
        radio_layout.addWidget(self.radio_select)
        self.radio_all.setChecked(True)  # Set 'All' as default selected option
        self.spatial_plot_layout.addLayout(radio_layout)

        # Second row: Axis label and input box
        axis_layout = QHBoxLayout()
        axis_label = QLabel('Axis :')
        self.axis_input = QLineEdit()
        self.axis_input.setText('0')  # Default value is 0
        axis_layout.addWidget(axis_label)
        axis_layout.addWidget(self.axis_input)
        self.spatial_plot_layout.addLayout(axis_layout)

        # Third row: Start button (aligned to the right)
        start_button = QPushButton('Start')
        self.spatial_plot_layout.addWidget(start_button)
        start_button.clicked.connect(self.plotting_clusters)        

        self.cluster2_group_box = QGroupBox("Other Clusters")
        self.cluster2_group_box.setStyleSheet(cluster2_styles)
        self.cluster2_group_box.setCheckable(True)
        self.cluster2_group_box.setChecked(False)
        self.cluster2_layout = QVBoxLayout()
        self.cluster2_list = QListWidget()
        self.cluster2_list.setStyleSheet(cluster2_styles)
        self.cluster2_layout.addWidget(self.cluster2_list)
        self.cluster2_group_box.setLayout(self.cluster2_layout)

        #self.all_spatial_button = QPushButton("All Spatial plot")
        #self.spatial_plot_button = QPushButton("Spatial plot")
        self.deg_analysis_button = QPushButton("DEG Analysis")
        self.save_deg_button = QPushButton("Save DEG Result")
        #self.all_spatial_button.setStyleSheet(cluster_all_styles)
        #self.spatial_plot_button.setStyleSheet(cluster_all_styles)
        self.deg_analysis_button.setStyleSheet(cluster_all_styles)
        self.save_deg_button.setStyleSheet(cluster_all_styles)
        #self.all_spatial_button.clicked.connect(self.all_spatial)
        #self.spatial_plot_button.clicked.connect(self.spatial_plot)
        self.deg_analysis_button.clicked.connect(self.deg_plot)
        self.save_deg_button.clicked.connect(self.save_file)
        self.save_deg_button.setEnabled(False)

        #self.cluster_all_layout.addWidget(self.all_spatial_button)
        self.cluster_all_layout.addWidget(self.cluster1_group_box)
        #self.cluster_all_layout.addWidget(self.spatial_plot_button)
        self.cluster_all_layout.addWidget(self.spatial_plot_group_box)
        self.cluster_all_layout.addWidget(self.cluster2_group_box)
        self.cluster_all_layout.addWidget(self.deg_analysis_button)
        self.cluster_all_layout.addWidget(self.save_deg_button)
        self.cluster_all_group_box.setLayout(self.cluster_all_layout)
        
        self.main_layout.addWidget(self.cluster_all_group_box)

        self.cluster1_list.itemChanged.connect(self.sync_lists)
        self.cluster2_list.itemChanged.connect(self.sync_lists)

    def create_gene_group_box(self):
        self.gene_group_box = QGroupBox("Gene")
        self.gene_group_box.setStyleSheet(gene_group_box_styles)
        self.gene_group_box.setFixedWidth(200)
        self.gene_layout = QVBoxLayout()
        self.gene_list = QListWidget()
        self.gene_list.setStyleSheet(gene_group_box_styles)
        self.gene_layout.addWidget(self.gene_list)
        self.gene_group_box.setLayout(self.gene_layout)

        self.analysis_group_box = QGroupBox("Plotting")
        self.analysis_group_box.setStyleSheet(analysis_group_box_styles)
        self.analysis_layout = QVBoxLayout()
        self.dot_plot_radio = QRadioButton("Dot plot")
        self.violin_plot_radio = QRadioButton("Violin plot")
        self.feature_plot_radio = QRadioButton("Feature plot")
        self.spatial_plot_radio = QRadioButton("Spatial plot")
        self.dot_plot_radio.setStyleSheet(analysis_group_box_styles)
        self.violin_plot_radio.setStyleSheet(analysis_group_box_styles)
        self.feature_plot_radio.setStyleSheet(analysis_group_box_styles)
        self.spatial_plot_radio.setStyleSheet(analysis_group_box_styles)
        self.start_plotting_button = QPushButton("Start Plotting")
        self.start_plotting_button.setStyleSheet(analysis_group_box_styles)
        self.start_plotting_button.clicked.connect(self.start_plotting)
        self.analysis_layout.addWidget(self.dot_plot_radio)
        self.analysis_layout.addWidget(self.violin_plot_radio)
        self.analysis_layout.addWidget(self.feature_plot_radio)
        self.analysis_layout.addWidget(self.spatial_plot_radio)
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
            self.update_cluster_group_widths()
            self.update_gene_group_widths()

    def update_list(self, list_widget, is_gene=False):
        list_widget.clear()
        items = self.adata.var_names if is_gene else sorted(self.adata.obs['cell_type_2'].unique())
        for item in items:
            list_item = QListWidgetItem(item)
            list_item.setFlags(list_item.flags() | Qt.ItemIsUserCheckable)
            list_item.setCheckState(Qt.Unchecked)
            list_widget.addItem(list_item)

    def update_cluster_group_widths(self):
        self.update_cluster_group_width(self.cluster_all_group_box, self.cluster1_list)

    def update_cluster_group_width(self, group_box, list_widget):
        fm = QFontMetrics(list_widget.font())
        max_width = max([fm.horizontalAdvance(item.text()) for item in [list_widget.item(i) for i in range(list_widget.count())]])
        group_box.setFixedWidth(max_width + 100)

    def update_gene_group_widths(self):
        self.update_gene_group_width(self.gene_group_box, self.gene_list)

    def update_gene_group_width(self, group_box, list_widget):
        fm = QFontMetrics(list_widget.font())
        max_width = max([fm.horizontalAdvance(item.text()) for item in [list_widget.item(i) for i in range(list_widget.count())]])
        group_box.setFixedWidth(max_width + 75)

    def load_data(self):
        handle_load_data(self)

    def sync_lists(self, item):
        handle_sync_lists(self, item)

    #def all_spatial(self):
    #    handle_plot_umap(self)

    def plot_umap(self):
        handle_plot_umap(self)

    def start_plotting(self):
        handle_start_plotting(self)

    def plotting_clusters(self):
        handle_plotting_clusters(self)

    def plot(self, genes, plot_type):
        if plot_type == 'dotplot':
            sc.pl.dotplot(self.adata, genes, groupby="cell_type_2", categories_order=sorted(self.adata.obs['cell_type_2'].unique()))
        elif plot_type == 'violin':
            sc.pl.stacked_violin(self.adata, genes, groupby="cell_type_2", categories_order=sorted(self.adata.obs['cell_type_2'].unique()))
        elif plot_type == 'feature':
            sc.pl.umap(self.adata, color=genes, ncols=4)
        elif plot_type == 'spatial':
            sq.pl.spatial_scatter(self.adata, shape=None, color=genes, frameon=False, title=None, cmap='Purples', size=1, ncols=4)
        plt.show()

    def plot_clusters(self, clusters, plot_type):
        num = float(self.axis_input.text())
        if num == 0:
            axis = 360
        else:
            axis = 180/num        

        if plot_type == 'all':
            ax = self.plot_canvas.figure.add_subplot(111)
            # Extract x and y coordinates from adata
            x_coords = self.adata.obsm['spatial'][:, 0]
            y_coords = self.adata.obsm['spatial'][:, 1]

            # Calculate the center of the plot
            center_x = (x_coords.max() + x_coords.min()) / 2
            center_y = (y_coords.max() + y_coords.min()) / 2

            # Apply rotation transformation around the center
            rotation = Affine2D().rotate_around(center_x, center_y, -np.pi/(axis))
            ax.transData = rotation + ax.transData

            # Plot using squidpy
            sq.pl.spatial_scatter(self.adata, shape=None, color="cell_type_2", ax=ax, legend_loc=None,frameon=False, title=None)

            # Set the title
            ax.set_title('')

            # Calculate new limits after rotation
            coords = rotation.transform(np.vstack([x_coords, y_coords]).T)
            x_min, y_min = coords.min(axis=0)
            x_max, y_max = coords.max(axis=0)

            # Set new limits
            ax.set_xlim(x_min, x_max)
            ax.set_ylim(y_min, y_max)
            
        elif plot_type == 'select':
            if not clusters:
                QMessageBox.warning(self, "Warning", "Please select at least one cluster in Cluster 1.")
                return
            else:
                ax = self.plot_canvas.figure.add_subplot(111)
                # Extract x and y coordinates from adata
                x_coords = self.adata.obsm['spatial'][:, 0]
                y_coords = self.adata.obsm['spatial'][:, 1]

                # Calculate the center of the plot
                center_x = (x_coords.max() + x_coords.min()) / 2
                center_y = (y_coords.max() + y_coords.min()) / 2

                # Apply rotation transformation around the center
                rotation = Affine2D().rotate_around(center_x, center_y, -np.pi/(axis))
                ax.transData = rotation + ax.transData

                # Plot using squidpy
                sq.pl.spatial_scatter(self.adata, shape=None, color="cell_type_2", groups=clusters, ax=ax, frameon=False, title=None, legend_fontsize=15)

                # Set the title
                ax.set_title('')

                # Calculate new limits after rotation
                coords = rotation.transform(np.vstack([x_coords, y_coords]).T)
                x_min, y_min = coords.min(axis=0)
                x_max, y_max = coords.max(axis=0)

                # Set new limits
                ax.set_xlim(x_min, x_max)
                ax.set_ylim(y_min, y_max)
        
        self.plot_canvas.draw()
        

    #def spatial_plot(self):
    #    handle_spatial_plot(self)

    def deg_plot(self):
        handle_deg_plot(self)

    def save_file(self):
        handle_save_file(self)
