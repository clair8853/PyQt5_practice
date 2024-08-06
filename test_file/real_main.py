import sys
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QMainWindow, QAction, QFileDialog, QMessageBox, QVBoxLayout, QWidget, QLabel, QGroupBox, QCheckBox, QListWidget, QListWidgetItem, QHBoxLayout, QRadioButton, QPushButton, QGridLayout, QDialog
from PyQt5.QtGui import QIcon
import scanpy as sc
import squidpy as sq
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd

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
        self.deg_analysis_button.clicked.connect(self.deg_plot)

        self.save_deg_button = QPushButton("Save DEG Result")
        self.save_deg_button.clicked.connect(self.save_file)
        self.save_deg_button.setEnabled(False)

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
            sq.pl.spatial_scatter(self.adata, shape=None, color=['cell_type_2'], legend_loc=None, ax=ax, frameon=False, title=None, size=5)
            ax.set_title('')
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
        else:            
            self.plot_canvas.figure.clear()
            ax = self.plot_canvas.figure.add_subplot(111)
            self.adata.obs['cell_type_spatial'] = self.adata.obs['cell_type_2']
            self.adata.obs['cell_type_spatial'] = self.adata.obs['cell_type_spatial'].astype(str)
            self.adata.obs.loc[~self.adata.obs['cell_type_spatial'].isin(selected_clusters1), 'cell_type_spatial'] = 'other cell type'
            self.adata.obs['cell_type_spatial'] = self.adata.obs['cell_type_spatial'].astype('category')
            categories = self.adata.obs['cell_type_2'].cat.categories
            category_colors = dict(zip(categories, self.adata.uns['cell_type_2_colors']))
            new_categories = self.adata.obs['cell_type_spatial'].cat.categories
            new_category_colors = {cat: category_colors[cat] if cat in category_colors else '#f8f9fa' for cat in new_categories}
            custom_cmap = mcolors.ListedColormap([new_category_colors[key] for key in new_category_colors])

            sq.pl.spatial_scatter(self.adata, shape=None, color="cell_type_spatial", palette=custom_cmap, ax=ax, frameon=False, title=None, size=5)
            
            ax.set_title('')
            self.plot_canvas.draw()

    def deg_plot(self):
        self.save_deg_button.setEnabled(True)
        selected_clusters1 = [item.text() for item in self.cluster1_list.findItems("*", Qt.MatchWildcard) if item.checkState() == Qt.Checked]
        selected_clusters2 = [item.text() for item in self.cluster2_list.findItems("*", Qt.MatchWildcard) if item.checkState() == Qt.Checked]
        if not selected_clusters1 or not selected_clusters2:
            QMessageBox.warning(self, "Warning", "Please select at least one cluster in both Cluster 1 and Cluster 2.")
            return
        else:
            merged_cluster_A = '_'.join(selected_clusters1)
            merged_cluster_B = '_'.join(selected_clusters2)
                       
            self.adata.obs['cell_type_3'] = self.adata.obs['cell_type_2']
            self.adata.obs['cell_type_3'] = self.adata.obs['cell_type_3'].astype(str)
            self.adata.obs['cell_type_3'][self.adata.obs['cell_type_3'].isin(selected_clusters1)] = merged_cluster_A
            self.adata.obs['cell_type_3'][self.adata.obs['cell_type_3'].isin(selected_clusters2)] = merged_cluster_B                        
            self.adata.obs['cell_type_3'] = self.adata.obs['cell_type_3'].astype('category')
            
            sc.tl.rank_genes_groups(self.adata, 'cell_type_3', groups=[merged_cluster_A], reference=merged_cluster_B, method='wilcoxon')            
            sc.pl.rank_genes_groups(self.adata, n_genes=25, sharey=False)

    def save_file(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        filePath, _ = QFileDialog.getSaveFileName(self, "Save CSV File", "", "CSV Files (*.csv);;All Files (*)", options=options)
        if 'rank_genes_groups' in self.adata.uns:            
            result = self.adata.uns['rank_genes_groups']
            groups = result['names'].dtype.names            
            result_df = pd.DataFrame(
                {group + '_' + key: result[key][group]
                for group in groups for key in ['names', 'scores', 'logfoldchanges', 'pvals', 'pvals_adj']}
            )
            if filePath:
                if not filePath.endswith('.csv'):
                    filePath += '.csv'
                result_df.to_csv(filePath, index=False)
                info = f"Data saved: {filePath}"
                self.statusBar().showMessage(info)                     
        else:
            QMessageBox.critical(self, "Error", f"No results to save. Please run the analysis first.")


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = MyApp()
    sys.exit(app.exec_())
