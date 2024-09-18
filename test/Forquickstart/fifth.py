import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QAction, QFileDialog, QMessageBox, QWidget, QVBoxLayout, QHBoxLayout, QGroupBox, QRadioButton, QLabel, QLineEdit, QPushButton, QListWidget, QListWidgetItem, QSplitter
from PyQt5.QtCore import Qt, QThread, pyqtSignal, pyqtSlot
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
import os

# 데이터 로딩을 위한 Worker 클래스 (데이터 경로만 백그라운드에서 처리)
class DataLoadingWorker(QThread):
    file_path_ready = pyqtSignal(str)
    error_occurred = pyqtSignal(str)

    def __init__(self, file_name):
        super().__init__()
        self.file_name = file_name

    def run(self):
        try:
            # 파일 경로만 전달 (scanpy는 메인 스레드에서 처리)
            self.file_path_ready.emit(self.file_name)
        except Exception as e:
            self.error_occurred.emit(str(e))

# 시각화를 위한 Worker 클래스 (결과만 메인 스레드에 전달)
class VisualizationWorker(QThread):
    visualization_done = pyqtSignal()

    def __init__(self):
        super().__init__()

    def run(self):
        # 시각화 완료 신호만 전달 (실제 시각화는 메인 스레드에서 처리)
        self.visualization_done.emit()

class MainWindow(QMainWindow):

    def __init__(self):
        super().__init__()
        self.adata = None  # Scanpy AnnData object to hold the loaded data
        self.initUI()

    def initUI(self):
        # Setting up the Menu bar
        menubar = self.menuBar()
        file_menu = menubar.addMenu('&File')

        # Adding 'File Open' action
        open_action = QAction('File Open', self)
        open_action.setShortcut('Ctrl+O')
        open_action.setStatusTip('Open Anndata File')
        open_action.triggered.connect(self.open_file)
        file_menu.addAction(open_action)

        # Setting up the Status bar
        self.statusbar = self.statusBar()
        self.statusbar.showMessage('Ready')

        # Setting up the Main Window Layout
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)

        # Create the main horizontal splitter
        main_splitter = QSplitter(Qt.Horizontal)
        self.central_widget.setLayout(QVBoxLayout())
        self.central_widget.layout().addWidget(main_splitter)

        # Initialize the UI sections
        self.setup_left_area(main_splitter)
        self.setup_central_area(main_splitter)
        self.setup_right_area(main_splitter)

        # Set the split ratios (1:3:1)
        main_splitter.setSizes([200, 850, 150])

        # Main window settings
        self.setWindowTitle('PyQt5 Application')
        self.setGeometry(100, 100, 1200, 800)

    def setup_left_area(self, parent_splitter):
        # Create the left area
        left_panel = QGroupBox()
        left_layout = QVBoxLayout(left_panel)

        # Create Clusters group box
        clusters_group_box = QGroupBox('Clusters')
        clusters_layout = QVBoxLayout(clusters_group_box)
        left_layout.addWidget(clusters_group_box)

        # Create the list widget for displaying Clusters
        self.clusters_list_widget = QListWidget()
        clusters_layout.addWidget(self.clusters_list_widget)

        # Adding the Spatial Plot group box under Clusters
        spatial_plot_group_box = QGroupBox('Spatial Plot')
        spatial_plot_layout = QVBoxLayout(spatial_plot_group_box)
        left_layout.addWidget(spatial_plot_group_box)

        # First row: All/Select radio buttons
        radio_layout = QHBoxLayout()
        self.radio_all = QRadioButton('All')
        self.radio_select = QRadioButton('Select')
        radio_layout.addWidget(self.radio_all)
        radio_layout.addWidget(self.radio_select)
        self.radio_all.setChecked(True)  # Set 'All' as default selected option
        spatial_plot_layout.addLayout(radio_layout)

        # Second row: Axis label and input box
        axis_layout = QHBoxLayout()
        axis_label = QLabel('Axis :')
        self.axis_input = QLineEdit()
        self.axis_input.setText('0')  # Default value is 0
        axis_layout.addWidget(axis_label)
        axis_layout.addWidget(self.axis_input)
        spatial_plot_layout.addLayout(axis_layout)

        # Third row: Start button (aligned to the right)
        start_button = QPushButton('Start')
        start_button.clicked.connect(self.plot_spatial_scatter)
        spatial_plot_layout.addWidget(start_button)

        # Create Other Clusters group box
        other_clusters_group_box = QGroupBox('Other Clusters')
        other_clusters_group_box.setCheckable(True)
        other_clusters_group_box.setChecked(False)
        other_clusters_layout = QVBoxLayout(other_clusters_group_box)
        left_layout.addWidget(other_clusters_group_box)

        # Create the list widget for displaying Clusters
        self.other_clusters_list_widget = QListWidget()
        other_clusters_layout.addWidget(self.other_clusters_list_widget)

        # Adding the DEG Plot Button
        deg_plot_button = QPushButton("DEG Analysis")
        deg_plot_button.clicked.connect(self.deg_plot)
        left_layout.addWidget(deg_plot_button)

        # Adding the Save DEG results Button
        self.deg_save_button = QPushButton("Save DEG Result")
        self.deg_save_button.clicked.connect(self.save_deg_result)
        self.deg_save_button.setEnabled(False)
        left_layout.addWidget(self.deg_save_button)

        # Sync the 'Clusters' and 'Other Clusters' list widget
        self.clusters_list_widget.itemChanged.connect(self.sync_lists)
        self.other_clusters_list_widget.itemChanged.connect(self.sync_lists)

        # Add the left panel to the parent splitter
        parent_splitter.addWidget(left_panel)

    def setup_right_area(self, parent_splitter):
        # Create the right area
        right_panel = QGroupBox()
        right_layout = QVBoxLayout(right_panel)

        # Create Genes group box
        genes_group_box = QGroupBox('Genes')
        genes_layout = QVBoxLayout(genes_group_box)
        right_layout.addWidget(genes_group_box)

        # Create the list widget for displaying genes
        self.genes_list_widget = QListWidget()
        genes_layout.addWidget(self.genes_list_widget)

        # Create Plotting radio button inside a group box
        genes_plotting_group_box = QGroupBox('Plotting')
        genes_plotting_layout = QVBoxLayout(genes_plotting_group_box)
        self.dot_plot_radio_button = QRadioButton('Dot plot')
        self.violin_plot_radio_button = QRadioButton('Violin plot')
        self.feature_plot_radio_button = QRadioButton('Feature plot')
        self.scatter_plot_radio_button = QRadioButton('Scatter plot')
        start_button = QPushButton('Start')
        start_button.clicked.connect(self.plot_genes)
        genes_plotting_group_box.setLayout(genes_plotting_layout)
        genes_plotting_layout.addWidget(self.dot_plot_radio_button)
        genes_plotting_layout.addWidget(self.violin_plot_radio_button)
        genes_plotting_layout.addWidget(self.feature_plot_radio_button)
        genes_plotting_layout.addWidget(self.scatter_plot_radio_button)
        genes_plotting_layout.addWidget(start_button)
        right_layout.addWidget(genes_plotting_group_box)

        # Add the right panel to the parent splitter
        parent_splitter.addWidget(right_panel)

    def setup_central_area(self, parent_splitter):
        # Create the central area
        central_panel = QWidget()
        self.central_layout = QVBoxLayout(central_panel)
        self.plot_canvas = FigureCanvas(plt.Figure())
        self.central_layout.addWidget(self.plot_canvas)
        parent_splitter.addWidget(central_panel)

    def open_file(self):
        if self.adata is not None:
            reply = QMessageBox.question(self, 'Warning', "Loading a new file will overwrite the existing data. Do you want to continue?", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)

            if reply == QMessageBox.No:
                return  # Exit the function if the user chooses not to continue

        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getOpenFileName(self, "Open HDF5 File", "", "HDF5 Files (*.h5ad *.hdf5);;All Files (*)", options=options)

        if file_name:
            self.start_data_loading(file_name)

    def start_data_loading(self, file_name):
        self.statusbar.showMessage('Preparing to load data...')
        
        # Worker 스레드를 시작하여 파일 경로 준비 (실제 로딩은 메인 스레드에서)
        self.worker = DataLoadingWorker(file_name)
        self.worker.file_path_ready.connect(self.on_file_path_ready)
        self.worker.error_occurred.connect(self.on_data_load_error)
        self.worker.start()

    @pyqtSlot(str)
    def on_file_path_ready(self, file_name):
        import scanpy as sc
        self.statusbar.showMessage('Loading data in the main thread...')

        try:
            # 메인 스레드에서 데이터를 로드합니다.
            self.adata = sc.read_h5ad(file_name)
            self.update_ui_after_data_loaded()
        except Exception as e:
            self.on_data_load_error(str(e))

    @pyqtSlot()
    def update_ui_after_data_loaded(self):
        # UI 업데이트 (데이터 로딩 후 처리)
        self.update_lists()
        self.show_spatial_scatter_plot()

        base_name = os.path.basename(self.worker.file_name)
        self.statusbar.showMessage(f"Loaded '{base_name}'")
        print(f'Data loaded from: {base_name}')

        # 스레드 종료 후 스레드 객체 삭제
        self.worker.deleteLater()

    @pyqtSlot(str)
    def on_data_load_error(self, error_message):
        QMessageBox.critical(self, "Error", f"Failed to load file: {error_message}")
        self.statusbar.showMessage('Ready')
        print(f"Error loading file: {error_message}")

        # 스레드 종료 후 스레드 객체 삭제
        self.worker.deleteLater()

    def start_visualization(self, groups=None, legend=None):
        self.statusbar.showMessage('Preparing to visualize data...')

        # VisualizationWorker가 완료 신호만 보냄 (실제 시각화는 메인 스레드에서)
        self.vis_worker = VisualizationWorker()
        self.vis_worker.visualization_done.connect(lambda: self.on_visualization_done(groups, legend))
        self.vis_worker.start()

    @pyqtSlot()
    def on_visualization_done(self, groups, legend):
        self.statusbar.showMessage('Visualizing data in the main thread...')
        self.perform_visualization(groups, legend)
        self.vis_worker.deleteLater()

    def perform_visualization(self, groups=None, legend=None):
        import squidpy as sq
        from matplotlib.transforms import Affine2D
        import numpy as np

        axis = 360
        num = float(self.axis_input.text())
        if num > 0:
            axis = 180 / num

        if self.adata is not None:
            # Clear any existing plot
            self.plot_canvas.figure.clear()

            # Create a new plot and add it to the central layout
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
            sq.pl.spatial_scatter(self.adata, shape=None, color="cell_type_2", groups=groups, ax=ax, legend_loc=legend, frameon=False, title=None)

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
            self.statusbar.showMessage('Visualization complete')

    def update_lists(self):
        # Update List
        if self.adata is not None:
            self.update_list(self.clusters_list_widget)
            self.update_list(self.other_clusters_list_widget)
            self.update_list(self.genes_list_widget, is_gene=True)

    def update_list(self, list_widget, is_gene=False):
        # Clear the list widget before adding new items
        list_widget.clear()

        # Extract and sort cluster or gene list
        items = sorted(self.adata.var_names.to_list()) if is_gene else sorted(self.adata.obs['cell_type_2'].unique())
        
        # Populate the list widget with sorted items and add checkboxes
        for item in items:
            list_item = QListWidgetItem(item)
            list_item.setFlags(list_item.flags() | Qt.ItemIsUserCheckable)
            list_item.setCheckState(Qt.Unchecked)
            list_widget.addItem(list_item)

    def sync_lists(self, item):
        self.clusters_list_widget.blockSignals(True)
        self.other_clusters_list_widget.blockSignals(True)

        if item.listWidget() == self.clusters_list_widget:
            other_list = self.other_clusters_list_widget
        else:
            other_list = self.clusters_list_widget

        for i in range(other_list.count()):
            other_item = other_list.item(i)
            if other_item.text() == item.text():
                other_item.setFlags(other_item.flags() & ~Qt.ItemIsEnabled if item.checkState() == Qt.Checked else other_item.flags() | Qt.ItemIsEnabled)

        self.clusters_list_widget.blockSignals(False)
        self.other_clusters_list_widget.blockSignals(False)

    def show_spatial_scatter_plot(self, groups=None, legend=None):
        # 시각화 작업을 메인 스레드에서 수행
        self.start_visualization(groups, legend)

    def plot_spatial_scatter(self):
        if self.radio_all.isChecked():
            self.show_spatial_scatter_plot(groups=None, legend=None)
        elif self.radio_select.isChecked():
            self.selected_clusters = [item.text() for item in self.clusters_list_widget.findItems("*", Qt.MatchWildcard) if item.checkState() == Qt.Checked]
            if not self.selected_clusters:
                QMessageBox.warning(self, "Warning", "Please select at least one cluster.")
                return
            else:
                self.show_spatial_scatter_plot(groups=self.selected_clusters, legend='right margin')

    def plot_genes(self):
        if self.adata is not None:
            selected_genes = [item.text() for item in self.genes_list_widget.findItems("*", Qt.MatchWildcard) if item.checkState() == Qt.Checked]
            
            if not selected_genes:
                QMessageBox.warning(self, "Warning", "Please select at least one gene.")
                return

            axis = 360
            num = float(self.axis_input.text())
            if num > 0:
                axis = 180 / num

            if self.scatter_plot_radio_button.isChecked():
                import squidpy as sq
                from matplotlib.transforms import Affine2D
                import numpy as np
                                    
                # Determine the number of subplots needed
                num_genes = len(selected_genes)
                fig, axes = plt.subplots(1, num_genes, figsize=(5 * num_genes, 5))

                # If only one gene is selected, wrap the axes in a list for consistency
                if num_genes == 1:
                    axes = [axes]

                for gene, ax in zip(selected_genes, axes):
                    # Extract x and y coordinates from adata
                    x_coords = self.adata.obsm['spatial'][:, 0]
                    y_coords = self.adata.obsm['spatial'][:, 1]

                    # Calculate the center of the plot
                    center_x = (x_coords.max() + x_coords.min()) / 2
                    center_y = (y_coords.max() + y_coords.min()) / 2

                    # Apply rotation transformation around the center
                    rotation = Affine2D().rotate_around(center_x, center_y, -np.pi / axis)
                    ax.transData = rotation + ax.transData

                    # Plot using squidpy
                    sq.pl.spatial_scatter(self.adata, shape=None, color=gene, ax=ax, frameon=False, title=None, cmap='Purples', size=1)

                    # Set the title
                    ax.set_title(f"{gene}")

                    # Calculate new limits after rotation
                    coords = rotation.transform(np.vstack([x_coords, y_coords]).T)
                    x_min, y_min = coords.min(axis=0)
                    x_max, y_max = coords.max(axis=0)

                    # Set new limits
                    ax.set_xlim(x_min-1000, x_max+1000)
                    ax.set_ylim(y_min-1000, y_max+1000)

            else:
                import scanpy as sc
                if self.dot_plot_radio_button.isChecked():
                    sc.pl.dotplot(self.adata, selected_genes, groupby="cell_type_2", categories_order=sorted(self.adata.obs['cell_type_2'].unique()))
                elif self.violin_plot_radio_button.isChecked():
                    sc.pl.stacked_violin(self.adata, selected_genes, groupby="cell_type_2", categories_order=sorted(self.adata.obs['cell_type_2'].unique()))
                elif self.feature_plot_radio_button.isChecked():
                    selected_genes.append('cell_type_2')
                    sc.pl.umap(self.adata, color=selected_genes, ncols=4)

            plt.show()

    def deg_plot(self):
        self.deg_save_button.setEnabled(True)

        self.selected_clusters = [item.text() for item in self.clusters_list_widget.findItems("*", Qt.MatchWildcard) if item.checkState() == Qt.Checked]
        self.selected_other_clusters = [item.text() for item in self.other_clusters_list_widget.findItems("*", Qt.MatchWildcard) if item.checkState() == Qt.Checked]

        if not self.selected_clusters or not self.selected_other_clusters:
            QMessageBox.warning(self, "Warning", "Please select at least one cluster in both Cluster 1 and Cluster 2.")
            return
        
        else:
            import scanpy as sc
            merged_cluster_A = '_'.join(self.selected_clusters)
            merged_cluster_B = '_'.join(self.selected_other_clusters)
                        
            self.adata.obs['cell_type_3'] = self.adata.obs['cell_type_2']
            self.adata.obs['cell_type_3'] = self.adata.obs['cell_type_3'].astype(str)
            self.adata.obs['cell_type_3'][self.adata.obs['cell_type_3'].isin(self.selected_clusters)] = merged_cluster_A
            self.adata.obs['cell_type_3'][self.adata.obs['cell_type_3'].isin(self.selected_other_clusters)] = merged_cluster_B                        
            self.adata.obs['cell_type_3'] = self.adata.obs['cell_type_3'].astype('category')
            
            sc.tl.rank_genes_groups(self.adata, 'cell_type_3', groups=[merged_cluster_A], reference=merged_cluster_B, method='wilcoxon')            
            sc.pl.rank_genes_groups(self.adata, n_genes=25, sharey=False)

    def save_deg_result(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        filePath, _ = QFileDialog.getSaveFileName(self, "Save CSV File", "", "CSV Files (*.csv);;All Files (*)", options=options)
        if 'rank_genes_groups' in self.adata.uns:
            import pandas as pd            
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
    from multiprocessing import freeze_support

    freeze_support()

    app = QApplication(sys.argv)
    main_window = MainWindow()
    main_window.show()
    sys.exit(app.exec_())
