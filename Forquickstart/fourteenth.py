import sys
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtWidgets import (QApplication, QMainWindow, QAction, QFileDialog, 
                             QMessageBox, QWidget, QVBoxLayout, QGridLayout, QGroupBox, 
                             QRadioButton, QLabel, QLineEdit, QPushButton, QProgressBar,
                             QListWidget, QListWidgetItem, QSplitter, QDialog)
import tempfile
from PyQt5.QtCore import Qt, QUrl, QThread, pyqtSignal
from numpy import array, cos, sin, dot, vstack, pi
import gc
from scanpy import read_h5ad

class DataProcessingThread(QThread):
    progress_changed = pyqtSignal(int)
    processing_finished = pyqtSignal(object)

    def __init__(self, file_name, parent=None):
        super().__init__(parent)
        self.file_name = file_name
        self.adata = None

    def run(self):
        self.progress_changed.emit(10)  # 파일 로드 시작
        try:
            self.adata = read_h5ad(self.file_name)  # 파일 로드
            self.progress_changed.emit(50)  # 파일 로드 완료
        except Exception as e:
            self.progress_changed.emit(0)
            self.processing_finished.emit(e)  # 에러 발생 시 처리
            return

        self.progress_changed.emit(80)  # 유전자 리스트 업데이트 중

        # 추가적인 데이터 처리 로직 수행
        self.progress_changed.emit(100)  # 작업 완료
        self.processing_finished.emit(self.adata)  # 성공적으로 데이터 처리 완료


class MainWindow(QMainWindow):

    def __init__(self):
        super().__init__()
        self.adata = None  # Scanpy AnnData object to hold the loaded data
        self.progress_bar = None
        self.initUI()

    def initUI(self):
        menubar = self.menuBar()
        file_menu = menubar.addMenu('&File')
        open_action = QAction('File Open', self)
        open_action.setShortcut('Ctrl+O')
        open_action.triggered.connect(self.open_file)
        file_menu.addAction(open_action)

        tool_menu = menubar.addMenu('&Tool')
        self.deg_action = QAction('DEG Analysis', self)
        self.deg_action.setShortcut('Ctrl+A')
        self.deg_action.triggered.connect(self.open_deg_analysis)
        self.deg_action.setEnabled(False)
        tool_menu.addAction(self.deg_action)

        self.statusbar = self.statusBar()
        self.statusbar.showMessage('Ready')

        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)

        main_splitter = QSplitter(Qt.Horizontal)
        self.central_widget.setLayout(QVBoxLayout())
        self.central_widget.layout().addWidget(main_splitter)

        self.setup_left_area(main_splitter)
        self.setup_central_area(main_splitter)
        
        main_splitter.setSizes([200, 1000])

        self.progress_bar = QProgressBar(self)
        self.statusbar.addPermanentWidget(self.progress_bar)
        self.progress_bar.setValue(0)

        self.setWindowTitle('PyQt5 Application')
        self.setGeometry(100, 100, 1200, 800)

    def setup_left_area(self, parent_splitter):
        left_panel = QGroupBox()
        left_layout = QGridLayout(left_panel)

        spatial_plot_group_box = QGroupBox('Spatial Plot')
        spatial_plot_layout = QGridLayout(spatial_plot_group_box)
        left_layout.addWidget(spatial_plot_group_box, 0, 0, 1, 2)

        axis_label = QLabel('Axis :')
        self.axis_input = QLineEdit()
        self.axis_input.setText('0')  # Default value is 0
        spatial_plot_layout.addWidget(axis_label, 0, 0)
        spatial_plot_layout.addWidget(self.axis_input, 0, 1)

        start_button = QPushButton('Start')
        start_button.clicked.connect(self.plot_spatial_scatter)
        spatial_plot_layout.addWidget(start_button, 1, 0, 1, 2)

        genes_group_box = QGroupBox('Genes')
        genes_layout = QGridLayout(genes_group_box)
        left_layout.addWidget(genes_group_box, 1, 0, 1, 2)

        self.genes_list_widget = QListWidget()
        genes_layout.addWidget(self.genes_list_widget, 0, 0, 1, 2)

        genes_plotting_group_box = QGroupBox('Plotting')
        genes_plotting_layout = QGridLayout(genes_plotting_group_box)
        self.dot_plot_radio_button = QRadioButton('Dot plot')
        self.violin_plot_radio_button = QRadioButton('Violin plot')
        self.feature_plot_radio_button = QRadioButton('Feature plot')
        self.scatter_plot_radio_button = QRadioButton('Scatter plot')

        start_button = QPushButton('Start')
        start_button.clicked.connect(self.plot_genes)

        genes_plotting_layout.addWidget(self.dot_plot_radio_button, 0, 0)
        genes_plotting_layout.addWidget(self.violin_plot_radio_button, 0, 1)
        genes_plotting_layout.addWidget(self.feature_plot_radio_button, 1, 0)
        genes_plotting_layout.addWidget(self.scatter_plot_radio_button, 1, 1)
        genes_plotting_layout.addWidget(start_button, 2, 0, 1, 2)

        left_layout.addWidget(genes_plotting_group_box, 2, 0, 1, 2)
        parent_splitter.addWidget(left_panel)

    def setup_central_area(self, parent_splitter):
        central_panel = QWidget()
        self.central_layout = QVBoxLayout(central_panel)
        self.web_view = QWebEngineView(self)
        self.central_layout.addWidget(self.web_view)
        parent_splitter.addWidget(central_panel)
    
    def open_file(self):
        if self.adata is not None:
            reply = QMessageBox.question(self, 'Warning', "Loading a new file will overwrite the existing data. Do you want to continue?", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if reply == QMessageBox.No:
                return

        file_name, _ = QFileDialog.getOpenFileName(self, "Open HDF5 File", "", "HDF5 Files (*.h5ad *.hdf5);;All Files (*)")
        if file_name:
            self.start_data_processing(file_name)

    def start_data_processing(self, file_name):
        self.progress_bar.setValue(0)
        self.data_thread = DataProcessingThread(file_name)
        self.data_thread.progress_changed.connect(self.update_progress)
        self.data_thread.processing_finished.connect(self.on_data_loaded)
        self.data_thread.start()

    def update_progress(self, value):
        self.progress_bar.setValue(value)

    def on_data_loaded(self, result):
        if isinstance(result, Exception):
            QMessageBox.critical(self, "Error", f"Failed to load data: {str(result)}")
        else:
            self.adata = result
            self.update_ui_after_data_loaded()

    def update_ui_after_data_loaded(self):
        self.update_lists()
        self.plot_spatial_scatter()
        self.deg_action.setEnabled(True)
        self.progress_bar.setValue(100)
        self.statusbar.showMessage("Data loaded successfully.")

    def update_lists(self):
        if self.adata is not None:
            self.update_list(self.genes_list_widget, is_gene=True)

    def update_list(self, list_widget, is_gene=False):
        list_widget.clear()
        items = sorted(self.adata.var_names.to_list()) if is_gene else sorted(self.adata.obs['cell_type_2'].unique())
        for item in items:
            list_item = QListWidgetItem(item)
            list_item.setFlags(list_item.flags() | Qt.ItemIsUserCheckable)
            list_item.setCheckState(Qt.Unchecked)
            list_widget.addItem(list_item)

    def plot_spatial_scatter(self):
        import plotly.graph_objects as go
        axis = 360
        num = float(self.axis_input.text())
        if num > 0:
            axis = 180 / num
        
        if self.adata is not None:
            try:
                color_categories = self.adata.obs['cell_type_2'].astype('category').cat.categories
                color_map = {category: color for category, color in zip(color_categories, self.adata.uns['cell_type_2_colors'])}

                cell_type_categories_sorted = sorted(self.adata.obs['cell_type_2'].unique(), key=str)

                x_coords = self.adata.obsm['spatial'][:, 0]
                y_coords = self.adata.obsm['spatial'][:, 1]

                center_x = (x_coords.max() + x_coords.min()) / 2
                center_y = (y_coords.max() + y_coords.min()) / 2

                rotation_matrix = array([
                    [cos(pi / axis), -sin(pi / axis)],
                    [sin(pi / axis), cos(pi / axis)]
                ])
                rotated_coords = dot(vstack([x_coords - center_x, y_coords - center_y]).T, rotation_matrix)
                rotated_x_coords = rotated_coords[:, 0] + center_x
                rotated_y_coords = rotated_coords[:, 1] + center_y

                fig = go.Figure()

                for category in cell_type_categories_sorted:
                    mask = self.adata.obs['cell_type_2'] == category
                    fig.add_trace(go.Scattergl(
                        x=rotated_x_coords[mask],
                        y=rotated_y_coords[mask],
                        mode='markers',
                        name=str(category),
                        marker=dict(size=2, color=color_map[category])
                    ))

                fig.update_layout(
                    xaxis=dict(scaleanchor="y", scaleratio=1, range=[rotated_x_coords.min(), rotated_x_coords.max()], showticklabels=False),
                    yaxis=dict(scaleanchor="x", scaleratio=1, range=[rotated_y_coords.min(), rotated_y_coords.max()], showticklabels=False),
                    legend=dict(font=dict(size=15), itemsizing="constant", title=dict(font=dict(size=30))),
                    plot_bgcolor="white", paper_bgcolor="white"
                )

                with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as temp_file:
                    fig.write_html(temp_file.name)
                    self.web_view.setUrl(QUrl.fromLocalFile(temp_file.name))

            except Exception as e:
                QMessageBox.warning(self, 'Error', f'Could not plot data:\n{str(e)}')
            finally:
                gc.collect()

    def plot_genes(self):
        from matplotlib.pyplot import subplots, colorbar, show
        if self.adata is not None:
            selected_genes = [item.text() for item in self.genes_list_widget.findItems("*", Qt.MatchWildcard) if item.checkState() == Qt.Checked]
            if not selected_genes:
                QMessageBox.warning(self, "Warning", "Please select at least one gene.")
                return

            axis = 360
            num = float(self.axis_input.text())
            if num > 0:
                axis = 180 / num

            x_coords = self.adata.obsm['spatial'][:, 0]
            y_coords = self.adata.obsm['spatial'][:, 1]

            center_x = (x_coords.max() + x_coords.min()) / 2
            center_y = (y_coords.max() + y_coords.min()) / 2

            rotation_matrix = array([
                [cos(pi / axis), -sin(pi / axis)],
                [sin(pi / axis), cos(pi / axis)]
            ])
            rotated_coords = dot(vstack([x_coords - center_x, y_coords - center_y]).T, rotation_matrix)
            rotated_x_coords = rotated_coords[:, 0] + center_x
            rotated_y_coords = rotated_coords[:, 1] + center_y

            num_genes = len(selected_genes)
            rows = 1 if num_genes <= 4 else (num_genes + 3) // 4
            cols = num_genes if num_genes <= 4 else 4
            
            fig, axs = subplots(rows, cols, figsize=(cols*5, rows*5))
            axs = axs.flat if num_genes > 1 else [axs]

            for i, gene in enumerate(selected_genes):
                ax = axs[i]
                gene_expression = self.adata[:, gene].X.toarray().flatten()
                scatter = ax.scatter(x=rotated_x_coords, y=rotated_y_coords, c=gene_expression, cmap='Purples', s=1)
                ax.set_title(f'{gene}')
                ax.axis('off')
                ax.set_aspect('equal', 'box')
                colorbar(scatter, ax=ax)

            for j in range(i+1, len(axs)):
                fig.delaxes(axs[j])

            show()
            gc.collect()

    def open_deg_analysis(self):
        if self.adata is None:
            QMessageBox.warning(self, 'Error', 'No data loaded to perform DEG analysis.')
            return
        
        dialog = DEGAnalysisDialog(self.adata, self)
        dialog.show()

class DEGAnalysisDialog(QDialog):
    def __init__(self, adata, parent=None):
        super().__init__(parent)
        self.adata = adata
        self.initUI()

    def initUI(self):
        self.setWindowTitle('DEG Analysis')
        self.setGeometry(200, 200, 600, 400)
        layout = QVBoxLayout()

        clusters_layout = QGridLayout()
        clusters_group_box = QGroupBox('Clusters')
        clusters_vbox = QVBoxLayout(clusters_group_box)
        self.clusters_list_widget = QListWidget()
        self.populate_list(self.clusters_list_widget)
        clusters_vbox.addWidget(self.clusters_list_widget)
        clusters_layout.addWidget(clusters_group_box, 0, 0)

        other_clusters_group_box = QGroupBox('Other Clusters')
        other_clusters_vbox = QVBoxLayout(other_clusters_group_box)
        self.other_clusters_list_widget = QListWidget()
        self.populate_list(self.other_clusters_list_widget)
        other_clusters_vbox.addWidget(self.other_clusters_list_widget)
        clusters_layout.addWidget(other_clusters_group_box, 0, 1)

        self.clusters_list_widget.itemChanged.connect(self.sync_lists)
        self.other_clusters_list_widget.itemChanged.connect(self.sync_lists)

        layout.addLayout(clusters_layout)

        clear_button = QPushButton('Clear')
        clear_button.clicked.connect(lambda: self.clear_selection(self.clusters_list_widget))
        clear_button.clicked.connect(lambda: self.clear_selection(self.other_clusters_list_widget))
        layout.addWidget(clear_button)
        
        start_button = QPushButton('Start')
        start_button.clicked.connect(self.run_deg_analysis)
        layout.addWidget(start_button)

        save_button = QPushButton('Save')
        save_button.clicked.connect(self.save_deg_result)
        layout.addWidget(save_button)

        self.setLayout(layout)

    def populate_list(self, list_widget):
        clusters = sorted(self.adata.obs['cell_type_2'].unique())
        for cluster in clusters:
            list_item = QListWidgetItem(cluster)
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

    def clear_selection(self, list_widget):
        for i in range(list_widget.count()):
            item = list_widget.item(i)
            item.setCheckState(Qt.Unchecked)

    def run_deg_analysis(self):
        from scanpy import tools as tl
        from scanpy import plotting as pl

        selected_clusters = [item.text() for item in self.clusters_list_widget.findItems("*", Qt.MatchWildcard) if item.checkState() == Qt.Checked]
        selected_other_clusters = [item.text() for item in self.other_clusters_list_widget.findItems("*", Qt.MatchWildcard) if item.checkState() == Qt.Checked]

        if not selected_clusters or not selected_other_clusters:
            QMessageBox.warning(self, "Warning", "Please select at least one cluster in both Clusters and Other Clusters.")
            return

        merged_cluster_A = '_'.join(selected_clusters)
        merged_cluster_B = '_'.join(selected_other_clusters)
                    
        self.adata.obs['cell_type_3'] = self.adata.obs['cell_type_2'].astype(str)
        self.adata.obs['cell_type_3'][self.adata.obs['cell_type_3'].isin(selected_clusters)] = merged_cluster_A
        self.adata.obs['cell_type_3'][self.adata.obs['cell_type_3'].isin(selected_other_clusters)] = merged_cluster_B                        
        self.adata.obs['cell_type_3'] = self.adata.obs['cell_type_3'].astype('category')
        
        tl.rank_genes_groups(self.adata, 'cell_type_3', groups=[merged_cluster_A], reference=merged_cluster_B, method='wilcoxon')            
        pl.rank_genes_groups(self.adata, n_genes=25, sharey=False)
        QMessageBox.information(self, "Success", "DEG Analysis completed successfully.")

    def save_deg_result(self):
        from pandas import DataFrame
        filePath, _ = QFileDialog.getSaveFileName(self, "Save CSV File", "", "CSV Files (*.csv);;All Files (*)")
        if 'rank_genes_groups' in self.adata.uns:            
            result = self.adata.uns['rank_genes_groups']
            groups = result['names'].dtype.names            
            result_df = DataFrame(
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
            QMessageBox.critical(self, "Error", "No results to save. Please run the analysis first.")

if __name__ == '__main__':
    from multiprocessing import freeze_support
    freeze_support()

    app = QApplication(sys.argv)
    main_window = MainWindow()
    main_window.show()
    sys.exit(app.exec_())
