import sys
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtWidgets import (QApplication, QMainWindow, QAction, QFileDialog, 
                             QMessageBox, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QGroupBox, 
                             QRadioButton, QLabel, QLineEdit, QPushButton,
                             QListWidget, QListWidgetItem, QSplitter, QDialog)
import tempfile
from PyQt5.QtCore import Qt, QUrl
from numpy import array, cos, sin, dot, vstack, pi
import os

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

        # Adding 'Tool' menu
        tool_menu = menubar.addMenu('&Tool')

        # Adding 'DEG Analysis' action
        self.deg_action = QAction('DEG Analysis', self)
        self.deg_action.setShortcut('Ctrl+A')
        self.deg_action.setStatusTip('Differential Gene Expression Analysis')
        self.deg_action.triggered.connect(self.open_deg_analysis)
        self.deg_action.setEnabled(False)
        tool_menu.addAction(self.deg_action)

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
        
        # Set the split ratios (1:3:1)
        main_splitter.setSizes([200, 1000])

        # Main window settings
        self.setWindowTitle('PyQt5 Application')
        self.setGeometry(100, 100, 1200, 800)

    def setup_left_area(self, parent_splitter):
        # Create the left area
        left_panel = QGroupBox()
        left_layout = QGridLayout(left_panel)  # QVBoxLayout -> QGridLayout으로 변경

        # Adding the Spatial Plot group box under Clusters
        spatial_plot_group_box = QGroupBox('Spatial Plot')
        spatial_plot_layout = QGridLayout(spatial_plot_group_box)
        left_layout.addWidget(spatial_plot_group_box, 0, 0, 1, 2)

        # Axis label and input box
        axis_label = QLabel('Axis :')
        self.axis_input = QLineEdit()
        self.axis_input.setText('0')  # Default value is 0
        spatial_plot_layout.addWidget(axis_label, 0, 0)
        spatial_plot_layout.addWidget(self.axis_input, 0, 1)

        # Start button
        start_button = QPushButton('Start')
        start_button.clicked.connect(self.plot_spatial_scatter)
        spatial_plot_layout.addWidget(start_button, 1, 0, 1, 2)

        # Create Genes group box
        genes_group_box = QGroupBox('Genes')
        genes_layout = QGridLayout(genes_group_box)
        left_layout.addWidget(genes_group_box, 1, 0, 1, 2)

        # Create the list widget for displaying genes
        self.genes_list_widget = QListWidget()
        genes_layout.addWidget(self.genes_list_widget, 0, 0, 1, 2)

        # Create Plotting radio buttons inside a group box
        genes_plotting_group_box = QGroupBox('Plotting')
        genes_plotting_layout = QGridLayout(genes_plotting_group_box)
        self.dot_plot_radio_button = QRadioButton('Dot plot')
        self.violin_plot_radio_button = QRadioButton('Violin plot')
        self.feature_plot_radio_button = QRadioButton('Feature plot')
        self.scatter_plot_radio_button = QRadioButton('Scatter plot')

        start_button = QPushButton('Start')
        start_button.clicked.connect(self.plot_genes)

        # Grid에 라디오 버튼 추가
        genes_plotting_layout.addWidget(self.dot_plot_radio_button, 0, 0)
        genes_plotting_layout.addWidget(self.violin_plot_radio_button, 0, 1)
        genes_plotting_layout.addWidget(self.feature_plot_radio_button, 1, 0)
        genes_plotting_layout.addWidget(self.scatter_plot_radio_button, 1, 1)
        genes_plotting_layout.addWidget(start_button, 2, 0, 1, 2)

        left_layout.addWidget(genes_plotting_group_box, 2, 0, 1, 2)

        # Add the left panel to the parent splitter
        parent_splitter.addWidget(left_panel)

    def setup_central_area(self, parent_splitter):
        # Create the central area
        central_panel = QWidget()
        self.central_layout = QVBoxLayout(central_panel)
        self.web_view = QWebEngineView(self)
        self.central_layout.addWidget(self.web_view)
        parent_splitter.addWidget(central_panel)
    
    def open_file(self):
        if self.adata is not None:
            reply = QMessageBox.question(
                self, 'Warning', 
                "Loading a new file will overwrite the existing data. Do you want to continue?", 
                QMessageBox.Yes | QMessageBox.No, 
                QMessageBox.No
            )
            if reply == QMessageBox.No:
                return  # 사용자가 'No'를 선택하면 파일 로딩을 중단합니다.

        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getOpenFileName(self, "Open HDF5 File", "", "HDF5 Files (*.h5ad *.hdf5);;All Files (*)", options=options)

        if file_name:
            self.load_data(file_name)

    def load_data(self, file_name):
        from scanpy import read_h5ad
        try:
            self.adata = read_h5ad(file_name)
            self.update_ui_after_data_loaded(file_name)
            self.deg_action.setEnabled(True)
            QMessageBox.information(self, "Success", "Data loaded successfully.")  # 사용자 피드백 추가
        except IOError:
            QMessageBox.critical(self, "Error", "File not found or unable to open.")
        except ValueError as ve:
            QMessageBox.critical(self, "Error", f"Data format error: {str(ve)}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Unexpected error: {str(e)}")

    def update_ui_after_data_loaded(self, file_name):
        self.update_lists()
        self.plot_spatial_scatter()
        base_name = os.path.basename(file_name)
        self.statusbar.showMessage(f"Loaded '{base_name}'")
        print(f'Data loaded from: {base_name}')

    def update_lists(self):
        if self.adata is not None:
            self.update_list(self.genes_list_widget, is_gene=True)
        else:
            self.disable_ui_elements()  # self.adata가 None일 때 UI 비활성화

    def update_list(self, list_widget, is_gene=False):
        list_widget.clear()
        items = sorted(self.adata.var_names.to_list()) if is_gene else sorted(self.adata.obs['cell_type_2'].unique())
        for item in items:
            list_item = QListWidgetItem(item)
            list_item.setFlags(list_item.flags() | Qt.ItemIsUserCheckable)
            list_item.setCheckState(Qt.Unchecked)
            list_widget.addItem(list_item)

    def disable_ui_elements(self):
        # UI 비활성화 (데이터 로드 실패 시)
        self.deg_action.setEnabled(False)
        self.statusbar.showMessage('Ready')

    def plot_spatial_scatter(self):
        import plotly.graph_objects as go
        import gc  # 메모리 관리용
               
        axis = 360
        num = float(self.axis_input.text())
        if num > 0:
            axis = 180/num
        
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

                x_range = [rotated_x_coords.min(), rotated_x_coords.max()]
                y_range = [rotated_y_coords.min(), rotated_y_coords.max()]

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
                    xaxis=dict(
                        scaleanchor="y",
                        scaleratio=1,
                        range=x_range,
                        showticklabels=False,
                        title = ''
                    ),
                    yaxis=dict(
                        scaleanchor="x",
                        scaleratio=1,
                        range=y_range,
                        showticklabels=False,
                        title = ''
                    ),
                    legend=dict(
                        font=dict(
                            size=15
                        ),
                        itemsizing="constant",
                        title=dict(
                            font=dict(
                                size=30
                            )
                        )
                    ),
                    legend_title_text = 'Clusters',
                    plot_bgcolor="white",
                    paper_bgcolor="white"
                )

                with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as temp_file:
                    fig.write_html(temp_file.name)
                    self.web_view.setUrl(QUrl.fromLocalFile(temp_file.name))

            except Exception as e:
                QMessageBox.warning(self, 'Error', f'Could not plot data:\n{str(e)}')
            finally:
                gc.collect()  # 메모리 관리

        else:
            QMessageBox.warning(self, 'Error', 'No data loaded to plot.')

    def plot_genes(self):
        from matplotlib.pyplot import subplots, colorbar, show
        import gc  # 메모리 관리용
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
                from math import ceil

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

                if num_genes <= 4:
                    rows = 1
                    cols = num_genes
                else:
                    rows = ceil(num_genes / 4)
                    cols = 4
                
                if num_genes > 1:
                    fig, axs = subplots(rows, cols, figsize=(cols*5, rows*5))
                    axs = axs.flat
                else:
                    fig, ax = subplots(1, 1, figsize=(5,5))
                    axs = [ax]
                
                for i, gene in enumerate(selected_genes):
                    ax = axs[i]
                    gene_expression = self.adata[:, gene].X.toarray().flatten()
                    scatter = ax.scatter(x=rotated_x_coords,
                                         y=rotated_y_coords,
                                         c=gene_expression,
                                         cmap='Purples', s=1)
                    ax.set_title(f'{gene}')
                    ax.axis('off')
                    ax.set_aspect('equal', 'box')
                    colorbar(scatter, ax=ax)

                for j in range(i+1, len(axs)):
                    fig.delaxes(axs[j]) 
                show()
                gc.collect()  # 메모리 관리
            else:
                from scanpy import plotting as pl
                if self.dot_plot_radio_button.isChecked():
                    pl.dotplot(self.adata, selected_genes, groupby="cell_type_2", categories_order=sorted(self.adata.obs['cell_type_2'].unique()))
                elif self.violin_plot_radio_button.isChecked():
                    pl.stacked_violin(self.adata, selected_genes, groupby="cell_type_2", categories_order=sorted(self.adata.obs['cell_type_2'].unique()))
                elif self.feature_plot_radio_button.isChecked():
                    selected_genes.append('cell_type_2')
                    pl.umap(self.adata, color=selected_genes, ncols=4)
                show()
                gc.collect()  # 메모리 관리

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

        # Horizontal layout for Clusters and Other Clusters
        clusters_layout = QGridLayout()  # QHBoxLayout -> QGridLayout으로 변경

        # Clusters selection
        clusters_group_box = QGroupBox('Clusters')
        clusters_vbox = QVBoxLayout(clusters_group_box)
        self.clusters_list_widget = QListWidget()
        self.populate_list(self.clusters_list_widget)
        clusters_vbox.addWidget(self.clusters_list_widget)
        clusters_layout.addWidget(clusters_group_box, 0, 0)

        # Other Clusters selection
        other_clusters_group_box = QGroupBox('Other Clusters')
        other_clusters_vbox = QVBoxLayout(other_clusters_group_box)
        self.other_clusters_list_widget = QListWidget()
        self.populate_list(self.other_clusters_list_widget)
        other_clusters_vbox.addWidget(self.other_clusters_list_widget)
        clusters_layout.addWidget(other_clusters_group_box, 0, 1)

        # Sync the 'Clusters' and 'Other Clusters' list widget
        self.clusters_list_widget.itemChanged.connect(self.sync_lists)
        self.other_clusters_list_widget.itemChanged.connect(self.sync_lists)

        layout.addLayout(clusters_layout)

        # Clear button
        clear_button = QPushButton('Clear')
        clear_button.clicked.connect(lambda: self.clear_selection(self.clusters_list_widget))
        clear_button.clicked.connect(lambda: self.clear_selection(self.other_clusters_list_widget))
        layout.addWidget(clear_button)
        
        # Start button
        start_button = QPushButton('Start')
        start_button.clicked.connect(self.run_deg_analysis)
        layout.addWidget(start_button)

        # Save button
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

        self.selected_clusters = [item.text() for item in self.clusters_list_widget.findItems("*", Qt.MatchWildcard) if item.checkState() == Qt.Checked]
        self.selected_other_clusters = [item.text() for item in self.other_clusters_list_widget.findItems("*", Qt.MatchWildcard) if item.checkState() == Qt.Checked]

        if not self.selected_clusters or not self.selected_other_clusters:
            QMessageBox.warning(self, "Warning", "Please select at least one cluster in both Clusters and Other Clusters.")
            return
        
        else:
            merged_cluster_A = '_'.join(self.selected_clusters)
            merged_cluster_B = '_'.join(self.selected_other_clusters)
                        
            self.adata.obs['cell_type_3'] = self.adata.obs['cell_type_2']
            self.adata.obs['cell_type_3'] = self.adata.obs['cell_type_3'].astype(str)
            self.adata.obs['cell_type_3'][self.adata.obs['cell_type_3'].isin(self.selected_clusters)] = merged_cluster_A
            self.adata.obs['cell_type_3'][self.adata.obs['cell_type_3'].isin(self.selected_other_clusters)] = merged_cluster_B                        
            self.adata.obs['cell_type_3'] = self.adata.obs['cell_type_3'].astype('category')
            
            tl.rank_genes_groups(self.adata, 'cell_type_3', groups=[merged_cluster_A], reference=merged_cluster_B, method='wilcoxon')            
            pl.rank_genes_groups(self.adata, n_genes=25, sharey=False)
            QMessageBox.information(self, "Success", "DEG Analysis completed successfully.")  # 사용자 피드백 추가
    
    def save_deg_result(self):
        from pandas import DataFrame
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        filePath, _ = QFileDialog.getSaveFileName(self, "Save CSV File", "", "CSV Files (*.csv);;All Files (*)", options=options)
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
            QMessageBox.critical(self, "Error", f"No results to save. Please run the analysis first.")

if __name__ == '__main__':
    from multiprocessing import freeze_support
    freeze_support()

    app = QApplication(sys.argv)
    main_window = MainWindow()
    main_window.show()
    sys.exit(app.exec_())
