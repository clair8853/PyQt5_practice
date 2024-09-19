import sys
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtWidgets import (QApplication, QMainWindow, QAction, QFileDialog, 
                             QMessageBox, QWidget, QVBoxLayout, QHBoxLayout, QGroupBox, 
                             QRadioButton, QLabel, QLineEdit, QPushButton,
                             QListWidget, QListWidgetItem, QSplitter)
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
        left_layout = QVBoxLayout(left_panel)

        # Adding the Spatial Plot group box under Clusters
        spatial_plot_group_box = QGroupBox('Spatial Plot')
        spatial_plot_layout = QVBoxLayout(spatial_plot_group_box)
        left_layout.addWidget(spatial_plot_group_box)

        # Second row: Axis label and input box
        axis_layout = QHBoxLayout()
        axis_label = QLabel('Axis :')
        self.axis_input = QLineEdit()
        self.axis_input.setText('0')  # Default value is 0
        axis_layout.addWidget(axis_label)
        axis_layout.addWidget(self.axis_input)
        spatial_plot_layout.addLayout(axis_layout)

        # Third row: Start button (aligned to the left)
        start_button = QPushButton('Start')
        start_button.clicked.connect(self.plot_spatial_scatter)
        spatial_plot_layout.addWidget(start_button)

        # Create Genes group box
        genes_group_box = QGroupBox('Genes')
        genes_layout = QVBoxLayout(genes_group_box)
        left_layout.addWidget(genes_group_box)

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
        left_layout.addWidget(genes_plotting_group_box)

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
            reply = QMessageBox.question(self, 'Warning', "Loading a new file will overwrite the existing data. Do you want to continue?", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)

            if reply == QMessageBox.No:
                return  # Exit the function if the user chooses not to continue

        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getOpenFileName(self, "Open HDF5 File", "", "HDF5 Files (*.h5ad *.hdf5);;All Files (*)", options=options)

        if file_name:
            self.load_data(file_name)

    def load_data(self, file_name):
        from scanpy import read_h5ad
        try:
            # Load the data using scanpy's read_h5ad function
            self.adata = read_h5ad(file_name)
            self.update_ui_after_data_loaded(file_name)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load file: {str(e)}")
            print(f"Error loading file: {str(e)}")

    def update_ui_after_data_loaded(self, file_name):
        # UI 업데이트 (데이터 로딩 후 처리)
        self.update_lists()
        self.plot_spatial_scatter()

        base_name = os.path.basename(file_name)
        self.statusbar.showMessage(f"Loaded '{base_name}'")
        print(f'Data loaded from: {base_name}')

    def update_lists(self):
        # Update List
        if self.adata is not None:
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

    def plot_spatial_scatter(self):
        import plotly.express as px
        
        axis = 360
        num = float(self.axis_input.text())
        if num > 0:
            axis = 180/num
        
        if self.adata is not None:
            try:
                # Plotting code using anndata data
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

                # 전체 데이터에 대한 x, y 축의 범위 설정
                x_range = [rotated_x_coords.min(), rotated_x_coords.max()]
                y_range = [rotated_y_coords.min(), rotated_y_coords.max()]

                fig = px.scatter(x=rotated_x_coords,
                                 y=rotated_y_coords,
                                 color=self.adata.obs['cell_type_2'],
                                 labels={'color': 'Clusters'})

                fig.update_traces(marker=dict(size=2))

                fig.update_layout(
                    xaxis=dict(
                        scaleanchor="y",
                        scaleratio=1,
                        range=x_range,  # x축 범위 고정
                        showticklabels=False  # x축 라벨과 수치 숨기기
                    ),
                    yaxis=dict(
                        scaleanchor="x",
                        scaleratio=1,
                        range=y_range,  # y축 범위 고정
                        showticklabels=False  # y축 라벨과 수치 숨기기
                    ),
                    plot_bgcolor="white",
                    paper_bgcolor="white"
                )

                # Plot을 HTML로 변환하여 QWebEngineView에 표시
                with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as temp_file:
                    fig.write_html(temp_file.name)
                    self.web_view.setUrl(QUrl.fromLocalFile(temp_file.name))

            except Exception as e:
                # Plotting 실패 시 경고 메시지 표시
                QMessageBox.warning(self, 'Error', f'Could not plot data:\n{str(e)}')
        else:
            QMessageBox.warning(self, 'Error', 'No data loaded to plot.')

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

            if self.scatter_plot_radio_button.isChecked():
                from math import ceil, sqrt

                # Plotting code using anndata data
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

                fig, axes = subplots(1, num_genes, figsize=(5*num_genes,5))
                
                if num_genes == 1:
                    axes = [axes]
                
                for gene, ax in zip(selected_genes, axes):
                    gene_expression = self.adata[:, gene].X.toarray().flatten()
                    scatter = ax.scatter(x=rotated_x_coords,
                                         y=rotated_y_coords,
                                         c=gene_expression,
                                         cmap='Purples', s=1)
                    ax.set_title(f'{gene}')
                    ax.axis('off')
                    colorbar(scatter, ax=ax)
                
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

if __name__ == '__main__':
    from multiprocessing import freeze_support
    freeze_support()

    app = QApplication(sys.argv)
    main_window = MainWindow()
    main_window.show()
    sys.exit(app.exec_())