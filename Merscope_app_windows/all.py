import sys
import matplotlib
matplotlib.use('Qt5Agg')
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtWidgets import (QApplication, QMainWindow, QAction, QFileDialog, 
                             QMessageBox, QWidget, QVBoxLayout, QGridLayout, QGroupBox, 
                             QRadioButton, QLabel, QLineEdit, QPushButton, QProgressBar,
                             QListWidget, QListWidgetItem, QSplitter, QDialog)
import tempfile
from PyQt5.QtCore import Qt, QUrl, QThread, pyqtSignal
from numpy import array, cos, sin, dot, vstack, pi


class DataProcessingThread(QThread):
    # Signal for progress update and processing finish
    progress_changed = pyqtSignal(int)
    processing_finished = pyqtSignal(object)

    def __init__(self, file_name, parent=None):
        super().__init__(parent)
        self.file_name = file_name
        self.adata = None

    def run(self):
        from scanpy import read_h5ad
        # Emit signal for the start of file loading
        self.progress_changed.emit(10)
        try:
            # Load the file and update progress
            self.adata = read_h5ad(self.file_name)
            self.progress_changed.emit(50)
        except Exception as e:
            # In case of an error, emit signal with error object
            self.progress_changed.emit(0)
            self.processing_finished.emit(e)
            return

        self.progress_changed.emit(80)
        # Emit signal when processing is completed
        self.progress_changed.emit(100)
        self.processing_finished.emit(self.adata)


class MainWindow(QMainWindow):

    def __init__(self):
        super().__init__()
        self.adata = None  # Holds the loaded data
        self.progress_bar = None  # Progress bar for file loading
        self.initUI()

    def initUI(self):
        # Initialize the main UI layout and components
        self.setup_menu_bar()
        self.setup_status_bar()
        self.setup_main_layout()
        self.setWindowTitle('PyQt5 Application')
        self.setGeometry(100, 100, 1200, 800)

    def setup_menu_bar(self):
        # Setup the menu bar with 'File' and 'Tool' menus
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

    def setup_status_bar(self):
        # Setup the status bar with a progress bar
        self.statusbar = self.statusBar()
        self.statusbar.showMessage('Ready')
        self.progress_bar = QProgressBar(self)
        self.statusbar.addPermanentWidget(self.progress_bar)
        self.progress_bar.setValue(0)

    def setup_main_layout(self):
        # Setup the main layout of the window with a horizontal splitter
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)

        main_splitter = QSplitter(Qt.Horizontal)
        self.central_widget.setLayout(QVBoxLayout())
        self.central_widget.layout().addWidget(main_splitter)

        self.setup_left_area(main_splitter)
        self.setup_central_area(main_splitter)
        
        main_splitter.setSizes([200, 1000])

    def setup_left_area(self, parent_splitter):
        # Setup the left panel with spatial plot and gene list options
        left_panel = QGroupBox()
        left_layout = QGridLayout(left_panel)

        # Add sections to the left layout
        self.create_group_box('Spatial Plot', left_layout, 0, 0, 1, 2, self.create_spatial_plot_controls)
        self.create_group_box('Genes', left_layout, 1, 0, 1, 2, self.create_genes_list)
        self.create_group_box('Plotting', left_layout, 2, 0, 1, 2, self.create_plotting_controls)

        parent_splitter.addWidget(left_panel)

    def setup_central_area(self, parent_splitter):
        # Setup the central area for displaying web views or plots
        central_panel = QWidget()
        self.central_layout = QVBoxLayout(central_panel)
        self.web_view = QWebEngineView(self)
        self.central_layout.addWidget(self.web_view)
        parent_splitter.addWidget(central_panel)

    def create_group_box(self, title, layout, row, col, rowspan, colspan, setup_func):
        # Helper function to create a group box in the layout
        group_box = QGroupBox(title)
        group_layout = QGridLayout(group_box)
        setup_func(group_layout)
        layout.addWidget(group_box, row, col, rowspan, colspan)

    def create_spatial_plot_controls(self, layout):
        # Controls for spatial plot configuration
        axis_label = QLabel('Axis :')
        self.axis_input = QLineEdit('0')
        layout.addWidget(axis_label, 0, 0)
        layout.addWidget(self.axis_input, 0, 1)

        start_button = QPushButton('Start')
        start_button.clicked.connect(self.plot_spatial_scatter)
        layout.addWidget(start_button, 1, 0, 1, 2)

    def create_genes_list(self, layout):
        # Gene list widget for selecting genes to analyze
        self.genes_list_widget = QListWidget()
        layout.addWidget(self.genes_list_widget, 0, 0, 1, 2)

        # Adding a "Clear" button under the Genes list
        clear_button = QPushButton('Clear')
        clear_button.clicked.connect(self.clear_gene_selection)
        layout.addWidget(clear_button, 1, 0, 1, 2)

    def create_plotting_controls(self, layout):
        # Plotting control options for different plot types
        self.dot_plot_radio_button = QRadioButton('Dot plot')
        self.violin_plot_radio_button = QRadioButton('Violin plot')
        self.feature_plot_radio_button = QRadioButton('Feature plot')
        self.scatter_plot_radio_button = QRadioButton('Scatter plot')

        start_button = QPushButton('Start')
        start_button.clicked.connect(self.plot_genes)

        layout.addWidget(self.dot_plot_radio_button, 0, 0)
        layout.addWidget(self.violin_plot_radio_button, 0, 1)
        layout.addWidget(self.feature_plot_radio_button, 1, 0)
        layout.addWidget(self.scatter_plot_radio_button, 1, 1)
        layout.addWidget(start_button, 2, 0, 1, 2)

    def open_file(self):
        # Open file dialog to load an h5ad file
        if self.adata is not None:
            reply = QMessageBox.question(self, 'Warning', "Loading a new file will overwrite the existing data. Do you want to continue?", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if reply == QMessageBox.No:
                return

        file_name, _ = QFileDialog.getOpenFileName(self, "Open HDF5 File", "", "HDF5 Files (*.h5ad *.hdf5);;All Files (*)")
        if file_name:
            self.start_data_processing(file_name)

    def start_data_processing(self, file_name):
        # Start the file loading and processing in a separate thread
        self.progress_bar.setValue(0)
        self.data_thread = DataProcessingThread(file_name)
        self.data_thread.progress_changed.connect(self.update_progress)
        self.data_thread.processing_finished.connect(self.on_data_loaded)
        self.data_thread.start()

    def update_progress(self, value):
        # Update the progress bar with current progress value
        self.progress_bar.setValue(value)

    def on_data_loaded(self, result):
        # Callback when data loading is finished
        if isinstance(result, Exception):
            QMessageBox.critical(self, "Error", f"Failed to load data: {str(result)}")
        else:
            self.adata = result
            self.update_ui_after_data_loaded()

    def update_ui_after_data_loaded(self):
        # Update the UI after data is loaded successfully
        self.update_list(self.genes_list_widget, is_gene=True)
        self.plot_spatial_scatter()
        self.deg_action.setEnabled(True)
        self.progress_bar.setValue(100)
        self.statusbar.showMessage("Data loaded successfully.")

    def update_list(self, list_widget, is_gene=False):
        # Update the list widget with gene or cluster data
        list_widget.clear()
        items = sorted(self.adata.var_names.to_list()) if is_gene else sorted(self.adata.obs['cell_type_2'].unique())
        for item in items:
            list_item = QListWidgetItem(item)
            list_item.setFlags(list_item.flags() | Qt.ItemIsUserCheckable)
            list_item.setCheckState(Qt.Unchecked)
            list_widget.addItem(list_item)
    
    def clear_gene_selection(self):
        for i in range(self.genes_list_widget.count()):
            item = self.genes_list_widget.item(i)
            item.setCheckState(Qt.Unchecked)

    def plot_spatial_scatter(self):
        import gc
        # Plot the spatial scatter plot using Plotly
        if self.adata is None:
            return
        
        import plotly.graph_objects as go
        axis = self.get_rotation_axis()

        try:
            # Retrieve data for plotting
            color_categories = self.adata.obs['cell_type_2'].astype('category').cat.categories
            color_map = {category: color for category, color in zip(color_categories, self.adata.uns['cell_type_2_colors'])}

            cell_type_categories_sorted = sorted(self.adata.obs['cell_type_2'].unique(), key=str)

            rotated_x_coords, rotated_y_coords = self.rotate_coords(axis)
            fig = go.Figure()

            # Plot each category as a separate trace
            for category in cell_type_categories_sorted:
                mask = self.adata.obs['cell_type_2'] == category
                fig.add_trace(go.Scattergl(
                    x=rotated_x_coords[mask], y=rotated_y_coords[mask],
                    mode='markers', name=str(category),
                    marker=dict(size=2, color=color_map[category])
                ))

            # Update plot layout
            fig.update_layout(
                    xaxis=dict(scaleanchor="y", scaleratio=1, range=[rotated_x_coords.min(), rotated_x_coords.max()], showticklabels=False),
                    yaxis=dict(scaleanchor="x", scaleratio=1, range=[rotated_y_coords.min(), rotated_y_coords.max()], showticklabels=False),
                    legend=dict(font=dict(size=15), itemsizing="constant", title=dict(font=dict(size=30))),
                    plot_bgcolor="white", paper_bgcolor="white"
                )

            self.display_plot(fig)

        except Exception as e:
            QMessageBox.warning(self, 'Error', f'Could not plot data:\n{str(e)}')
        finally:
            gc.collect()

    def plot_genes(self):
        # Plot selected genes using matplotlib
        from matplotlib.pyplot import subplots, colorbar, show
        import gc

        if self.adata is None:
            return

        selected_genes = [item.text() for item in self.genes_list_widget.findItems("*", Qt.MatchWildcard) if item.checkState() == Qt.Checked]
        if not selected_genes:
            QMessageBox.warning(self, "Warning", "Please select at least one gene.")
            return

        if self.scatter_plot_radio_button.isChecked():
            axis = self.get_rotation_axis()
            rotated_x_coords, rotated_y_coords = self.rotate_coords(axis)

            num_genes = len(selected_genes)
            rows, cols = (1, num_genes) if num_genes <= 4 else ((num_genes + 3) // 4, 4)

            # Create subplots for multiple genes
            fig, axs = subplots(rows, cols, figsize=(cols * 5, rows * 5))
            axs = axs.flat if num_genes > 1 else [axs]

            # Plot each gene in a separate subplot
            for i, gene in enumerate(selected_genes):
                gene_expression = self.adata[:, gene].X.toarray().flatten()
                ax = axs[i]
                scatter = ax.scatter(x=rotated_x_coords, y=rotated_y_coords, c=gene_expression, cmap='Purples', s=0.5, edgecolors='none')
                ax.set_title(gene)
                ax.axis('off')
                ax.set_aspect('equal', 'box')
                colorbar(scatter, ax=ax)

            for j in range(i + 1, len(axs)):
                fig.delaxes(axs[j])
        else:
            # Use Scanpy for other plot types (dot, violin, feature plot)
            from scanpy import plotting as pl
            if self.dot_plot_radio_button.isChecked():
                pl.dotplot(self.adata, selected_genes, groupby="cell_type_2", categories_order=sorted(self.adata.obs['cell_type_2'].unique()))
            elif self.violin_plot_radio_button.isChecked():
                pl.stacked_violin(self.adata, selected_genes, groupby="cell_type_2", categories_order=sorted(self.adata.obs['cell_type_2'].unique()))
            elif self.feature_plot_radio_button.isChecked():
                selected_genes.append('cell_type_2')
                pl.umap(self.adata, color=selected_genes, ncols=4)
        show()
        gc.collect()

    def get_rotation_axis(self):
        # Get rotation axis for spatial scatter plot
        num = float(self.axis_input.text())
        return 180 / num if num > 0 else 360

    def rotate_coords(self, axis):
        # Rotate coordinates for spatial scatter plot
        x_coords = self.adata.obsm['spatial'][:, 0]
        y_coords = self.adata.obsm['spatial'][:, 1]
        center_x = (x_coords.max() + x_coords.min()) / 2
        center_y = (y_coords.max() + y_coords.min()) / 2

        rotation_matrix = array([[cos(pi / axis), -sin(pi / axis)], [sin(pi / axis), cos(pi / axis)]])
        rotated_coords = dot(vstack([x_coords - center_x, y_coords - center_y]).T, rotation_matrix)
        return rotated_coords[:, 0] + center_x, rotated_coords[:, 1] + center_y

    def display_plot(self, fig):
        # Display the generated plot in the central web view
        with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as temp_file:
            fig.write_html(temp_file.name)
            self.web_view.setUrl(QUrl.fromLocalFile(temp_file.name))

    def open_deg_analysis(self):
        # Open the Differential Expression Gene (DEG) Analysis dialog
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
        # Initialize the UI for DEG Analysis
        self.setWindowTitle('DEG Analysis')
        self.setGeometry(200, 200, 600, 400)

        layout = QVBoxLayout()
        clusters_layout = QGridLayout()

        self.clusters_list_widget = QListWidget()
        self.other_clusters_list_widget = QListWidget()

        # Create group boxes for clusters and other clusters selection
        self.create_group_box('Clusters', clusters_layout, 0, 0, self.populate_clusters_list, self.clusters_list_widget)
        self.create_group_box('Other Clusters', clusters_layout, 0, 1, self.populate_other_clusters_list, self.other_clusters_list_widget)

        layout.addLayout(clusters_layout)
        layout.addWidget(self.create_button('Clear', self.clear_selection))
        layout.addWidget(self.create_button('Start', self.run_deg_analysis))
        layout.addWidget(self.create_button('Save', self.save_deg_result))
        self.setLayout(layout)

    def create_group_box(self, title, layout, row, col, populate_func, list_widget):
        # Create group box and populate it with clusters or other clusters
        group_box = QGroupBox(title)
        group_layout = QVBoxLayout(group_box)
        populate_func(list_widget)
        group_layout.addWidget(list_widget)
        layout.addWidget(group_box, row, col)

    def create_button(self, text, callback):
        # Helper function to create buttons
        button = QPushButton(text)
        button.clicked.connect(callback)
        return button

    def populate_clusters_list(self, list_widget):
        # Populate the clusters list with unique cluster values
        self.populate_list(list_widget, self.adata.obs['cell_type_2'].unique())

    def populate_other_clusters_list(self, list_widget):
        # Populate the other clusters list with unique cluster values
        self.populate_list(list_widget, self.adata.obs['cell_type_2'].unique())

    def populate_list(self, list_widget, items):
        # Helper function to populate a list widget with items
        list_widget.clear()
        for item in sorted(items):
            list_item = QListWidgetItem(item)
            list_item.setFlags(list_item.flags() | Qt.ItemIsUserCheckable)
            list_item.setCheckState(Qt.Unchecked)
            list_widget.addItem(list_item)

    def clear_selection(self):
        # Clear the selection in both clusters and other clusters lists
        self.clear_list(self.clusters_list_widget)
        self.clear_list(self.other_clusters_list_widget)

    def clear_list(self, list_widget):
        # Clear the selection in a specific list widget
        for i in range(list_widget.count()):
            item = list_widget.item(i)
            item.setCheckState(Qt.Unchecked)

    def run_deg_analysis(self):     
        # Run the Differential Expression Gene (DEG) analysis
        selected_clusters = self.get_selected_items(self.clusters_list_widget)
        selected_other_clusters = self.get_selected_items(self.other_clusters_list_widget)

        if not selected_clusters or not selected_other_clusters:
            QMessageBox.warning(self, "Warning", "Please select at least one cluster in both Clusters and Other Clusters.")
            return

        self.run_deg_analysis_logic(selected_clusters, selected_other_clusters)
        QMessageBox.information(self, "Success", "DEG Analysis completed successfully.")

    def run_deg_analysis_logic(self, selected_clusters, selected_other_clusters):
        # Execute the actual DEG analysis logic using Scanpy
        from scanpy import tools as tl
        from scanpy import plotting as pl

        merged_cluster_A = '_'.join(selected_clusters)
        merged_cluster_B = '_'.join(selected_other_clusters)

        self.adata.obs['cell_type_3'] = self.adata.obs['cell_type_2'].astype(str)
        self.adata.obs['cell_type_3'][self.adata.obs['cell_type_3'].isin(selected_clusters)] = merged_cluster_A
        self.adata.obs['cell_type_3'][self.adata.obs['cell_type_3'].isin(selected_other_clusters)] = merged_cluster_B                        
        self.adata.obs['cell_type_3'] = self.adata.obs['cell_type_3'].astype('category')

        tl.rank_genes_groups(self.adata, 'cell_type_3', groups=[merged_cluster_A], reference=merged_cluster_B, method='wilcoxon')            
        pl.rank_genes_groups(self.adata, n_genes=25, sharey=False)

    def get_selected_items(self, list_widget):
        # Get selected items from the list widget
        return [item.text() for item in list_widget.findItems("*", Qt.MatchWildcard) if item.checkState() == Qt.Checked]

    def save_deg_result(self):
        # Save the DEG analysis result to a CSV file
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
                self.statusBar().showMessage(f"Data saved: {filePath}")
        else:
            QMessageBox.critical(self, "Error", "No results to save. Please run the analysis first.")


if __name__ == '__main__':
    from multiprocessing import freeze_support
    freeze_support()

    app = QApplication(sys.argv)
    main_window = MainWindow()
    main_window.show()
    sys.exit(app.exec_())
