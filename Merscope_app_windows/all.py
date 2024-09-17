import sys
import matplotlib
matplotlib.use('Qt5Agg')
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtWidgets import (QApplication, QMainWindow, QAction, QFileDialog, QHBoxLayout,
                             QMessageBox, QWidget, QVBoxLayout, QGridLayout, QGroupBox, 
                             QRadioButton, QLabel, QLineEdit, QPushButton, QProgressBar,
                             QListWidget, QListWidgetItem, QSplitter, QDialog, QFontDialog,
                             QTextBrowser, QComboBox)
from PyQt5.QtCore import Qt, QUrl, QThread, pyqtSignal
from PIL import Image, ImageOps
from numpy import array, cos, sin, dot, vstack, pi


class DataProcessingThread(QThread):
    # Thread for loading and processing data in the background
    progress_changed = pyqtSignal(int)  # Signal to update progress bar
    processing_finished = pyqtSignal(object)  # Signal emitted when processing is done

    def __init__(self, file_name, parent=None):
        super().__init__(parent)
        self.file_name = file_name
        self.adata = None

    def run(self):
        # Method that runs the data loading
        from scanpy import read_h5ad
        self.progress_changed.emit(10)  # Initial progress update
        try:
            # Load the file and update progress
            self.adata = read_h5ad(self.file_name)
            self.progress_changed.emit(50)
        except Exception as e:
            # Emit error signal if loading fails
            self.progress_changed.emit(0)
            self.processing_finished.emit(e)
            return
        self.progress_changed.emit(80)
        self.progress_changed.emit(100)
        self.processing_finished.emit(self.adata)


class MainWindow(QMainWindow):
    # Main application window

    def __init__(self):
        super().__init__()
        self.adata = None  # Loaded data
        self.progress_bar = None  # Progress bar for loading
        self.selected_metadata = 'cell_type_2'
        self.font_settings = {"family": "Arial", "size": 15}  # Default font for plots
        self.initUI()

    def initUI(self):
        # Initialize UI components
        self.setup_menu_bar()
        self.setup_status_bar()
        self.setup_main_layout()
        self.setWindowTitle('Merscope Visualizer')
        self.setGeometry(100, 100, 1600, 900)

    def setup_menu_bar(self):
        # Setup the menu bar
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

        self.image_action = QAction('RGB Image', self)
        self.image_action.setShortcut('Ctrl+I')
        self.image_action.triggered.connect(self.open_image_merge_dialog)
        tool_menu.addAction(self.image_action)

        self.comparison_action = QAction('DEG other Samples', self)
        self.comparison_action.setShortcut('Ctrl+P')
        self.comparison_action.triggered.connect(self.comparison_other_samples_dialog)
        self.comparison_action.setEnabled(False)
        tool_menu.addAction(self.comparison_action)

        font_menu = menubar.addMenu('&Font')

        layout_font_action = QAction('Main Layout Font', self)
        layout_font_action.triggered.connect(self.change_font_for_layout)
        font_menu.addAction(layout_font_action)

        font_action = QAction('Main Plot Font', self)
        font_action.triggered.connect(self.change_font_for_plot)
        font_menu.addAction(font_action)

    def setup_status_bar(self):
        # Setup the status bar with a progress bar
        self.statusbar = self.statusBar()
        self.statusbar.showMessage('Ready')
        self.progress_bar = QProgressBar(self)
        self.statusbar.addPermanentWidget(self.progress_bar)
        self.progress_bar.setValue(0)

    def setup_main_layout(self):
        # Setup the main layout of the window
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)

        main_splitter = QSplitter(Qt.Horizontal)
        self.central_widget.setLayout(QVBoxLayout())
        self.central_widget.layout().addWidget(main_splitter)

        self.setup_left_area(main_splitter)
        self.setup_central_area(main_splitter)
        
        main_splitter.setSizes([200, 1600])

    def setup_left_area(self, parent_splitter):
        # Setup the left panel with controls for plotting
        left_panel = QGroupBox()
        left_layout = QGridLayout(left_panel)

        self.create_group_box('Spatial Plot', left_layout, 0, 0, 1, 2, self.create_spatial_plot_controls)
        self.create_group_box('Genes', left_layout, 1, 0, 1, 2, self.create_genes_list)
        self.create_group_box('Plotting', left_layout, 2, 0, 1, 2, self.create_plotting_controls)

        parent_splitter.addWidget(left_panel)

    def setup_central_area(self, parent_splitter):
        # Setup the central area to display plots or web views
        central_panel = QWidget()
        self.central_layout = QVBoxLayout(central_panel)
        self.web_view = QWebEngineView(self)
        self.central_layout.addWidget(self.web_view)
        parent_splitter.addWidget(central_panel)

    def create_group_box(self, title, layout, row, col, rowspan, colspan, setup_func):
        # Create group boxes for UI sections
        group_box = QGroupBox(title)
        group_layout = QGridLayout(group_box)
        setup_func(group_layout)
        layout.addWidget(group_box, row, col, rowspan, colspan)

    def create_spatial_plot_controls(self, layout):
        # Create controls for spatial plot options
        axis_label = QLabel('Axis :')
        self.axis_input = QLineEdit('0')
        layout.addWidget(axis_label, 0, 0)
        layout.addWidget(self.axis_input, 0, 1)

        start_button = QPushButton('Start')
        start_button.clicked.connect(self.plot_spatial_scatter)
        layout.addWidget(start_button, 1, 0, 1, 2)

    def create_genes_list(self, layout):
        # Create gene list for selection
        self.search_bar = QLineEdit(self)
        self.search_bar.setPlaceholderText("Enter the target gene")
        self.search_bar.textChanged.connect(self.search_list)
        layout.addWidget(self.search_bar, 0, 0, 1, 2)

        self.genes_list_widget = QListWidget()
        layout.addWidget(self.genes_list_widget, 1, 0, 1, 2)

        clear_button = QPushButton('Clear')
        clear_button.clicked.connect(self.clear_gene_selection)
        layout.addWidget(clear_button, 2, 0, 1, 2)

    def create_plotting_controls(self, layout):
        # Create controls for selecting plot types
        self.dot_plot_radio_button = QRadioButton('Dot plot')
        self.violin_plot_radio_button = QRadioButton('Violin plot')
        self.feature_plot_radio_button = QRadioButton('Feature plot')
        self.scatter_plot_radio_button = QRadioButton('Scatter plot')

        start_button = QPushButton('Start')
        start_button.clicked.connect(self.plot_genes)

        layout.addWidget(self.dot_plot_radio_button, 0, 0)
        layout.addWidget(self.violin_plot_radio_button, 1, 0)
        layout.addWidget(self.feature_plot_radio_button, 2, 0)
        layout.addWidget(self.scatter_plot_radio_button, 3, 0)
        layout.addWidget(start_button, 4, 0, 1, 2)
    
    def change_font_for_layout(self):
        # Open dialog to change font for the layout
        font, ok = QFontDialog.getFont()
        if ok:
            self.setFont(font)
            self.update_layout_font(self.central_widget, font)
            self.metadata_selector.setFont(font)
    
    def update_layout_font(self, widget, font):
        # Apply the selected font to all layout widgets
        widget.setFont(font)
        for child in widget.findChildren(QWidget):
            child.setFont(font)

    def change_font_for_plot(self):
        # Open dialog to change font for plot legends
        font, ok = QFontDialog.getFont()
        if ok:
            self.font_settings["family"] = font.family()
            self.font_settings["size"] = font.pointSize()
            self.plot_spatial_scatter()

    def open_file(self):
        # Open file dialog to load a new file
        if self.adata is not None:
            reply = QMessageBox.question(self, 'Warning', "Loading a new file will overwrite the existing data. Do you want to continue?", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if reply == QMessageBox.No:
                return

        file_name, _ = QFileDialog.getOpenFileName(self, "Open HDF5 File", "", "HDF5 Files (*.h5ad *.hdf5);;All Files (*)")
        if file_name:
            self.current_file_path = file_name
            self.start_data_processing(file_name)

    def start_data_processing(self, file_name):
        # Start processing the file in a separate thread
        self.progress_bar.setValue(0)
        self.data_thread = DataProcessingThread(file_name)
        self.data_thread.progress_changed.connect(self.update_progress)
        self.data_thread.processing_finished.connect(self.on_data_loaded)
        self.data_thread.start()

    def update_progress(self, value):
        # Update the progress bar with current value
        self.progress_bar.setValue(value)

    def on_data_loaded(self, result):
        # Handle data after it is loaded
        if isinstance(result, Exception):
            QMessageBox.critical(self, "Error", f"Failed to load data: {str(result)}")
        else:
            self.adata = result
            self.update_ui_after_data_loaded()

    def update_ui_after_data_loaded(self):
        # Update the UI after data is successfully loaded
        self.all_genes = self.adata.var_names.to_list()
        self.update_list(self.genes_list_widget, is_gene=True)
        self.plot_spatial_scatter()
        self.deg_action.setEnabled(True)
        self.comparison_action.setEnabled(True)
        self.progress_bar.setValue(100)
        self.statusbar.showMessage("Data loaded successfully.")
        self.setup_metadata_selector()
    
    def search_list(self):
        # Filter genes in the list based on search input
        search_text = self.search_bar.text().lower()

        if not search_text:
            self.update_list_with_checked_items(self.genes_list_widget, self.all_genes)
            return

        filtered_items = [item for item in self.all_genes if search_text in item.lower()]
        self.update_list_with_checked_items(self.genes_list_widget, filtered_items)
    
    def update_list_with_checked_items(self, list_widget, items):
        # Update list widget with checked items
        checked_items = {list_widget.item(i).text(): list_widget.item(i).checkState() for i in range(list_widget.count())}
        
        list_widget.clear()
        for item in items:
            list_item = QListWidgetItem(item)
            list_item.setFlags(list_item.flags() | Qt.ItemIsUserCheckable)
            list_item.setCheckState(checked_items.get(item, Qt.Unchecked))
            list_widget.addItem(list_item)

    def update_list(self, list_widget, is_gene=False):
        # Update gene or cluster list in the widget
        list_widget.clear()
        items = sorted(self.adata.var_names.to_list()) if is_gene else sorted(self.adata.obs[self.selected_metadata].unique())
        for item in items:
            list_item = QListWidgetItem(item)
            list_item.setFlags(list_item.flags() | Qt.ItemIsUserCheckable)
            list_item.setCheckState(Qt.Unchecked)
            list_widget.addItem(list_item)
    
    def clear_gene_selection(self):
        # Clear all selected genes in the list
        for i in range(self.genes_list_widget.count()):
            item = self.genes_list_widget.item(i)
            item.setCheckState(Qt.Unchecked)
    
    def setup_metadata_selector(self):
        # Setup a dropdown to select metadata for plotting
        if self.adata is None:
            return
        
        if hasattr(self, 'metadata_selector'):
            self.central_widget.layout().removeWidget(self.metadata_selector)
            self.metadata_selector.deleteLater()
            self.metadata_selector = None

        metadata_layout = self.central_widget.layout()
        self.metadata_selector = QComboBox(self)
        self.metadata_selector.addItems(self.adata.obs.columns)
        self.metadata_selector.currentTextChanged.connect(self.on_metadata_selected)
        metadata_layout.insertWidget(0, self.metadata_selector)
        self.metadata_selector.setFont(self.font())

    def on_metadata_selected(self, selected_metadata):
        # Handle metadata selection changes
        self.selected_metadata = selected_metadata
        self.plot_spatial_scatter()

    def plot_spatial_scatter(self):
        # Plot spatial scatter plot using Plotly
        import gc
        import plotly.graph_objects as go
        from scanpy.plotting._utils import _set_default_colors_for_categorical_obs
        if self.adata is None:
            return
        
        axis = self.get_rotation_axis()

        try:
            if self.adata.obs[self.selected_metadata].dtype.name == 'category':
                _set_default_colors_for_categorical_obs(self.adata, self.selected_metadata)
            color_categories = self.adata.obs[self.selected_metadata].astype('category').cat.categories
            color_uns = self.selected_metadata + "_colors"
            color_map = {category: color for category, color in zip(color_categories, self.adata.uns[color_uns])}

            cell_type_categories_sorted = sorted(self.adata.obs[self.selected_metadata].unique(), key=str)

            rotated_x_coords, rotated_y_coords = self.rotate_coords(axis)
            fig = go.Figure()

            for category in cell_type_categories_sorted:
                mask = self.adata.obs[self.selected_metadata] == category
                fig.add_trace(go.Scattergl(
                    x=rotated_x_coords[mask], y=rotated_y_coords[mask],
                    mode='markers', name=str(category),
                    marker=dict(size=2, color=color_map[category])
                ))

            fig.update_layout(
                xaxis=dict(scaleanchor="y", scaleratio=1, range=[rotated_x_coords.min(), rotated_x_coords.max()], showticklabels=False),
                yaxis=dict(scaleanchor="x", scaleratio=1, range=[rotated_y_coords.min(), rotated_y_coords.max()], showticklabels=False),
                legend=dict(
                    font=dict(family=self.font_settings["family"], size=self.font_settings["size"]),
                    itemsizing="constant"
                ),
                plot_bgcolor="white", paper_bgcolor="white"
            )

            self.display_plot(fig)

        except Exception as e:
            QMessageBox.warning(self, 'Error', f'Could not plot data:\n{str(e)}')
        finally:
            gc.collect()
    
    def plot_genes(self):
        # Plot selected genes using Matplotlib or Scanpy
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

            fig, axs = matplotlib.pyplot.subplots(rows, cols, figsize=(cols * 5, rows * 5))
            axs = axs.flat if num_genes > 1 else [axs]

            for i, gene in enumerate(selected_genes):
                gene_expression = self.adata[:, gene].X.toarray().flatten()
                ax = axs[i]
                scatter = ax.scatter(x=rotated_x_coords, y=rotated_y_coords, c=gene_expression, cmap='Purples', s=0.1, edgecolors='none')
                ax.set_title(gene)
                ax.axis('off')
                ax.set_aspect('equal', 'box')
                matplotlib.pyplot.colorbar(scatter, ax=ax)

            for j in range(i + 1, len(axs)):
                fig.delaxes(axs[j])
        else:
            from scanpy import plotting as pl
            if self.dot_plot_radio_button.isChecked():
                pl.dotplot(self.adata, selected_genes, groupby=self.selected_metadata, categories_order=sorted(self.adata.obs[self.selected_metadata].unique()))
            elif self.violin_plot_radio_button.isChecked():
                pl.stacked_violin(self.adata, selected_genes, groupby=self.selected_metadata, categories_order=sorted(self.adata.obs[self.selected_metadata].unique()))
            elif self.feature_plot_radio_button.isChecked():
                selected_genes.append(self.selected_metadata)
                pl.umap(self.adata, color=selected_genes, ncols=4)
        matplotlib.pyplot.show()
        gc.collect()

    def get_rotation_axis(self):
        # Get the axis for rotating spatial coordinates
        num = float(self.axis_input.text())
        return 180 / num if num > 0 else 360

    def rotate_coords(self, axis):
        # Rotate the spatial coordinates for plotting
        x_coords = self.adata.obsm['spatial'][:, 0]
        y_coords = self.adata.obsm['spatial'][:, 1]
        center_x = (x_coords.max() + x_coords.min()) / 2
        center_y = (y_coords.max() + y_coords.min()) / 2

        rotation_matrix = array([[cos(pi / axis), -sin(pi / axis)], [sin(pi / axis), cos(pi / axis)]])
        rotated_coords = dot(vstack([x_coords - center_x, y_coords - center_y]).T, rotation_matrix)
        return rotated_coords[:, 0] + center_x, rotated_coords[:, 1] + center_y

    def display_plot(self, fig):
        # Display the generated plot in the web view
        import tempfile
        with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as temp_file:
            fig.write_html(temp_file.name)
            self.web_view.setUrl(QUrl.fromLocalFile(temp_file.name))

    def open_deg_analysis(self):
        # Open dialog for Differential Expression Gene (DEG) analysis
        if self.adata is None:
            QMessageBox.warning(self, 'Error', 'No data loaded to perform DEG analysis.')
            return

        dialog = DEGAnalysisDialog(self.adata, self.selected_metadata, self)
        dialog.show()

    def open_image_merge_dialog(self):
        # Open dialog for merging RGB images
        dialog = ImageMergeDialog(self)
        dialog.show()
    
    def comparison_other_samples_dialog(self):
        # Open dialog for comparing DEG across samples
        dialog = OtherSamplesDialog(self.adata, self.selected_metadata, self)
        dialog.show()

class DEGAnalysisDialog(QDialog):
    # Dialog for Differential Expression Gene (DEG) analysis
    def __init__(self, adata, selected_metadata, parent=None):
        super().__init__(parent)
        self.adata = adata
        self.selected_metadata = selected_metadata
        self.initUI()

    def initUI(self):
        # Initialize UI for DEG analysis dialog
        self.setWindowTitle('DEG Analysis')
        self.setGeometry(50, 50, 600, 600)

        layout = QVBoxLayout()
        clusters_layout = QGridLayout()

        self.clusters_list_widget = QListWidget()
        self.other_clusters_list_widget = QListWidget()

        self.create_group_box('Clusters', clusters_layout, 0, 0, self.populate_clusters_list, self.clusters_list_widget)
        self.create_group_box('Other Clusters', clusters_layout, 0, 1, self.populate_other_clusters_list, self.other_clusters_list_widget)

        layout.addLayout(clusters_layout)
        layout.addWidget(self.create_button('Clear', self.clear_selection))
        layout.addWidget(self.create_button('Start', self.run_deg_analysis))
        layout.addWidget(self.create_button('Save', self.save_deg_result))

        self.clusters_list_widget.itemChanged.connect(self.sync_lists)
        self.other_clusters_list_widget.itemChanged.connect(self.sync_lists)

        self.setLayout(layout)

    def create_group_box(self, title, layout, row, col, populate_func, list_widget):
        # Create a group box and populate it
        group_box = QGroupBox(title)
        group_layout = QVBoxLayout(group_box)
        populate_func(list_widget)
        group_layout.addWidget(list_widget)
        layout.addWidget(group_box, row, col)

    def create_button(self, text, callback):
        # Create a button and connect it to the callback
        button = QPushButton(text)
        button.clicked.connect(callback)
        return button

    def populate_clusters_list(self, list_widget):
        # Populate the clusters list with unique cluster values
        self.populate_list(list_widget, self.adata.obs[self.selected_metadata].unique())

    def populate_other_clusters_list(self, list_widget):
        # Populate the other clusters list
        self.populate_list(list_widget, self.adata.obs[self.selected_metadata].unique())

    def populate_list(self, list_widget, items):
        # Helper function to populate a list widget with items
        list_widget.clear()
        for item in sorted(items):
            list_item = QListWidgetItem(item)
            list_item.setFlags(list_item.flags() | Qt.ItemIsUserCheckable)
            list_item.setCheckState(Qt.Unchecked)
            list_widget.addItem(list_item)
    
    def sync_lists(self, item):
        # Synchronize the selections in clusters and other clusters lists
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

    def clear_selection(self):
        # Clear selection in both clusters lists
        self.clear_list(self.clusters_list_widget)
        self.clear_list(self.other_clusters_list_widget)

    def clear_list(self, list_widget):
        # Clear all items in the list
        for i in range(list_widget.count()):
            item = list_widget.item(i)
            item.setCheckState(Qt.Unchecked)

    def run_deg_analysis(self):     
        # Perform DEG analysis
        selected_clusters = self.get_selected_items(self.clusters_list_widget)
        selected_other_clusters = self.get_selected_items(self.other_clusters_list_widget)

        if not selected_clusters or not selected_other_clusters:
            QMessageBox.warning(self, "Warning", "Please select at least one cluster in both Clusters and Other Clusters.")
            return

        self.run_deg_analysis_logic(selected_clusters, selected_other_clusters)
        QMessageBox.information(self, "Success", "DEG Analysis completed successfully.")

    def run_deg_analysis_logic(self, selected_clusters, selected_other_clusters):
        # Core logic for DEG analysis using Scanpy
        from scanpy import tools as tl
        from scanpy import plotting as pl

        merged_cluster_A = '_'.join(selected_clusters)
        merged_cluster_B = '_'.join(selected_other_clusters)

        self.adata.obs['cell_type_3'] = self.adata.obs[self.selected_metadata].astype(str)
        self.adata.obs['cell_type_3'][self.adata.obs['cell_type_3'].isin(selected_clusters)] = merged_cluster_A
        self.adata.obs['cell_type_3'][self.adata.obs['cell_type_3'].isin(selected_other_clusters)] = merged_cluster_B                        
        self.adata.obs['cell_type_3'] = self.adata.obs['cell_type_3'].astype('category')

        tl.rank_genes_groups(self.adata, 'cell_type_3', groups=[merged_cluster_A], reference=merged_cluster_B, method='wilcoxon')            
        pl.rank_genes_groups(self.adata, n_genes=25, sharey=False)

    def get_selected_items(self, list_widget):
        # Get selected items from the list
        return [item.text() for item in list_widget.findItems("*", Qt.MatchWildcard) if item.checkState() == Qt.Checked]

    def save_deg_result(self):
        # Save DEG analysis results to CSV
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
                QMessageBox.information(self, "Success", f"Data saved: {filePath}")
            self.raise_()
            self.activateWindow()
        else:
            QMessageBox.critical(self, "Error", "No results to save. Please run the analysis first.")

class ImageMergeDialog(QDialog):
    # Dialog for merging RGB images
    def __init__(self, parent=None):
        super().__init__(parent)
        self.save_directory = None  # Initialize save directory
        self.initUI()

    def initUI(self):
        # Initialize UI for image merge dialog
        self.setWindowTitle('RGB Image')
        self.setGeometry(50, 50, 700, 300)

        sub_layout = QVBoxLayout()

        self.red_group = self.create_group_box("Red")
        self.green_group = self.create_group_box("Green")
        self.blue_group = self.create_group_box("Blue")
        self.blue_group.setCheckable(True)
        self.blue_group.setChecked(False)

        self.blue_group.toggled.connect(self.toggle_blue)

        sub_layout.addWidget(self.red_group)
        sub_layout.addWidget(self.green_group)
        sub_layout.addWidget(self.blue_group)

        # Save folder group box
        self.save_group_box = self.create_save_group_box("Save Folder")

        # Buttons
        clear_button = QPushButton("Clear", self)
        clear_button.clicked.connect(self.clear_file_paths)

        start_button = QPushButton("Start", self)
        start_button.clicked.connect(self.process_images)

        sub_layout.addWidget(self.save_group_box)  # Add the save folder group box
        sub_layout.addWidget(clear_button)
        sub_layout.addWidget(start_button)

        self.setLayout(sub_layout)

    def create_group_box(self, title):
        # Create a group box for selecting image channels
        group_box = QGroupBox(title)

        h_layout = QHBoxLayout()

        text_browser = QTextBrowser(self)
        text_browser.setPlaceholderText(f"Please select {title} file path")
        text_browser.setFixedHeight(30)

        button = QPushButton("...", self)
        button.clicked.connect(lambda: self.show_file_dialog(text_browser))

        h_layout.addWidget(text_browser)
        h_layout.addWidget(button)

        group_box.setLayout(h_layout)

        group_box.text_browser = text_browser
        return group_box

    def create_save_group_box(self, title):
        # Create a group box for selecting save folder and displaying path
        group_box = QGroupBox(title)

        h_layout = QHBoxLayout()

        # Text browser to display selected save folder path
        self.save_folder_browser = QTextBrowser(self)
        self.save_folder_browser.setPlaceholderText("Selected save folder will be displayed here")
        self.save_folder_browser.setFixedHeight(30)

        # Button to open file dialog to select folder
        select_folder_button = QPushButton("Select Save Folder", self)
        select_folder_button.clicked.connect(self.select_save_directory)

        h_layout.addWidget(self.save_folder_browser)
        h_layout.addWidget(select_folder_button)

        group_box.setLayout(h_layout)
        return group_box

    def show_file_dialog(self, text_browser):
        # Show file dialog to select image files
        options = QFileDialog.Options()
        file_path, _ = QFileDialog.getOpenFileName(self, "Select PNG File", "", "PNG Files (*.png)", options=options)

        if file_path:
            text_browser.setText(file_path)

    def select_save_directory(self):
        # Open a dialog to select the directory where images will be saved
        self.save_directory = QFileDialog.getExistingDirectory(self, "Select Save Folder")
        if self.save_directory:
            # Display the selected folder path in the text browser
            self.save_folder_browser.setText(self.save_directory)

    def toggle_blue(self, checked):
        # Toggle the blue channel selection
        if not checked:
            self.blue_group.text_browser.clear()

    def clear_file_paths(self):
        # Clear file paths from the text browsers
        self.red_group.text_browser.clear()
        self.green_group.text_browser.clear()
        self.blue_group.text_browser.clear()

    def process_images(self):
        # Process and merge selected images
        if not self.save_directory:
            QMessageBox.warning(self, "Warning", "Please select a folder to save the images.")
            return

        red_file = self.red_group.text_browser.toPlainText()
        green_file = self.green_group.text_browser.toPlainText()
        blue_file = self.blue_group.text_browser.toPlainText() if self.blue_group.isChecked() else None

        if not red_file or not green_file:
            QMessageBox.warning(self, "Warning", "Please select Red and Green files.")
            return

        if blue_file:
            self.merge_images([red_file, green_file, blue_file])
        else:
            self.merge_images([red_file, green_file])

    def process_images_common(self, image_path):
        # Convert image to grayscale and invert
        fixed_size = (500,500)
        image = Image.open(image_path)
        img_resized = image.resize(fixed_size, Image.Resampling.LANCZOS)
        bw_image = img_resized.convert('L')
        inverted_image = ImageOps.invert(bw_image)
        return inverted_image.convert('RGB')

    def channel_image(self, inverted_image, color):
        # Extract the selected color channel from the image
        red_channel = green_channel = blue_channel = Image.new('L', inverted_image.size, 0)

        if color == 'red':
            red_channel = inverted_image.split()[0]
        elif color == 'green':
            green_channel = inverted_image.split()[1]
        elif color == 'blue':
            blue_channel = inverted_image.split()[2]
        
        channel_image = Image.merge('RGB', (red_channel, green_channel, blue_channel))

        new_size = (int(channel_image.width * 10), int(channel_image.height * 10))
        high_res_image = channel_image.resize(new_size, Image.Resampling.LANCZOS)

        if self.save_directory:
            high_res_image.save(f'{self.save_directory}/{color}_image.jpg', quality=100)

    def merge_images(self, image_paths):
        # Merge red, green, and optionally blue images into an RGB image
        red_image = self.process_images_common(image_paths[0])
        green_image = self.process_images_common(image_paths[1])

        if len(image_paths) == 3:
            blue_image = self.process_images_common(image_paths[2])
            channels = [red_image.split()[0], green_image.split()[1], blue_image.split()[2]]
            self.channel_image(blue_image, 'blue')
        else:
            channels = [red_image.split()[0], green_image.split()[1], Image.new('L', red_image.size, 0)]

        merged_image = Image.merge('RGB', channels)

        self.channel_image(red_image, 'red')
        self.channel_image(green_image, 'green')

        new_size = (int(merged_image.width * 10), int(merged_image.height * 10))
        high_res_image = merged_image.resize(new_size, Image.Resampling.LANCZOS)
        
        if self.save_directory:
            high_res_image.save(f'{self.save_directory}/merged_image.png', quality=100)

        QMessageBox.information(self, "Success", f"Images saved in {self.save_directory}")


class OtherSamplesDialog(QDialog):
    # Dialog for comparing DEG across multiple samples
    def __init__(self, adata, selected_metadata, parent=None):
        super().__init__(parent)
        self.adata = adata
        self.selected_metadata = selected_metadata
        self.sample_files = [None]
        self.sample_data = {0: adata}
        self.initUI()

    def initUI(self):
        # Initialize UI for comparing samples dialog
        self.setWindowTitle('DEG Other Samples')
        self.setGeometry(50, 50, 1000, 800)

        self.main_layout = QVBoxLayout(self)

        self.label = QLabel("Select the number of Samples:")
        self.combo_box = QComboBox(self)
        self.combo_box.addItems([str(i) for i in range(2, 7)])
        self.select_button = QPushButton("Select", self)
        self.select_button.clicked.connect(self.create_sample_inputs)

        combo_layout = QHBoxLayout()
        combo_layout.addWidget(self.label)
        combo_layout.addWidget(self.combo_box)
        combo_layout.addWidget(self.select_button)

        self.main_layout.addLayout(combo_layout)

        self.sample1_label = QLabel(f"Sample1: {self.parent().current_file_path}")
        self.main_layout.addWidget(self.sample1_label)

        self.sample_input_layout = QVBoxLayout()
        self.main_layout.addLayout(self.sample_input_layout)

        self.show_clusters_button = QPushButton("Show the Clusters list of Samples", self)
        self.show_clusters_button.clicked.connect(self.show_clusters_for_all_samples)
        self.main_layout.addWidget(self.show_clusters_button)

        self.sample_groupbox_layout = QHBoxLayout()
        self.sample1_groupbox = self.create_sample_groupbox(self.adata)
        self.sample_groupbox_layout.addWidget(self.sample1_groupbox)
        self.main_layout.addLayout(self.sample_groupbox_layout)

        self.clear_button = QPushButton("Clear", self)
        self.clear_button.clicked.connect(self.clear_all_selection)
        self.main_layout.addWidget(self.clear_button)

        self.start_button = QPushButton("Start", self)
        self.start_button.clicked.connect(self.compare_gene_expression)
        self.main_layout.addWidget(self.start_button)

        self.save_button = QPushButton("Save", self)
        self.save_button.clicked.connect(self.save_file)
        self.main_layout.addWidget(self.save_button)

    def create_sample_inputs(self):
        # Create input fields for additional samples
        num_samples = int(self.combo_box.currentText())

        for i in reversed(range(self.sample_input_layout.count())):
            item = self.sample_input_layout.itemAt(i)
            if item is not None and item.layout() is not None:
                for j in reversed(range(item.layout().count())):
                    sub_item = item.layout().itemAt(j)
                    if sub_item is not None and sub_item.widget() is not None:
                        sub_item.widget().deleteLater()
                item.layout().deleteLater()

        for i in reversed(range(self.sample_groupbox_layout.count())):
            if i == 0:
                continue
            item = self.sample_groupbox_layout.itemAt(i)
            if item is not None and item.widget() is not None:
                item.widget().deleteLater()

        self.sample_files = [None] * num_samples

        for sample_idx in range(2, num_samples + 1):
            file_label = QLabel(f"Sample{sample_idx}:")
            file_browser = QTextBrowser(self)
            file_browser.setFixedHeight(20)
            file_button = QPushButton("File", self)
            file_button.clicked.connect(lambda _, idx=sample_idx-1, browser=file_browser: self.select_file(idx, browser))

            file_layout = QHBoxLayout()
            file_layout.addWidget(file_label)
            file_layout.addWidget(file_browser)
            file_layout.addWidget(file_button)

            self.sample_input_layout.addLayout(file_layout)

            groupbox = QGroupBox(f"Sample{sample_idx} Clusters")
            layout = QVBoxLayout(groupbox)
            list_widget = QListWidget()
            layout.addWidget(list_widget)
            self.sample_groupbox_layout.addWidget(groupbox)

    def select_file(self, idx, browser):
        # Select a file for each additional sample
        file_path, _ = QFileDialog.getOpenFileName(self, "Select HDF5 File", "", "HDF5 Files (*.h5ad *.hdf5)")
        if file_path:
            browser.setText(file_path)
            self.sample_files[idx] = file_path

    def show_clusters_for_all_samples(self):
        # Load and display clusters for all selected samples
        from scanpy import read_h5ad
        for idx, file_path in enumerate(self.sample_files):
            if idx == 0:
                adata = self.adata
            else:
                try:
                    adata = read_h5ad(file_path)
                    self.sample_data[idx] = adata
                except Exception as e:
                    QMessageBox.critical(self, "Error", f"Failed to load file: {str(e)}")
                    continue

            groupbox_layout = self.sample_groupbox_layout.itemAt(idx).widget().layout()
            list_widget = groupbox_layout.itemAt(0).widget()
            self.update_sample_groupbox(list_widget, adata)

    def create_sample_groupbox(self, adata):
        # Create a group box to display clusters for a sample
        group_box = QGroupBox("Sample1 Clusters")
        layout = QVBoxLayout(group_box)

        cluster_list_widget = QListWidget()
        clusters = sorted(adata.obs[self.selected_metadata].unique())
        for cluster in clusters:
            item = QListWidgetItem(cluster)
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            item.setCheckState(Qt.Unchecked)
            cluster_list_widget.addItem(item)

        layout.addWidget(cluster_list_widget)
        return group_box

    def update_sample_groupbox(self, list_widget, adata):
        # Update the clusters list in the group box for each sample
        list_widget.clear()
        clusters = sorted(adata.obs[self.selected_metadata].unique())
        for cluster in clusters:
            item = QListWidgetItem(cluster)
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            item.setCheckState(Qt.Unchecked)
            list_widget.addItem(item)

    def clear_all_selection(self):
        # Clear all selections in cluster lists
        for i in range(self.sample_groupbox_layout.count()):
            groupbox_layout = self.sample_groupbox_layout.itemAt(i).widget().layout()
            list_widget = groupbox_layout.itemAt(0).widget()
            for j in range(list_widget.count()):
                list_widget.item(j).setCheckState(Qt.Unchecked)

    def compare_gene_expression(self):
        # Compare gene expression between selected samples
        import anndata
        from scanpy import tools as tl
        from scanpy import plotting as pl
        combined_adata_list = []
        sample_names = []

        for idx, adata in self.sample_data.items():
            groupbox_layout = self.sample_groupbox_layout.itemAt(idx).widget().layout()
            list_widget = groupbox_layout.itemAt(0).widget()

            selected_clusters = []
            for j in range(list_widget.count()):
                item = list_widget.item(j)
                if item.checkState() == Qt.Checked:
                    selected_clusters.append(item.text())

            if not selected_clusters:
                continue

            adata_sub = adata[adata.obs[self.selected_metadata].isin(selected_clusters)].copy()

            adata_sub.obs_names = [f"sample{idx+1}_{i}" for i in adata_sub.obs_names]
            combined_adata_list.append(adata_sub)
            sample_names.append(f"sample{idx+1}")

        if not combined_adata_list:
            QMessageBox.warning(self, "Warning", "No clusters selected for comparison.")
            return

        self.combined_adata = anndata.concat(combined_adata_list, label="sample", keys=sample_names, join="outer")

        tl.rank_genes_groups(self.combined_adata, groupby='sample', method='wilcoxon')
        pl.rank_genes_groups(self.combined_adata, n_genes=25, sharey=False)
    
    def save_file(self):
        # Save the DEG comparison results to CSV
        from pandas import DataFrame
        filePath, _ = QFileDialog.getSaveFileName(self, "Save CSV File", "", "CSV Files (*.csv);;All Files (*)")
        
        if not filePath:
            return

        if 'rank_genes_groups' in self.combined_adata.uns:
            result = self.combined_adata.uns['rank_genes_groups']
            groups = result['names'].dtype.names
            result_df = DataFrame(
                {group + '_' + key: result[key][group]
                for group in groups for key in ['names', 'scores', 'logfoldchanges', 'pvals', 'pvals_adj']}
            )
            if not filePath.endswith('.csv'):
                filePath += '.csv'
            result_df.to_csv(filePath, index=False)
            QMessageBox.information(self, "Success", f"Data saved: {filePath}")

            self.raise_()
            self.activateWindow()
        else:
            QMessageBox.critical(self, "Error", "No results to save. Please run the analysis first.")

        
if __name__ == '__main__':
    from multiprocessing import freeze_support
    freeze_support()

    app = QApplication(sys.argv)
    main_window = MainWindow()
    main_window.show()
    sys.exit(app.exec_())
