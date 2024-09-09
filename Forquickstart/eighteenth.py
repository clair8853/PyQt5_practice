import sys
import matplotlib
matplotlib.use('Qt5Agg')
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtWidgets import (QApplication, QMainWindow, QAction, QFileDialog, QHBoxLayout,
                             QMessageBox, QWidget, QVBoxLayout, QGridLayout, QGroupBox, 
                             QRadioButton, QLabel, QLineEdit, QPushButton, QProgressBar,
                             QListWidget, QListWidgetItem, QSplitter, QDialog, QFontDialog,
                             QTextBrowser)
from PyQt5.QtCore import Qt, QUrl, QThread, pyqtSignal
from PIL import Image, ImageOps
from numpy import array, cos, sin, dot, vstack, pi
import tempfile
import gc
import plotly.graph_objects as go


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
        self.font_settings = {"family": "Arial", "size": 15}  # Default font settings for plotly legend
        self.initUI()

    def initUI(self):
        # Initialize the main UI layout and components
        self.setup_menu_bar()
        self.setup_status_bar()
        self.setup_main_layout()
        self.setWindowTitle('PyQt5 Application')
        self.setGeometry(100, 100, 1200, 800)

    def setup_menu_bar(self):
        # Setup the menu bar with 'File', 'Tool', and 'Font' menus
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

        self.image_action = QAction('In-Situ Hybridization', self)
        self.image_action.setShortcut('Ctrl+I')
        self.image_action.triggered.connect(self.open_image_merge)
        tool_menu.addAction(self.image_action)

        # Add Font menu for changing font in Plotly plots and main layout
        font_menu = menubar.addMenu('&Font')
        
        # Font for Plotly
        font_action = QAction('Main App Window', self)
        font_action.triggered.connect(self.change_font_for_plot)
        font_menu.addAction(font_action)

        # Font for Main Layout
        layout_font_action = QAction('Main Layout', self)
        layout_font_action.triggered.connect(self.change_font_for_layout)
        font_menu.addAction(layout_font_action)

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
        
        main_splitter.setSizes([100, 1100])

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
        layout.addWidget(self.violin_plot_radio_button, 1, 0)
        layout.addWidget(self.feature_plot_radio_button, 2, 0)
        layout.addWidget(self.scatter_plot_radio_button, 3, 0)
        layout.addWidget(start_button, 4, 0, 1, 2)
    
    def change_font_for_layout(self):
        # Open QFontDialog to select font
        font, ok = QFontDialog.getFont()
        if ok:
            # Apply the font to the entire layout
            self.setFont(font)
            self.update_layout_font(self.central_widget, font)
    
    def update_layout_font(self, widget, font):
        # Recursively apply the font to all child widgets
        widget.setFont(font)
        for child in widget.findChildren(QWidget):
            child.setFont(font)

    def change_font_for_plot(self):
        # Open QFontDialog to select font
        font, ok = QFontDialog.getFont()
        if ok:
            # Update the font settings for plotly legend
            self.font_settings["family"] = font.family()
            self.font_settings["size"] = font.pointSize()
            self.plot_spatial_scatter()

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
        # Plot the spatial scatter plot using Plotly
        if self.adata is None:
            return
        
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

            # Apply custom font to legend
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
                scatter = ax.scatter(x=rotated_x_coords, y=rotated_y_coords, c=gene_expression, cmap='Purples', s=0.1, edgecolors='none')
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

    def open_image_merge(self):
        dialog = ImageMergeDialog(self)
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

class ImageMergeDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.initUI()

    def initUI(self):
        self.setWindowTitle('RGB Image')
        self.setGeometry(100, 100, 600, 300)

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

        clear_button = QPushButton("Clear", self)
        clear_button.clicked.connect(self.clear_file_paths)

        start_button = QPushButton("Start", self)
        start_button.clicked.connect(self.process_images)

        sub_layout.addWidget(clear_button)
        sub_layout.addWidget(start_button)

        self.setLayout(sub_layout)

    def create_group_box(self, title):
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
    
    def show_file_dialog(self, text_browser):
        options = QFileDialog.Options()
        file_path, _ = QFileDialog.getOpenFileName(self, "Select PNG File", "", "PNG Files (*.png)", options=options)

        if file_path:
            text_browser.setText(file_path)

    def toggle_blue(self, checked):
        if not checked:
            self.blue_group.text_browser.clear()

    def clear_file_paths(self):
        self.red_group.text_browser.clear()
        self.green_group.text_browser.clear()
        self.blue_group.text_browser.clear()

    def process_images(self):
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
        image = Image.open(image_path)
        bw_image = image.convert('L')
        inverted_image = ImageOps.invert(bw_image)
        return inverted_image.convert('RGB')
    
    def channel_image(self, inverted_image, color):
        if color == 'red':
            red_channel = inverted_image.split()[0]
            green_channel = Image.new('L', inverted_image.size, 0)
            blue_channel = Image.new('L', inverted_image.size, 0)
        elif color == 'green':
            green_channel = inverted_image.split()[1]
            red_channel = Image.new('L', inverted_image.size, 0)
            blue_channel = Image.new('L', inverted_image.size, 0)
        elif color == 'blue':
            blue_channel = inverted_image.split()[2]
            green_channel = Image.new('L', inverted_image.size, 0)
            red_channel = Image.new('L', inverted_image.size, 0)
        
        channel_image = Image.merge('RGB', (red_channel, green_channel, blue_channel))

        new_size = (int(channel_image.width * 10), int(channel_image.height * 10))
        high_res_image = channel_image.resize(new_size, Image.Resampling.LANCZOS)

        high_res_image.show()

    def merge_images(self, image_paths):
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
        high_res_image.show()   


if __name__ == '__main__':
    from multiprocessing import freeze_support
    freeze_support()

    app = QApplication(sys.argv)
    main_window = MainWindow()
    main_window.show()
    sys.exit(app.exec_())
