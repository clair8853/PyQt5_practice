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
import tempfile
import gc
import plotly.graph_objects as go


class DataProcessingThread(QThread):
    progress_changed = pyqtSignal(int)
    processing_finished = pyqtSignal(object)

    def __init__(self, file_name, parent=None):
        super().__init__(parent)
        self.file_name = file_name
        self.adata = None

    def run(self):
        from scanpy import read_h5ad
        self.progress_changed.emit(10)
        try:
            self.adata = read_h5ad(self.file_name)
            self.progress_changed.emit(50)
        except Exception as e:
            self.progress_changed.emit(0)
            self.processing_finished.emit(e)
            return

        self.progress_changed.emit(80)
        self.progress_changed.emit(100)
        self.processing_finished.emit(self.adata)


class MainWindow(QMainWindow):

    def __init__(self):
        super().__init__()
        self.adata = None
        self.progress_bar = None
        self.selected_metadata = 'cell_type_2'
        self.font_settings = {"family": "Arial", "size": 15}
        self.initUI()

    def initUI(self):
        self.setup_menu_bar()
        self.setup_status_bar()
        self.setup_main_layout()
        self.setWindowTitle('PyQt5 Application')
        self.setGeometry(100, 100, 1200, 800)

    def setup_menu_bar(self):
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
        self.statusbar = self.statusBar()
        self.statusbar.showMessage('Ready')
        self.progress_bar = QProgressBar(self)
        self.statusbar.addPermanentWidget(self.progress_bar)
        self.progress_bar.setValue(0)

    def setup_main_layout(self):
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)

        main_splitter = QSplitter(Qt.Horizontal)
        self.central_widget.setLayout(QVBoxLayout())
        self.central_widget.layout().addWidget(main_splitter)

        self.setup_left_area(main_splitter)
        self.setup_central_area(main_splitter)
        
        main_splitter.setSizes([100, 1100])

    def setup_left_area(self, parent_splitter):
        left_panel = QGroupBox()
        left_layout = QGridLayout(left_panel)

        self.create_group_box('Spatial Plot', left_layout, 0, 0, 1, 2, self.create_spatial_plot_controls)
        self.create_group_box('Genes', left_layout, 1, 0, 1, 2, self.create_genes_list)
        self.create_group_box('Plotting', left_layout, 2, 0, 1, 2, self.create_plotting_controls)
        parent_splitter.addWidget(left_panel)

    def setup_central_area(self, parent_splitter):
        central_panel = QWidget()
        self.central_layout = QVBoxLayout(central_panel)
        self.web_view = QWebEngineView(self)
        self.central_layout.addWidget(self.web_view)
        parent_splitter.addWidget(central_panel)

    def create_group_box(self, title, layout, row, col, rowspan, colspan, setup_func):
        group_box = QGroupBox(title)
        group_layout = QGridLayout(group_box)
        setup_func(group_layout)
        layout.addWidget(group_box, row, col, rowspan, colspan)

    def create_spatial_plot_controls(self, layout):
        axis_label = QLabel('Axis :')
        self.axis_input = QLineEdit('0')
        layout.addWidget(axis_label, 0, 0)
        layout.addWidget(self.axis_input, 0, 1)

        start_button = QPushButton('Start')
        start_button.clicked.connect(self.plot_spatial_scatter)
        layout.addWidget(start_button, 1, 0, 1, 2)

    def create_genes_list(self, layout):
        self.genes_list_widget = QListWidget()
        layout.addWidget(self.genes_list_widget, 0, 0, 1, 2)

        clear_button = QPushButton('Clear')
        clear_button.clicked.connect(self.clear_gene_selection)
        layout.addWidget(clear_button, 1, 0, 1, 2)

    def create_plotting_controls(self, layout):
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
        self.deg_action.setEnabled(True)
        self.progress_bar.setValue(100)
        self.statusbar.showMessage("Data loaded successfully.")
        self.setup_metadata_selector()

    def setup_metadata_selector(self):
        if self.adata is None:
            return

        metadata_layout = self.central_widget.layout()
        self.metadata_selector = QComboBox(self)
        self.metadata_selector.addItems(self.adata.obs.columns)
        self.metadata_selector.currentTextChanged.connect(self.on_metadata_selected)
        metadata_layout.insertWidget(0, self.metadata_selector)

    def on_metadata_selected(self, selected_metadata):
        self.selected_metadata = selected_metadata
        QMessageBox.information(self, "Metadata Selected", f"You selected {self.selected_metadata}")

    def clear_gene_selection(self):
        for i in range(self.genes_list_widget.count()):
            item = self.genes_list_widget.item(i)
            item.setCheckState(Qt.Unchecked)

    def plot_spatial_scatter(self):
        if self.adata is None:
            return

        axis = self.get_rotation_axis()

        try:
            color_categories = self.adata.obs[self.selected_metadata].astype('category').cat.categories
            color_map = {category: color for category, color in zip(color_categories, self.adata.uns[self.selected_metadata + '_colors'])}

            categories_sorted = sorted(self.adata.obs[self.selected_metadata].unique(), key=str)

            rotated_x_coords, rotated_y_coords = self.rotate_coords(axis)
            fig = go.Figure()

            for category in categories_sorted:
                mask = self.adata.obs[self.selected_metadata] == category
                fig.add_trace(go.Scattergl(
                    x=rotated_x_coords[mask], y=rotated_y_coords[mask],
                    mode='markers', name=str(category),
                    marker=dict(size=2, color=color_map[category])
                ))

            fig.update_layout(
                xaxis=dict(scaleanchor="y", scaleratio=1, range=[rotated_x_coords.min(), rotated_x_coords.max()], showticklabels=False),
                yaxis=dict(scaleanchor="x", scaleratio=1, range=[rotated_y_coords.min(), rotated_y_coords.max()], showticklabels=False),
                legend=dict(font=dict(family=self.font_settings["family"], size=self.font_settings["size"]), itemsizing="constant"),
                plot_bgcolor="white", paper_bgcolor="white"
            )

            self.display_plot(fig)

        except Exception as e:
            QMessageBox.warning(self, 'Error', f'Could not plot data:\n{str(e)}')
        finally:
            gc.collect()

    def get_rotation_axis(self):
        num = float(self.axis_input.text())
        return 180 / num if num > 0 else 360

    def rotate_coords(self, axis):
        x_coords = self.adata.obsm['spatial'][:, 0]
        y_coords = self.adata.obsm['spatial'][:, 1]
        center_x = (x_coords.max() + x_coords.min()) / 2
        center_y = (y_coords.max() + y_coords.min()) / 2

        rotation_matrix = array([[cos(pi / axis), -sin(pi / axis)], [sin(pi / axis), cos(pi / axis)]])
        rotated_coords = dot(vstack([x_coords - center_x, y_coords - center_y]).T, rotation_matrix)
        return rotated_coords[:, 0] + center_x, rotated_coords[:, 1] + center_y

    def display_plot(self, fig):
        with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as temp_file:
            fig.write_html(temp_file.name)
            self.web_view.setUrl(QUrl.fromLocalFile(temp_file.name))
