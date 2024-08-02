from PyQt5.QtWidgets import QFileDialog
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
import scanpy as sc

def open_file_dialog():
    options = QFileDialog.Options()
    file_name, _ = QFileDialog.getOpenFileName(None, "QFileDialog.getOpenFileName()", "", "All Files (*);;Python Files (*.py)", options=options)
    return file_name

def plot_umap(adata):
    ncols = 2
    nrows = 1
    figsize = 3  # figsize를 조금 더 크게 설정
    wspace = 0.4
    hspace = 0.4  # 세로 간격을 추가
    fig, axs = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(ncols * figsize, nrows * figsize)
    )
    plt.subplots_adjust(wspace=wspace, hspace=hspace)
    
    sc.pl.umap(adata, color=["cell_type_1"], legend_loc="on data", title='UMAP - Cell Type 1', ax=axs[0], frameon=False, show=False)
    sc.pl.umap(adata, color=["cell_type_2"], title='UMAP - Cell Type 2',  ax=axs[1], frameon=False, show=False)
    
    # 범례 위치 조정
    handles, labels = axs[1].get_legend_handles_labels()
    labels_sorted, handles_sorted = zip(*sorted(zip(labels, handles)))
    axs[1].legend(handles_sorted, labels_sorted, loc='upper left', bbox_to_anchor=(1.05, 1), ncol=1, frameon=False)
    
    canvas = FigureCanvas(fig)
    return canvas
