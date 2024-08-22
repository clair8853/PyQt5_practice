from PyQt5.QtWidgets import QFileDialog, QMessageBox
from PyQt5.QtCore import Qt
import scanpy as sc
import squidpy as sq
from matplotlib.transforms import Affine2D
import numpy as np
import pandas as pd


def handle_load_data(app):
    if app.adata is not None:
        reply = QMessageBox.question(app, 'Message',
                                     "A file is already loaded. Do you want to load a new file and discard the current data?",
                                     QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.No:
            return
    
    options = QFileDialog.Options()
    file_name, _ = QFileDialog.getOpenFileName(app, "Open Anndata File", "", "Anndata Files (*.hdf5);;All Files (*)", options=options)
    if file_name:
        try:
            app.adata = sc.read_h5ad(file_name)
            app.statusBar().showMessage(f'Loaded {file_name}')
            app.update_data_info()
            app.update_lists()
            app.plot_umap()
        except Exception as e:
            QMessageBox.critical(app, "Error", f"Could not load file:\n{str(e)}")

def handle_sync_lists(app, item):
    app.cluster1_list.blockSignals(True)
    app.cluster2_list.blockSignals(True)

    if item.listWidget() == app.cluster1_list:
        other_list = app.cluster2_list
    else:
        other_list = app.cluster1_list

    for i in range(other_list.count()):
        other_item = other_list.item(i)
        if other_item.text() == item.text():
            other_item.setFlags(other_item.flags() & ~Qt.ItemIsEnabled if item.checkState() == Qt.Checked else other_item.flags() | Qt.ItemIsEnabled)

    app.cluster1_list.blockSignals(False)
    app.cluster2_list.blockSignals(False)

def handle_plot_umap(app):
    if app.adata is not None:
        app.plot_canvas.figure.clear()
        ax = app.plot_canvas.figure.add_subplot(111)
        # Extract x and y coordinates from adata
        x_coords = app.adata.obsm['spatial'][:, 0]
        y_coords = app.adata.obsm['spatial'][:, 1]

        # Calculate the center of the plot
        center_x = (x_coords.max() + x_coords.min()) / 2
        center_y = (y_coords.max() + y_coords.min()) / 2

        # Apply rotation transformation around the center
        rotation = Affine2D().rotate_around(center_x, center_y, 0)
        ax.transData = rotation + ax.transData

        # Plot using squidpy
        sq.pl.spatial_scatter(app.adata, shape=None, color="cell_type_2", ax=ax, legend_loc=None,frameon=False, title=None)

        # Set the title
        ax.set_title('')

        # Calculate new limits after rotation
        coords = rotation.transform(np.vstack([x_coords, y_coords]).T)
        x_min, y_min = coords.min(axis=0)
        x_max, y_max = coords.max(axis=0)

        # Set new limits
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)

        app.plot_canvas.draw()

def handle_start_plotting(app):
    if app.adata is not None:
        selected_genes = [item.text() for item in app.gene_list.findItems("*", Qt.MatchWildcard) if item.checkState() == Qt.Checked]
        
        if not selected_genes:
            QMessageBox.warning(app, "Warning", "Please select at least one gene.")
            return

        plot_type = None
        if app.dot_plot_radio.isChecked():
            plot_type = 'dotplot'
        elif app.violin_plot_radio.isChecked():
            plot_type = 'violin'
        elif app.feature_plot_radio.isChecked():
            plot_type = 'feature'
        elif app.spatial_plot_radio.isChecked():
            plot_type = 'spatial'

        if plot_type is not None:
            app.plot(selected_genes, plot_type)

def handle_plotting_clusters(app):
    selected_clusters1 = [item.text() for item in app.cluster1_list.findItems("*", Qt.MatchWildcard) if item.checkState() == Qt.Checked]
                
    app.plot_canvas.figure.clear()
    plot_type = None
    if app.radio_all.isChecked():
        plot_type = 'all'
    elif app.radio_select.isChecked():
        plot_type = 'select'
        
    if plot_type is not None:
        app.plot_clusters(selected_clusters1, plot_type)  

def handle_spatial_plot(app):
    selected_clusters1 = [item.text() for item in app.cluster1_list.findItems("*", Qt.MatchWildcard) if item.checkState() == Qt.Checked]
    if not selected_clusters1:
        QMessageBox.warning(app, "Warning", "Please select at least one cluster in Cluster 1.")
        return
    else:            
        app.plot_canvas.figure.clear()
        
        ax = app.plot_canvas.figure.add_subplot(111)

        sq.pl.spatial_scatter(app.adata, shape=None, color="cell_type_2", groups=selected_clusters1, ax=ax, frameon=False, title=None, legend_fontsize=15)

        x_coords = app.adata.obsm['spatial'][:, 0]
        y_coords = app.adata.obsm['spatial'][:, 1]
        ax.set_xlim(x_coords.min(), x_coords.max())
        ax.set_ylim(y_coords.min(), y_coords.max())
        ax.set_title('')
        app.plot_canvas.draw()

def handle_deg_plot(app):
    app.save_deg_button.setEnabled(True)
    selected_clusters1 = [item.text() for item in app.cluster1_list.findItems("*", Qt.MatchWildcard) if item.checkState() == Qt.Checked]
    selected_clusters2 = [item.text() for item in app.cluster2_list.findItems("*", Qt.MatchWildcard) if item.checkState() == Qt.Checked]
    if not selected_clusters1 or not selected_clusters2:
        QMessageBox.warning(app, "Warning", "Please select at least one cluster in both Cluster 1 and Cluster 2.")
        return
    else:
        merged_cluster_A = '_'.join(selected_clusters1)
        merged_cluster_B = '_'.join(selected_clusters2)
                    
        app.adata.obs['cell_type_3'] = app.adata.obs['cell_type_2']
        app.adata.obs['cell_type_3'] = app.adata.obs['cell_type_3'].astype(str)
        app.adata.obs['cell_type_3'][app.adata.obs['cell_type_3'].isin(selected_clusters1)] = merged_cluster_A
        app.adata.obs['cell_type_3'][app.adata.obs['cell_type_3'].isin(selected_clusters2)] = merged_cluster_B                        
        app.adata.obs['cell_type_3'] = app.adata.obs['cell_type_3'].astype('category')
        
        sc.tl.rank_genes_groups(app.adata, 'cell_type_3', groups=[merged_cluster_A], reference=merged_cluster_B, method='wilcoxon')            
        sc.pl.rank_genes_groups(app.adata, n_genes=25, sharey=False)

def handle_save_file(app):
    options = QFileDialog.Options()
    options |= QFileDialog.DontUseNativeDialog
    filePath, _ = QFileDialog.getSaveFileName(app, "Save CSV File", "", "CSV Files (*.csv);;All Files (*)", options=options)
    if 'rank_genes_groups' in app.adata.uns:            
        result = app.adata.uns['rank_genes_groups']
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
            app.statusBar().showMessage(info)                     
    else:
        QMessageBox.critical(app, "Error", f"No results to save. Please run the analysis first.")
