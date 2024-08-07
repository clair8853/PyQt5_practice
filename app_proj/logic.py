from PyQt5.QtWidgets import QFileDialog, QMessageBox
from PyQt5.QtCore import Qt
import scanpy as sc
import squidpy as sq
import matplotlib.colors as mcolors
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
        sq.pl.spatial_scatter(app.adata, shape=None, color=['cell_type_2'], legend_loc=None, ax=ax, frameon=False, title=None, size=5)
        ax.set_title('')
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

        if plot_type is not None:
            app.plot(selected_genes, plot_type)

def handle_spatial_plot(app):
    selected_clusters1 = [item.text() for item in app.cluster1_list.findItems("*", Qt.MatchWildcard) if item.checkState() == Qt.Checked]
    if not selected_clusters1:
        QMessageBox.warning(app, "Warning", "Please select at least one cluster in Cluster 1.")
        return
    else:            
        app.plot_canvas.figure.clear()
        ax = app.plot_canvas.figure.add_subplot(111)
        app.adata.obs['cell_type_spatial'] = app.adata.obs['cell_type_2']
        app.adata.obs['cell_type_spatial'] = app.adata.obs['cell_type_spatial'].astype(str)
        app.adata.obs.loc[~app.adata.obs['cell_type_spatial'].isin(selected_clusters1), 'cell_type_spatial'] = 'other cell type'
        app.adata.obs['cell_type_spatial'] = app.adata.obs['cell_type_spatial'].astype('category')
        categories = app.adata.obs['cell_type_2'].cat.categories
        category_colors = dict(zip(categories, app.adata.uns['cell_type_2_colors']))
        new_categories = app.adata.obs['cell_type_spatial'].cat.categories
        new_category_colors = {cat: category_colors[cat] if cat in category_colors else '#f8f9fa' for cat in new_categories}
        custom_cmap = mcolors.ListedColormap([new_category_colors[key] for key in new_category_colors])

        sq.pl.spatial_scatter(app.adata, shape=None, color="cell_type_spatial", palette=custom_cmap, ax=ax, frameon=False, title=None, size=5)
        
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
