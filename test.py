import scanpy as sc

adata = sc.read_h5ad(r"C:\Users\visiu\Documents\pyqt5_proj\P4_visual.hdf5")

import plotly.express as px
import numpy as np

# 좌표 데이터 가져오기
x_coords = adata.obsm['spatial'][:, 0]
y_coords = adata.obsm['spatial'][:, 1]

# 중심 좌표 계산
center_x = (x_coords.max() + x_coords.min()) / 2
center_y = (y_coords.max() + y_coords.min()) / 2

# 회전 행렬을 적용하여 좌표 회전
rotation_matrix = np.array([
    [np.cos(np.pi / 2), -np.sin(np.pi / 2)],
    [np.sin(np.pi / 2), np.cos(np.pi / 2)]
])
rotated_coords = np.dot(np.vstack([x_coords - center_x, y_coords - center_y]).T, rotation_matrix)
rotated_x_coords = rotated_coords[:, 0] + center_x
rotated_y_coords = rotated_coords[:, 1] + center_y

# Plotly를 사용하여 플롯 생성
fig = px.scatter(x=rotated_x_coords,
                 y=rotated_y_coords,
                 color=adata.obs['leiden'],
                 labels={'x': 'Rotated X', 'y': 'Rotated Y'},
                 title='Rotated Spatial Plot')

fig.update_traces(marker=dict(size=2))

fig.update_layout(
    xaxis=dict(
        scaleanchor="y",  
        scaleratio=1      
    ),
    yaxis=dict(
        scaleanchor="x",  
        scaleratio=1      
    ),
    plot_bgcolor="white",  
    paper_bgcolor="white"  
)

# 그래프 표시
fig.show()
