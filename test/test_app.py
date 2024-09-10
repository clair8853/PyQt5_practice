import sys
import numpy as np
import plotly.express as px
from PyQt5.QtWidgets import QApplication, QMainWindow, QAction, QFileDialog, QMessageBox, QStatusBar, QPushButton, QVBoxLayout, QWidget
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtCore import QUrl
import anndata as ad
import tempfile

class AnnDataApp(QMainWindow):
    def __init__(self):
        super().__init__()

        # anndata 객체를 저장할 변수
        self.adata = None

        # 윈도우 설정
        self.setWindowTitle('AnnData Loader and Plotter')
        self.setGeometry(100, 100, 800, 600)

        # 메뉴바 설정
        menubar = self.menuBar()
        file_menu = menubar.addMenu('File')

        # 'Open File' 액션 추가
        open_file_action = QAction('Open File', self)
        open_file_action.triggered.connect(self.open_file)
        file_menu.addAction(open_file_action)

        # 스테이터스바 설정
        self.statusbar = QStatusBar()
        self.setStatusBar(self.statusbar)

        # 플롯 버튼 설정
        self.plot_button = QPushButton('Plot Data', self)
        self.plot_button.clicked.connect(self.plot_data)
        self.plot_button.setEnabled(False)  # 데이터가 로드될 때까지 비활성화

        # QWebEngineView 설정
        self.web_view = QWebEngineView(self)

        # 중앙 위젯 및 레이아웃 설정
        central_widget = QWidget()
        layout = QVBoxLayout()
        layout.addWidget(self.plot_button)
        layout.addWidget(self.web_view)
        central_widget.setLayout(layout)
        self.setCentralWidget(central_widget)

    def open_file(self):
        # 파일 열기 다이얼로그 표시
        file_name, _ = QFileDialog.getOpenFileName(self, 'Open File', '', 'AnnData Files (*.h5ad *.hdf5);;All Files (*)')

        if file_name:
            try:
                # anndata 파일 로드
                self.adata = ad.read_h5ad(file_name)  # 또는 ad.read_hdf(file_name) for .hdf5
                self.statusbar.showMessage(f'Loaded {file_name}')
                self.plot_button.setEnabled(True)  # 파일이 로드되면 플롯 버튼 활성화
            except Exception as e:
                # 파일 로드 실패 시 경고 메시지 표시
                QMessageBox.warning(self, 'Error', f'Could not open file: {file_name}\n{str(e)}')

    def plot_data(self):
        if self.adata is not None:
            try:
                # Plotting code using anndata data
                x_coords = self.adata.obsm['spatial'][:, 0]
                y_coords = self.adata.obsm['spatial'][:, 1]

                center_x = (x_coords.max() + x_coords.min()) / 2
                center_y = (y_coords.max() + y_coords.min()) / 2

                rotation_matrix = np.array([
                    [np.cos(np.pi / 2), -np.sin(np.pi / 2)],
                    [np.sin(np.pi / 2), np.cos(np.pi / 2)]
                ])
                rotated_coords = np.dot(np.vstack([x_coords - center_x, y_coords - center_y]).T, rotation_matrix)
                rotated_x_coords = rotated_coords[:, 0] + center_x
                rotated_y_coords = rotated_coords[:, 1] + center_y

                # 전체 데이터에 대한 x, y 축의 범위 설정
                x_range = [rotated_x_coords.min(), rotated_x_coords.max()]
                y_range = [rotated_y_coords.min(), rotated_y_coords.max()]

                fig = px.scatter(x=rotated_x_coords,
                                 y=rotated_y_coords,
                                 color=self.adata.obs['leiden'])

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

def main():
    app = QApplication(sys.argv)
    window = AnnDataApp()
    window.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
