import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QHBoxLayout, QComboBox, QGroupBox, QPushButton, QLabel, QWidget


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        # 메인 레이아웃 설정
        self.central_widget = QWidget(self)
        self.setCentralWidget(self.central_widget)
        self.main_layout = QVBoxLayout(self.central_widget)

        # ComboBox와 버튼을 포함할 레이아웃 설정 (가로로 배치)
        self.control_layout = QHBoxLayout()
        self.main_layout.addLayout(self.control_layout)

        # ComboBox 생성 및 초기화
        self.combo_box = QComboBox(self)
        self.combo_box.addItems([str(i) for i in range(1, 6)])  # 1부터 5까지 숫자 추가
        self.control_layout.addWidget(self.combo_box)

        # Select 버튼 생성
        self.select_button = QPushButton("Select", self)
        self.select_button.clicked.connect(self.on_button_clicked)
        self.control_layout.addWidget(self.select_button)

        # GroupBox를 추가할 레이아웃 (가로로 배치)
        self.groupbox_layout = QHBoxLayout()
        self.main_layout.addLayout(self.groupbox_layout)

        # 항상 존재하는 GroupBox1 생성
        self.group_box1 = QGroupBox("GroupBox 1", self)
        self.groupbox_layout.addWidget(self.group_box1)

        group_box1_layout = QVBoxLayout()
        group_box1_layout.addWidget(QLabel("Content 1"))
        self.group_box1.setLayout(group_box1_layout)

    def on_button_clicked(self):
        # ComboBox에서 선택된 숫자에 맞춰 groupbox를 업데이트
        num_groupboxes = int(self.combo_box.currentText())

        # GroupBox1을 제외한 기존의 groupbox들을 모두 제거
        for i in reversed(range(1, self.groupbox_layout.count())):  # 첫 번째 GroupBox는 유지
            widget = self.groupbox_layout.itemAt(i).widget()
            if widget is not None:
                widget.setParent(None)

        # 새로운 groupbox들을 추가
        for i in range(num_groupboxes):
            group_box = QGroupBox(f"GroupBox {i+2}", self)  # GroupBox 번호는 2부터 시작
            group_layout = QVBoxLayout()
            group_layout.addWidget(QLabel(f"Content {i+2}"))
            group_box.setLayout(group_layout)
            self.groupbox_layout.addWidget(group_box)


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
