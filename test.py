import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QFileDialog, QTextBrowser, QVBoxLayout, QHBoxLayout, QGroupBox, QWidget, QCheckBox
from PyQt5.QtCore import Qt
from PIL import Image, ImageOps
import os

class FileDialogDemo(QMainWindow):
    def __init__(self):
        super().__init__()

        self.initUI()

    def initUI(self):
        # 메인 레이아웃
        main_layout = QVBoxLayout()

        # Red, Green, Blue 그룹박스 추가
        self.red_group = self.create_group_box("Red")
        self.green_group = self.create_group_box("Green")
        self.blue_group = self.create_group_box("Blue")
        self.blue_group.setEnabled(False)  # Blue는 처음에 비활성화

        # Blue 옆에 활성화 체크박스 추가
        blue_check_layout = QHBoxLayout()
        self.blue_checkbox = QCheckBox("Enable Blue", self)
        self.blue_checkbox.stateChanged.connect(self.toggle_blue_group)
        blue_check_layout.addWidget(self.blue_checkbox)

        # 그룹박스와 체크박스를 메인 레이아웃에 추가
        main_layout.addWidget(self.red_group)
        main_layout.addWidget(self.green_group)
        main_layout.addLayout(blue_check_layout)
        main_layout.addWidget(self.blue_group)

        # Start 버튼 추가
        self.start_button = QPushButton("Start", self)
        self.start_button.clicked.connect(self.process_images)
        main_layout.addWidget(self.start_button)

        # 중앙 위젯 설정
        container = QWidget()
        container.setLayout(main_layout)
        self.setCentralWidget(container)

        self.setWindowTitle("File Dialog Example")
        self.setGeometry(100, 100, 400, 300)

    def create_group_box(self, title):
        # 그룹박스 생성
        group_box = QGroupBox(title)

        # 그룹박스 내의 가로 레이아웃
        h_layout = QHBoxLayout()

        # 파일 경로를 표시할 QTextBrowser (한 줄 크기로 설정)
        text_browser = QTextBrowser(self)
        text_browser.setPlaceholderText(f"Selected {title} file path will be shown here")
        text_browser.setFixedHeight(30)  # 한 줄 정도의 높이로 고정

        # 파일 선택 버튼
        button = QPushButton("button", self)
        button.setFixedWidth(50)  # 버튼의 너비를 50으로 설정
        button.clicked.connect(lambda: self.show_file_dialog(text_browser))

        # 텍스트 브라우저와 버튼을 가로 레이아웃에 추가
        h_layout.addWidget(text_browser)
        h_layout.addWidget(button)

        # 그룹박스에 가로 레이아웃 설정
        group_box.setLayout(h_layout)

        # 그룹박스에 있는 텍스트 브라우저를 반환
        group_box.text_browser = text_browser
        return group_box

    def show_file_dialog(self, text_browser):
        # 파일 탐색기 열기, PNG 파일만 허용
        options = QFileDialog.Options()
        file_path, _ = QFileDialog.getOpenFileName(self, "Select PNG File", "", "PNG Files (*.png)", options=options)

        if file_path:
            # 선택된 파일 경로를 QTextBrowser에 표시
            text_browser.setText(file_path)

    def toggle_blue_group(self):
        # Blue 그룹박스 활성화/비활성화
        self.blue_group.setEnabled(self.blue_checkbox.isChecked())

    def process_images(self):
        # 활성화된 그룹박스에서 선택된 파일을 가져옴
        red_file = self.red_group.text_browser.toPlainText()
        green_file = self.green_group.text_browser.toPlainText()
        blue_file = self.blue_group.text_browser.toPlainText() if self.blue_group.isEnabled() else None

        # 선택된 파일이 있는지 확인
        if not red_file or not green_file:
            print("Red and Green files must be selected.")
            return

        if blue_file:
            # Red, Green, Blue 이미지 처리 (3색 처리)
            self.merge_images([red_file, green_file, blue_file], "merged_image_3colors.jpg")
        else:
            # Red, Green 이미지 처리 (2색 처리)
            self.merge_images([red_file, green_file], "merged_image_2colors.jpg")

    def process_images_common(self, image_path):
        """공통 이미지 처리: 이미지 로드, 흑백 변환 및 반전"""
        image = Image.open(image_path)
        bw_image = image.convert('L')
        inverted_image = ImageOps.invert(bw_image)
        return inverted_image.convert('RGB')

    def save_channel_image(self, inverted_image, color, output_filename):
        """각 채널별로 이미지를 저장하는 함수"""
        if color == 'red':
            red_channel = inverted_image.split()[0]
            green_channel = Image.new('L', inverted_image.size, 0)
            blue_channel = Image.new('L', inverted_image.size, 0)
        elif color == 'green':
            red_channel = Image.new('L', inverted_image.size, 0)
            green_channel = inverted_image.split()[1]
            blue_channel = Image.new('L', inverted_image.size, 0)
        elif color == 'blue':
            red_channel = Image.new('L', inverted_image.size, 0)
            green_channel = Image.new('L', inverted_image.size, 0)
            blue_channel = inverted_image.split()[2]

        channel_image = Image.merge('RGB', (red_channel, green_channel, blue_channel))

        # 해상도를 10배로 확대
        new_size = (int(channel_image.width * 10), int(channel_image.height * 10))
        high_res_image = channel_image.resize(new_size, Image.Resampling.LANCZOS)

        # 고해상도로 저장
        high_res_image.save(output_filename, quality=100)
        print(f"{output_filename} saved successfully.")

    def merge_images(self, image_paths, output_filename):
        """이미지를 병합하고 고해상도로 저장하는 함수"""
        # 각 이미지를 채널로 변환
        channels = [self.process_images_common(image_path).split()[i] for i, image_path in enumerate(image_paths)]
        
        # Blue 채널이 없으면 빈 채널 추가
        if len(channels) == 2:
            channels.append(Image.new('L', channels[0].size, 0))  # 빈 Blue 채널 추가

        # 새로운 RGB 이미지 생성
        merged_image = Image.merge('RGB', channels)

        # Red, Green, Blue 각각의 이미지를 저장
        self.save_channel_image(self.process_images_common(image_paths[0]), 'red', 'red_image.jpg')
        self.save_channel_image(self.process_images_common(image_paths[1]), 'green', 'green_image.jpg')
        if len(image_paths) == 3:
            self.save_channel_image(self.process_images_common(image_paths[2]), 'blue', 'blue_image.jpg')

        # 병합된 이미지를 고해상도로 저장
        new_size = (int(merged_image.width * 10), int(merged_image.height * 10))
        high_res_image = merged_image.resize(new_size, Image.Resampling.LANCZOS)
        high_res_image.save(output_filename, quality=100)
        print(f"{output_filename} saved successfully.")

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = FileDialogDemo()
    ex.show()
    sys.exit(app.exec_())
