import sys
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QLineEdit, QListWidget

class SearchApp(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.layout = QVBoxLayout()

        self.search_bar = QLineEdit(self)
        self.search_bar.setPlaceholderText("검색어를 입력하세요")
        self.search_bar.textChanged.connect(self.filter_list)
        self.layout.addWidget(self.search_bar)

        self.list_widget = QListWidget(self)
        self.items = ["apple", "banana", "cherry", "date", "elderberry", "fig", "grape"]
        self.list_widget.addItems(self.items)
        self.layout.addWidget(self.list_widget)

        self.setLayout(self.layout)
        self.setWindowTitle('리스트 검색 예제')
        self.setGeometry(300, 300, 300, 200)
        self.show()

    def filter_list(self):
        filter_text = self.search_bar.text().lower()
        self.list_widget.clear()
        filtered_items = [item for item in self.items if filter_text in item.lower()]
        self.list_widget.addItems(filtered_items)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = SearchApp()
    sys.exit(app.exec_())
