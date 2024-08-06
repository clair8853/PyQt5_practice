import sys
from PyQt5.QtWidgets import QApplication, QWidget, QLabel, QDateTimeEdit, QVBoxLayout
from PyQt5.QtCore import QDateTime

class MyApp(QWidget):

    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        lbl = QLabel('QDateTimeEdit')

        DateTimeedit = QDateTimeEdit(self)
        DateTimeedit.setDateTime(QDateTime.currentDateTime())
        DateTimeedit.setDateTimeRange(QDateTime(1900, 1, 1, 00, 00, 00), QDateTime(2100, 1, 1, 00, 00, 00))
        DateTimeedit.setDisplayFormat('yyyy.MM.dd hh:mm:ss')

        vbox = QVBoxLayout()
        vbox.addWidget(lbl)
        vbox.addWidget(DateTimeedit)
        vbox.addStretch()

        self.setLayout(vbox)

        self.setWindowTitle('QDateTimeEdit')
        self.setGeometry(300,300,300,200)
        self.show()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = MyApp()
    sys.exit(app.exec_())