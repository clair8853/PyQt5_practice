import sys
from PyQt5.QtWidgets import QApplication
from ui import MyApp

def main():
    app = QApplication(sys.argv)
    ex = MyApp()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
