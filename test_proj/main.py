import sys
from PyQt5.QtWidgets import QApplication
from ui import MyApp
from config_reader import load_config

def main():
    config = load_config()
    app = QApplication(sys.argv)
    ex = MyApp(config)
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
