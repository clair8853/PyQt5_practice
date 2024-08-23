from setuptools import setup
import sys

sys.setrecursionlimit(10000)

APP = ['/Users/mmbd7267/Development/PyQt5_practice/Marmoset_app_mac/main.py']
DATA_FILES = []
OPTIONS = {
    'argv_emulation': True, 
    'packages': ['PyQt5', 'matplotlib', 'numpy', 'pandas', 'scanpy', 'squidpy'], 
    'includes': ['matplotlib.backends.backend_qt5agg'], 
}

setup(
    app=APP,
    data_files=DATA_FILES,
    options={'py2app': OPTIONS},
    setup_requires=['py2app'],
)
