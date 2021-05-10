import sys
from PyQt5 import QtWidgets, uic

class Ui(QtWidgets.QDialog):
    def __init__(self):
        super(Ui, self).__init__() # Call the inherited classes __init__ method
        uic.loadUi("hes_off_gui.ui", self) # Load the .ui file
        self.setWindowTitle("HES-OFF research-based prototype tool")
        self.show() # Show the GUI

app = QtWidgets.QApplication(sys.argv)
window = Ui()
app.exec_()