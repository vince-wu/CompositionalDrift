import sys 	
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *


import static.images_qr
from modules.MainForm import Ui_MainWindow
from modules.generateUI import setupDynamicUi
from modules.graph import plotData
from modules.simulate import run_simulation
from modules.visualize import setup_scene, draw_polymers
from modules.save import save_state, load_state
from modules.export import exportPolymerArray, exportImage
import modules.parse as parse


class MainWindow(QMainWindow, Ui_MainWindow):

	def __init__(self, *args, **kwargs):
		super(MainWindow, self).__init__(*args, **kwargs)


		self.setVars()
		self.setupUi(self)
		self.setWindowTitle("Compositional Drift v{}".format(self.version))
		setupDynamicUi(self, MainWindow)

		self.connectEvents()

	def setVars(self):
		self.version = "2.0"
		self.simulated = False
		self.simulation_running = False
		self.ratioLabelList = []
		self.ratioDoubleSpinBoxList = []
		self.formLayoutList = []
		

	def connectEvents(self):
		self.systemComboBox.currentTextChanged.connect(self.on_systemComboBox_changed)
		self.modelComboBox.currentTextChanged.connect(self.on_modelComboBox_changed)
		self.graphComboBox.currentTextChanged.connect(self.on_graphComboBox_changed)
		self.animateComboBox.currentTextChanged.connect(self.on_animateComboBox_changed)
		self.runLengthSpinBox.valueChanged.connect(self.on_runLengthComboBox_changed)
		self.rowsSpinBox.valueChanged.connect(self.on_rowsSpinBox_changed)
		self.legendComboBox.currentTextChanged.connect(self.on_legendComboBox_changed)
		self.simulateButton.clicked.connect(self.simulate)
		self.saveButton.clicked.connect(self.save)
		self.loadButton.clicked.connect(self.load)
		self.visualizeWindow.setContextMenuPolicy(Qt.CustomContextMenu)
		self.visualizeWindow.customContextMenuRequested.connect(self.openMenu)
		self.exportButton.clicked.connect(self.exportPolymers)


	
	def on_systemComboBox_changed(self, value):
		setupDynamicUi(self, MainWindow)


	def on_modelComboBox_changed(self, value):
		setupDynamicUi(self, MainWindow)

	def on_graphComboBox_changed(self, value):
		self.graphType = value
		self.runLengthMonomer = self.runLengthSpinBox.value()
		if self.simulated:
			plotData(self)

	def on_runLengthComboBox_changed(self, value):
		self.runLengthMonomer = value
		if self.simulated and self.graphType == "Run Length":
			plotData(self)

	def on_rowsSpinBox_changed(self, value):
		if self.simulated or self.simulation_running:
			rowDiff = value - self.numRows
			if rowDiff > 0:
				for i in range(rowDiff):
					self.polymer_animated_tracker.append(0)
				self.numRows = value
				draw_polymers(self)
			else:
				self.numRows = value
				setup_scene(self)
				draw_polymers(self)
		
	def on_animateComboBox_changed(self, value):
		self.animate = parse.get_animate_status(self)

	def on_legendComboBox_changed(self, value):
		self.show_legend = parse.get_legend_status(self)
		if self.simulated or self.simulation_running:
			plotData(self)

	def simulate(self):
		self.interrupt = True
		if not self.simulation_running:
			run_simulation(self)

	def save(self):
		save_state(self)

	def load(self):
		load_state(self)	

	def save_image(self):
		if self.simulated or self.simulation_running:
			exportImage(self)

	def exportPolymers(self):
		if self.simulated or self.simulation_running:
			exportPolymerArray(self)

	def keyPressEvent(self,event): 
	    if event.key()== Qt.Key_Return: 
	        self.simulate()

	def openMenu(self, position):
		menu = QMenu()
		saveImage_action = menu.addAction("Save Image")
		action = menu.exec_(self.visualizeWindow.mapToGlobal(position))
		if action == saveImage_action:
			self.save_image()
		saveImage_action.triggered.connect(self.save_image)


if __name__ == '__main__':

	app = QApplication([])
	app.setWindowIcon(QIcon(':/static/app.ico'))
	QApplication.setStyle(QStyleFactory.create('WindowsVista'))
	window = MainWindow()
	window.show()
	sys.exit(app.exec_())