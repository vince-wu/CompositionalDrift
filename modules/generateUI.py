from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from modules.export import exportImage

def setupDynamicUi(self, MainWindow):
	#define local values for numMonomers and model so as to not interfere with any currently running simulations
	model = self.modelComboBox.currentText()
	numMonomers = get_numMonomer(self.systemComboBox.currentText())

	#Set limits for runLengthSpinBox
	self.runLengthSpinBox.setMaximum(numMonomers)

	#for Penultimate, only 2-monomer system is allowed
	if model == "Penultimate":
		self.systemComboBox.setProperty("currentIndex", 1)
		numMonomers = 2

	#Remove old QtWidgets
	for labelWidget in self.ratioLabelList:
		self.monomerRatio_layout.removeWidget(labelWidget)
		labelWidget.deleteLater()
		labelWidget = None

	for doubleSpinBoxWidget in self.ratioDoubleSpinBoxList:
		self.monomerRatio_layout.removeWidget(doubleSpinBoxWidget)
		doubleSpinBoxWidget.deleteLater()
		doubleSpinBoxWidget = None

	for formLayout in self.formLayoutList:
		clearLayoutItems(self, formLayout)
		self.horizontalLayout.removeItem(formLayout)
		formLayout.deleteLater()
		formLayout = None

	#create lists to hold widget objects
	self.ratioLabelList = []
	self.ratioDoubleSpinBoxList = []
	self.formLayoutList = []
	self.rrLabel2DList = []
	self.rrDoubleSpinBox2DList = []
	self.rrLabel3DList = []
	self.rrDoubleSpinBox3DList = []




	#generate dynamic input widgets
	_translate = QCoreApplication.translate

	#generate input widgets for monomer ratios, set up layouts for reactivity ratios
	for i in range(numMonomers):
		self.ratioLabelList.append(QtWidgets.QLabel(self.inputLayout))
		self.ratioLabelList[i].setObjectName("monomerRatioLabel" + str(i))
		self.monomerRatio_layout.setWidget(i, QtWidgets.QFormLayout.LabelRole, self.ratioLabelList[i])
		self.ratioDoubleSpinBoxList.append(QtWidgets.QDoubleSpinBox(self.inputLayout))
		self.ratioDoubleSpinBoxList[i].setObjectName("monomerRatioBox" + str(i))
		self.ratioDoubleSpinBoxList[i].setProperty("value", 1)
		self.monomerRatio_layout.setWidget(i, QtWidgets.QFormLayout.FieldRole, self.ratioDoubleSpinBoxList[i])
		self.ratioLabelList[i].setText(_translate("MainWindow", "Monomer {} Ratio".format(i+1)))

		#set up layouts for reactivity ratios
		self.formLayoutList.append(QtWidgets.QFormLayout())
		self.formLayoutList[i].setObjectName("formLayoutRR{}".format(i+1))

	#Case for Mayo-Lewis Model
	if model == "Mayo-Lewis":
		for i in range(numMonomers):
			self.rrLabel2DList.append([])
			self.rrDoubleSpinBox2DList.append([])
			#generate input widgets for reactivity ratios
			jIndex = 0
			for j in range(numMonomers):
				if j != i:
					self.rrLabel2DList[i].append(QtWidgets.QLabel(self.inputLayout))
					self.rrLabel2DList[i][jIndex].setObjectName("rr{}{}Label".format(i+1,j+1))
					self.rrLabel2DList[i][jIndex].setText(_translate("MainWindow", "rr{}{}".format(i+1,j+1)))
					self.formLayoutList[i].setWidget(jIndex, QtWidgets.QFormLayout.LabelRole, self.rrLabel2DList[i][jIndex])
					self.rrDoubleSpinBox2DList[i].append(QtWidgets.QDoubleSpinBox(self.inputLayout))
					self.rrDoubleSpinBox2DList[i][jIndex].setObjectName("rr{}{}DoubleSpinBox".format(i+1,j+1))
					self.rrDoubleSpinBox2DList[i][jIndex].setProperty("value", 1)
					self.rrDoubleSpinBox2DList[i][jIndex].setMaximum(100)
					self.formLayoutList[i].setWidget(jIndex, QtWidgets.QFormLayout.FieldRole, self.rrDoubleSpinBox2DList[i][jIndex])
					jIndex += 1
			self.horizontalLayout.addLayout(self.formLayoutList[i])

	#case for Penultimate Model
	elif model == "Penultimate":
		for i in range(numMonomers):
			self.rrLabel3DList.append([])
			self.rrDoubleSpinBox3DList.append([])
			positionIndex = 0
			for j in range(numMonomers):
				self.rrLabel3DList[i].append([])
				self.rrDoubleSpinBox3DList[i].append([])
				for k in range(numMonomers):
					self.rrLabel3DList[i][j].append(QtWidgets.QLabel(self.inputLayout))
					self.rrLabel3DList[i][j][k].setObjectName("rr{}{}{}Label".format(i+1,j+1,k+1))
					self.rrLabel3DList[i][j][k].setText(_translate("MainWindow", "rr{}{}{}".format(i+1,j+1,k+1)))
					self.formLayoutList[i].setWidget(positionIndex, QtWidgets.QFormLayout.LabelRole, self.rrLabel3DList[i][j][k])
					self.rrDoubleSpinBox3DList[i][j].append(QtWidgets.QDoubleSpinBox(self.inputLayout))
					self.rrDoubleSpinBox3DList[i][j][k].setObjectName("rr{}{}{}DoubleSpinBox".format(i+1,j+1,k+1))
					self.rrDoubleSpinBox3DList[i][j][k].setProperty("value", 1)
					self.rrDoubleSpinBox3DList[i][j][k].setMaximum(100)
					self.formLayoutList[i].setWidget(positionIndex, QtWidgets.QFormLayout.FieldRole, self.rrDoubleSpinBox3DList[i][j][k])
					positionIndex += 1
			self.horizontalLayout.addLayout(self.formLayoutList[i])

#deletes all items within an layout
def clearLayoutItems(self, layout):
	if layout is not None:
	    while layout.count():
	        item = layout.takeAt(0)
	        widget = item.widget()
	        if widget is not None:
	            widget.deleteLater()
	        else:
	            self.clearLayout(item.layout())


def displayMessage(self, type, text):
	msg = QMessageBox()

	if type == "Error":
		msg.setIcon(QMessageBox.Warning)
		msg.setWindowTitle("Error")
	msg.setText(text)
	msg.setStandardButtons(QMessageBox.Ok)
	msg.exec_()

def get_numMonomer(str):
	if str[:2] == "10":
		return 10
	else:
		return int(str[0])
