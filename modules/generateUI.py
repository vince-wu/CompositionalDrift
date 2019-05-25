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
		self.ratioLabelList[i].setToolTip("Initial molar ratio of monomer {}".format(i+1))

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
					if numMonomers == 2:
						self.rrLabel2DList[i][jIndex].setToolTip("Proportional to monomer {} homopolymerization probability".format(i+1))
					else:
						self.rrLabel2DList[i][jIndex].setToolTip("Inversely proportional to probability of monomer {} binding to monomer {}".format(i+1, j+1))
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

def rr_setupDynamicUi(self):
	num_data_sets = self.numSetsSpinBox.value()

	for formLayout in self.rr_formLayoutList:
		clearLayoutItems(self, formLayout)
		self.verticalLayout_10.removeItem(formLayout)
		formLayout.deleteLater()
		formLayout = None

	for line in self.lineList:
		self.verticalLayout_10.removeWidget(line)
		line.deleteLater()
		line = None

	if self.spacerItem:
		self.verticalLayout_10.removeItem(self.spacerItem)
		self.spacerItem = None


	self.rr_formLayoutList = []
	self.setDataLabelList = []
	self.setDataDoubleSpinBoxList = []
	self.lineList = []



	for i in range(num_data_sets):

		self.setDataLabelList.append([])
		self.setDataDoubleSpinBoxList.append([])

		self.rr_formLayoutList.append(QtWidgets.QFormLayout())

		self.setDataLabelList[i].append(QtWidgets.QLabel(self.scrollAreaWidgetContents_5))
		self.setDataLabelList[i][0].setText("Set {} Conversion".format(i+1))
		self.rr_formLayoutList[i].setWidget(0, QtWidgets.QFormLayout.LabelRole, self.setDataLabelList[i][0])
		self.setDataDoubleSpinBoxList[i].append(QtWidgets.QDoubleSpinBox(self.scrollAreaWidgetContents_5))
		self.rr_formLayoutList[i].setWidget(0, QtWidgets.QFormLayout.FieldRole, self.setDataDoubleSpinBoxList[i][0])

		self.setDataLabelList[i].append(QtWidgets.QLabel(self.scrollAreaWidgetContents_5))
		self.setDataLabelList[i][1].setText("Set {} Monomer Fraction".format(i+1))
		self.rr_formLayoutList[i].setWidget(1, QtWidgets.QFormLayout.LabelRole, self.setDataLabelList[i][1])
		self.setDataDoubleSpinBoxList[i].append(QtWidgets.QDoubleSpinBox(self.scrollAreaWidgetContents_5))
		self.rr_formLayoutList[i].setWidget(1, QtWidgets.QFormLayout.FieldRole, self.setDataDoubleSpinBoxList[i][1])

		self.setDataLabelList[i].append(QtWidgets.QLabel(self.scrollAreaWidgetContents_5))
		self.setDataLabelList[i][2].setText("Set {} Initial Monomer Fraction".format(i+1))
		self.rr_formLayoutList[i].setWidget(2, QtWidgets.QFormLayout.LabelRole, self.setDataLabelList[i][2])
		self.setDataDoubleSpinBoxList[i].append(QtWidgets.QDoubleSpinBox(self.scrollAreaWidgetContents_5))
		self.rr_formLayoutList[i].setWidget(2, QtWidgets.QFormLayout.FieldRole, self.setDataDoubleSpinBoxList[i][2])

		self.verticalLayout_10.addLayout(self.rr_formLayoutList[i])

		self.lineList.append(QtWidgets.QFrame(self.scrollAreaWidgetContents_5))
		self.lineList[i].setFrameShape(QtWidgets.QFrame.HLine)
		self.lineList[i].setFrameShadow(QtWidgets.QFrame.Sunken)
		self.lineList[i].setObjectName("line_2")
		self.verticalLayout_10.addWidget(self.lineList[i])	

	self.spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
	self.verticalLayout_10.addItem(self.spacerItem)

def clearLayout(layout):
  while layout.count():
    child = layout.takeAt(0)
    if child.widget():
      child.widget().deleteLater()

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
