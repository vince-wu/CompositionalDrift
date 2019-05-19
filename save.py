import pickle
from PyQt5.QtWidgets import QFileDialog
from parse import parseUI_Inputs
from generateUI import displayMessage

def save_state(self):
	parseUI_Inputs(self)
	state = {}
	state['Polymer System'] = self.numMonomers
	state['Model'] = self.model 
	state['Hold Composition'] = self.holdComposition
	state['Average DP'] = self.averageDP
	state['Percent Conversion'] = self.percentConversion
	state['Monomer Pool Size'] = self.poolSize
	state['Monomer Ratios'] = self.monomerRatios
	state['Reactivity Ratios'] = self.rrLoadList
	state['Graph Type'] = self.graphType
	state['Run Length Monomer'] = self.runLengthMonomer
	state['Rows to Show'] = self.numRows
	state['Animation'] = self.animate
	state['Legend'] = self.show_legend
	file = str(QFileDialog.getSaveFileName(self, "Select Directory")[0])
	if file:
		try:
			pickle_out = open(file, 'wb')
			pickle.dump(state, pickle_out)
			pickle_out.close()

		except Exception as e:
			text = "The filename you selected is invalid.	"
			displayMessage(self, "Error", text)


def load_state(self):
	file = str(QFileDialog.getOpenFileName(self, "Select Directory")[0])
	if file:
		try:
			pickle_in = open(file, 'rb')
			state = pickle.load(pickle_in)
		except Exception as e:
			text = "Invalid load file.		"
			displayMessage(self, "Error", text)	
			return
	else:
		return

	self.modelComboBox.setProperty('currentText', state['Model'])
	self.systemComboBox.setProperty('currentIndex', getSystemIndex(state['Polymer System']))
	self.holdCheckBox.setProperty('checked', state['Hold Composition'])
	self.dpSpinBox.setProperty('value', state['Average DP'])
	self.conversionDoubleSpinBox.setProperty('value', state['Percent Conversion'])
	self.poolSpinBox.setProperty('value', state['Monomer Pool Size'])
	self.graphComboBox.setProperty('currentText', state['Graph Type'])
	self.runLengthSpinBox.setProperty('value', state['Run Length Monomer'])
	self.rowsSpinBox.setProperty('value', state['Rows to Show'])
	self.animateComboBox.setProperty('currentText', getState(state['Animation']))
	self.legendComboBox.setProperty('currentText', getState(state['Legend']))
	setMonomerCompositions(self, state)
	setReactivities(self, state)

def setMonomerCompositions(self, state):
	numMonomers = state['Polymer System']
	monomerRatios = state['Monomer Ratios']
	for i in range(numMonomers):
		self.ratioDoubleSpinBoxList[i].setProperty('value', monomerRatios[i])

def setReactivities(self, state):
	numMonomers = state['Polymer System']
	rrLoadList = state['Reactivity Ratios']
	model = state['Model']

	if model == "Mayo-Lewis":
		for i in range(numMonomers):
			for j in range(numMonomers-1):
				self.rrDoubleSpinBox2DList[i][j].setProperty('value', rrLoadList[i][j])

	elif model == "Penultimate":
		for i in range(numMonomers):
			for j in range(numMonomers):
				for k in range(numMonomers):
					self.rrDoubleSpinBox3DList[i][j][k].setProperty('value', rrLoadList[i][j][k])


def getSystemIndex(text):
	return int(text) - 1

def getState(status):
	if status:
		return 'On'
	else:
		return 'Off'

