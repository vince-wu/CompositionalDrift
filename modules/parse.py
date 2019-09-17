import math
from modules.generateUI import displayMessage

#Extracts all values needed to simulate polymerization from GUI
def parseUI_Inputs(self):
	self.model = self.modelComboBox.currentText()
	self.numMonomers = get_numMonomer(self.systemComboBox.currentText())
	self.chainTransferPercentage = self.transferSpinBox.value()
	self.holdComposition = getHoldComposition(self)
	self.rrList = getReactivityRatios(self)
	self.rrLoadList = getLoadReactivityRatios(self)
	self.rateConstantList = getRateConstants(self)
	self.monomerRatios = getMonomerRatios(self)
	self.averageDP = self.dpSpinBox.value()
	self.percentConversion = self.conversionDoubleSpinBox.value()
	self.poolSize = self.poolSpinBox.value()
	self.graphType = self.graphComboBox.currentText()
	self.animate = get_animate_status(self.animateComboBox)
	self.numRows = self.rowsSpinBox.value()	
	self.runLengthMonomer = self.runLengthSpinBox.value()
	self.cached_runLengthData = [0]*self.numMonomers
	self.show_legend = get_legend_status(self)
	self.interrupt = False
	self.size_set = False

def parseRR_Inputs(self):
	self.rr1_low = self.r1_low_doubleSpinBox.value()
	self.rr1_high = self.r1_high_doubleSpinBox.value()
	self.rr2_low = self.r2_low_doubleSpinBox.value()
	self.rr2_high = self.r2_high_doubleSpinBox.value()
	self.stepDensity = self.densitySpinBox.value()
	self.animateRR = get_animate_status(self.rr_animateComboBox)
	self.data = getDataSets(self)
	self.abort_rr = False




#Returns a nested array of reaqctvity ratio numbers (double) extracted from GUI inputs
#These values are used for ONLY initiation in Mayo Lewis (NOT PROPAGATION), and for


def getReactivityRatios(self):
	
	rrList = []

	#case for Mayo Lewis
	if self.model == "Mayo-Lewis":

		for i in range(self.numMonomers):
			rrList.append([])
			jIndex = 0

			for j in range(self.numMonomers):

				if i == j:
					rrList[i].append(1)

				else:
					rrValue = self.rrDoubleSpinBox2DList[i][jIndex].value()
					rrList[i].append(rrValue)
					jIndex += 1
	#case for Penultimate
	elif self.model == "Penultimate":

		for i in range(self.numMonomers):
			rrList.append([])

			for j in range(self.numMonomers):
				rrList[i].append([])

				for k in range(self.numMonomers):
					rrValue = self.rrDoubleSpinBox3DList[i][j][k].value()
					rrList[i][j].append(rrValue)

	#print("rrlist: ", rrList)
	return rrList

#Returns a nested array of rate constants extracted form GUI inputs
#These values are used in the PROPAGATION step
def getRateConstants(self):
	rateConstantList = []

	#case for Mayo Lewis
	if self.model == "Mayo-Lewis":

		for i in range(self.numMonomers):
			rateConstantList.append([])
			jIndex = 0

			for j in range(self.numMonomers):

				if i == j:
					rateConstantList[i].append(1)

				else:
					rrValue = self.rrDoubleSpinBox2DList[i][jIndex].value()
					if rrValue == 0:
						rateConstantList[i].append(math.inf)
					else:
						rateConstantList[i].append(1/rrValue)
					jIndex += 1

	#case for Penultimate
	elif self.model == "Penultimate":

		for i in range(self.numMonomers):
			rateConstantList.append([])

			for j in range(self.numMonomers):
				rateConstantList[i].append([])

				for k in range(self.numMonomers):
					rrValue = self.rrDoubleSpinBox3DList[i][j][k].value()
					rateConstantList[i][j].append(rrValue)

	#print("rateConstantList: ", rateConstantList)
	return rateConstantList

def getLoadReactivityRatios(self):
	rrLoadList = []

	if self.model == "Mayo-Lewis":
		for i in range(self.numMonomers):
			rrLoadList.append([])

			for j in range(self.numMonomers-1):
				rrValue = self.rrDoubleSpinBox2DList[i][j].value()
				rrLoadList[i].append(rrValue)

	elif self.model == "Penultimate":
		for i in range(self.numMonomers):
			rrLoadList.append([])

			for j in range(self.numMonomers):
				rrLoadList[i].append([])

				for k in range(self.numMonomers):
					rrValue = self.rrDoubleSpinBox3DList[i][j][k].value()
					rrLoadList[i][j].append(rrValue)

	return rrLoadList
#returns a list of monomer ratios (doubles) extracted from GUI inputs

def getMonomerRatios(self):

	monomerRatios = []

	for i in range(self.numMonomers):
		ratio = self.ratioDoubleSpinBoxList[i].value()
		monomerRatios.append(ratio)

	#print("ratio: ", monomerRatios)
	return monomerRatios

def getHoldComposition(self):
	if self.holdCheckBox.isChecked():
		return True
	else:
		return False

def get_numMonomer(str):
	if str[:2] == "10":
		return 10
	else:
		return int(str[0])

def get_animate_status(comboBox):
	value = comboBox.currentText()
	if value == "On":
		return True
	elif value == "Off":
		return False

def get_legend_status(self):
	value = self.legendComboBox.currentText()
	if value == "On":
		return True
	elif value == "Off":
		return False

def getDataSets(self):
	dataSet = []
	index = 0
	for data in self.setDataDoubleSpinBoxList:
		dataSet.append([])
		dataSet[index].append(data[0].value())
		dataSet[index].append(data[1].value())
		dataSet[index].append(data[2].value())
		index += 1
	return dataSet

#eqn 14
def weight(p2,f2):
        return (4 * (1 - p2) ** 2) / (1 + (1 - 2 * f2) ** 2)

def Distance2(f2,p2,f1,p1):
        return ((p1 - p2) ** 2) + ((f1 - f2) ** 2) * weight(p2,f2)

def testAssertions(self):
	inputsValid = True	

	if self.poolSize < self.averageDP*self.percentConversion/100:
		displayMessage(self, "Error", "Monomer Pool Size is too small!	")
		inputsValid = False

	elif self.numPolymers < 1:
		displayMessage(self, "Error", "DP or Monomer Pool Size is too small!	")
		inputsValid = False

	elif int(self.averageDP * self.percentConversion / 100) <= 0:
		displayMessage(self, "Error", "DP or Conversion cannot be zero!		")
		inputsValid = False

	elif min(self.monomerRatios) == 0:
		displayMessage(self, "Error", "Monomer Ratios cannot be zero!	")
		inputsValid = False

	return inputsValid

def rr_testAssertions(self):
	inputsValid = True

	for data in self.data:
		conversion = data[0]
		fraction = data[1]
		init_fraction = data[2]

		if conversion <= 0.0001:
			displayMessage(self, "Error", "Conversion too close to zero!	")
			inputsValid = False
			break

		elif init_fraction == 0:
			displayMessage(self, "Error", "Initial Monomer Fraction cannot be zero!	")
			inputsValid = False
			break

	if self.rr1_low >= self.rr1_high or self.rr2_low >= self.rr1_high:
		displayMessage(self, "Error", "Lower bound must be less than upper bound!	")
		inputsValid = False

	return inputsValid








