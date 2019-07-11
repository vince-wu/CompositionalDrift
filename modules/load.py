import csv
from PyQt5.QtWidgets import QFileDialog
from modules.generateUI import displayMessage

def load_data_sets(self):
	file = str(QFileDialog.getOpenFileName(self, "Select Directory")[0])
	dataSets = []
	if file:
		try:
			with open(file) as csv_file:
				csv_reader = csv.reader(csv_file, delimiter=',')
				for row in csv_reader:
					dataSets.append([float(row[0]), float(row[1]), float(row[2])])

				copy_data(self, dataSets)

		except Exception as e:
			text = "Invalid load file.		"	
			displayMessage(self, "Error", text)	
			return
	else:
		return

def copy_data(self, dataSets):
	numSets = len(dataSets)

	self.numSetsSpinBox.setProperty('value', numSets)

	for i in range(numSets):
		self.setDataDoubleSpinBoxList[i][0].setProperty('value', dataSets[i][0])
		self.setDataDoubleSpinBoxList[i][1].setProperty('value', dataSets[i][1])
		self.setDataDoubleSpinBoxList[i][2].setProperty('value', dataSets[i][2])


