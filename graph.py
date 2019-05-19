import sys 	
import pyqtgraph as pg 
import numpy as np
import analysis as analysis
from PyQt5.QtGui import *


from MainForm import Ui_MainWindow

def plotData(self):

	if self.graphType == "Monomer Occurrence":
		plotMonomerOccurrence(self)

	elif self.graphType == "Monomer Usage":
		plotMonomerUsage(self)

	elif self.graphType == "Run Length":
		plotRunLength(self, self.runLengthMonomer)

	elif self.graphType == "DP Distribution":
		plotDP(self)

	updateValues(self)

def plot_line_scatter(self, xList, yList, xLabel, yLabel, lgd_offset):
	clearGraph(self)

	self.graphWindow.setLabel('bottom', text=xLabel)
	self.graphWindow.setLabel('left', text=yLabel)

	if self.show_legend:
		self.lgd = self.graphWindow.addLegend(offset=lgd_offset)
		compositionList = analysis.getComposition(self, self.polymerArray)

	for i in range(self.numMonomers):

		x = xList[i]
		y = yList[i]

		self.graphWindow.setYRange(0, 1)
		self.graphWindow.setXRange(0, self.percentConversion)

		plt = self.graphWindow.plot(x, y, pen=(i,self.numMonomers))

		if self.show_legend:
			self.lgd.addItem(plt, "Monomer {} ({}% of polymer)".format(i+1, compositionList[i]))


	#Force the plot to update
	QApplication.processEvents()

"***Plotting Monomer Occurrence***"
def plotMonomerOccurrence(self):

	xList = []
	yList = []
	xLabel = "Conversion"
	yLabel = "Normalized Monomer Occurrence"
	lgd_offset = (10,10)

	for i in range(self.numMonomers):
		#extract the plotting data point for a single monomer
		single_monomerOccurrenceData = self.monomerOccurrenceList[i]

		#calculate the conversion index (x-axis array) 
		conversionIndex = range(len(single_monomerOccurrenceData))

		#adjustedPoolSize is the full size of the monomer pool, and i*numPolymers is how many monomers used cumulatively
		#at the data point. Dividing the latter by the former will give the conversion index
		conversionIndex = [(i+1)*self.numPolymers/self.adjustedPoolSize*100 for i in conversionIndex]

		xList.append(conversionIndex)
		yList.append(single_monomerOccurrenceData)

	plot_line_scatter(self, xList, yList, xLabel, yLabel, lgd_offset)


"***Plot Monomer Usage***"
def plotMonomerUsage(self):

	xList = []
	yList = []
	xLabel = "Conversion"
	yLabel = "Normalized Monomer Occurrence"
	lgd_offset = (-10,10)
	
	for i in range(self.numMonomers):
		#extract the plotting data point for a single monomer
		single_monomerRemainingData = self.monomerRemainingList[i]


		#calculate the conversion index (x-axis array) 
		conversionIndex = range(len(single_monomerRemainingData))

		#adjustedPoolSize is the full size of the monomer pool, and i*numPolymers is how many monomers used cumulatively
		#at the data point. Dividing the latter by the former will give the conversion index
		conversionIndex = [i*self.numPolymers/self.adjustedPoolSize*100 for i in conversionIndex]

		xList.append(conversionIndex)
		yList.append(single_monomerRemainingData)

	plot_line_scatter(self, xList, yList, xLabel, yLabel, lgd_offset)


"***Plot Run Length***"
def plotRunLength(self, monomerID):

	clearGraph(self)

	self.graphWindow.setLabel('bottom', text="Monomer {} Run Length".format(monomerID))
	self.graphWindow.setLabel('left', text="Normalized Counts")

	self.graphWindow.enableAutoRange('xy', True)
	#Use cached data if available
	if not self.simulation_running and self.cached_runLengthData[monomerID-1] != 0:
		x, y = self.cached_runLengthData[monomerID-1]

	else:
		x, y = analysis.getRunLengthHistData(self, monomerID)
		#caching data
		self.cached_runLengthData[monomerID-1] = [x,y]

	if len(y) == 0:
		return
	plt = self.graphWindow.plot(x, y, stepMode=True, fillLevel=0, brush=(monomerID-1, self.numMonomers))

	#Force the plot to update
	QApplication.processEvents()

"***Plot DP Distribution***"
def plotDP(self):

	clearGraph(self)

	self.graphWindow.setLabel('bottom', text="DP")
	self.graphWindow.setLabel('left', text="Counts")

	self.graphWindow.enableAutoRange('xy', True)

	x, y = analysis.DP_Distribution(self)

	if len(y) == 0:
		return

	plt = self.graphWindow.plot(x, y, stepMode=True, fillLevel=0, brush="y")


	#Force the plot to update
	QApplication.processEvents()

"***Update Number Avg DP, Weight Avg DP, and Dispersity Index***"
def updateValues(self):

	numberAverageDP = analysis.numberAverageDP(self)
	weightAverageDP = analysis.weightAverageDP(self, numberAverageDP)
	dispersityIndex = analysis.dispersityIndex(numberAverageDP, weightAverageDP)

	self.numAvDp_DoubleSpinBox.setProperty("value", numberAverageDP)
	self.weightAvDp_DoubleSpinBox.setProperty("value", weightAverageDP)
	self.dispersityDoubleSpinBox.setProperty("value", dispersityIndex)



"***Remove Legend Item***"
def removeLegend(self):

	try:
		self.lgd.scene().removeItem(self.lgd) 
	except Exception as e:
		return

"***Clear Graph, Including legend***"
def clearGraph(self):
	self.graphWindow.clear()
	removeLegend(self)