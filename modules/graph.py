import sys 	
import pyqtgraph as pg 
import numpy as np
from PyQt5.QtGui import *
from PyQt5 import QtCore, QtGui, QtWidgets

import modules.analysis as analysis
from modules.Heatmap import Heatmap
from modules.Legend import Legend
from modules.MainForm import Ui_MainWindow
from modules.CustomViewBox import CustomViewBox


"***Set Up Graph Screen***"
def setUpGraph(self):
	self.graphWindow = pg.PlotWidget(self.tab)
	sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
	sizePolicy.setHorizontalStretch(70)
	sizePolicy.setVerticalStretch(8)
	sizePolicy.setHeightForWidth(self.graphWindow.sizePolicy().hasHeightForWidth())
	self.graphWindow.setSizePolicy(sizePolicy)
	self.graphWindow.setFrameShape(QtWidgets.QFrame.NoFrame)
	self.graphWindow.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
	self.graphWindow.setDragMode(QtWidgets.QGraphicsView.ScrollHandDrag)
	self.horizontalLayout_6.insertWidget(0, self.graphWindow)
	self.horizontalLayout_6.removeWidget(self.displayView)
	self.displayView.deleteLater()
	self.displayView = None

"***Main Plot Function***"

def plotData(self):

	if self.displayView:
		setUpGraph(self)

	self.graphWindow.showAxis('left', True)
	self.graphWindow.showAxis('bottom', True)

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
		# compositionList = analysis.getComposition(self, self.polymerArray)
		monomer_dist = self.reaction.monomer_distribution()

	for i in range(self.numMonomers):

		x = xList[i]
		y = yList[i]

		self.graphWindow.setYRange(0, 1)
		self.graphWindow.setXRange(0, self.percentConversion)

		plt = self.graphWindow.plot(x, y, pen=(i,self.numMonomers))

		if self.show_legend:
			self.lgd.addItem(plt, "Monomer {} ({}% of polymer)".format(i+1, monomer_dist[i]))


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
		single_monomerOccurrenceData = self.monomer_species_occurrences[i]

		#calculate the conversion index (x-axis array) 
		conversionIndex = range(len(single_monomerOccurrenceData))

		#adjustedPoolSize is the full size of the monomer pool, and i*numPolymers is how many monomers used cumulatively
		#at the data point. Dividing the latter by the former will give the conversion index
		conversionIndex = [(i+1)*self.numPolymers/self.adjustedPoolSize*100 for i in conversionIndex]

		xList.append(self.conversion_index)
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
		single_monomerRemainingData = self.monomer_species_remaining[i]


		#calculate the conversion index (x-axis array) 
		conversionIndex = range(len(single_monomerRemainingData))

		#adjustedPoolSize is the full size of the monomer pool, and i*numPolymers is how many monomers used cumulatively
		#at the data point. Dividing the latter by the former will give the conversion index
		conversionIndex = [(i+1)*self.numPolymers/self.adjustedPoolSize*100 for i in conversionIndex]

		xList.append(self.conversion_index)
		yList.append(single_monomerRemainingData)


	plot_line_scatter(self, xList, yList, xLabel, yLabel, lgd_offset)


"***Plot Run Length***"
def plotRunLength(self, monomerID):

	clearGraph(self)

	self.graphWindow.setLabel('bottom', text="Monomer {} Run Length".format(monomerID))
	self.graphWindow.setLabel('left', text="Normalized Counts")

	self.graphWindow.enableAutoRange('xy', True)
	#Use cached data if available
	if not self.simulation_running and self.cached_runLengthData[monomerID] != 0:
		x, y = self.cached_runLengthData[monomerID]

	else:
		x, y = analysis.getRunLengthHistData(self, monomerID)
		#caching data
		self.cached_runLengthData[monomerID] = [x,y]

	if len(y) == 0:
		return
	plt = self.graphWindow.plot(x, y, stepMode=True, fillLevel=0, brush=(monomerID, self.numMonomers))

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

	if self.numMonomers != 2:
		self.lambdaDoubleSpinBox.setProperty("value", -2)

	elif self.numMonomers == 2 and min([len(polymer) for polymer in self.polymerArray]) > 1:
		lambdaValue = analysis.calculate_theta(self)
		self.lambdaDoubleSpinBox.setProperty("value", lambdaValue)

	else:
		self.lambdaDoubleSpinBox.setProperty("value", 0)

"***Remove Legend Item***"
def removeLegend(lgd):

	try:
		lgd.scene().removeItem(lgd) 
	except Exception as e:
		return

"***Clear Graph, Including legend***"
def clearGraph(self):
	self.graphWindow.clear()
	removeLegend(self.lgd)

"***Set up RR PlotWidget***"
def setUpRRPlot(self):
	self.vb = CustomViewBox()
	self.rrPlotWidget = pg.PlotWidget(viewBox=self.vb, enableMenu=False)
	self.rrPlotWidget.enableAutoRange('xy', False)
	self.verticalLayout_7.insertWidget(0, self.rrPlotWidget)
	self.verticalLayout_7.removeWidget(self.rr_textBrowser)
	self.rr_textBrowser.deleteLater()

"***Plot RR Heatmap ***"
def plot_rr(self, data, x_counts, y_counts, rr1_max, rr2_max):

	#Find scaling factors for heatmap plot
	self.scaleX = x_counts/(rr1_max)
	self.scaleY = y_counts/(rr2_max)

	#Find the min and max of each axis
	x_min = self.rr1_low*self.scaleX
	y_min = self.rr2_low*self.scaleY
	x_max = rr1_max*self.scaleX
	y_max = rr2_max*self.scaleY

	#Create a binding rectangle bases on dimensions
	bindingRect = QtCore.QRectF(x_min, y_min, x_max - x_min, y_max - y_min)

	if not self.htmap:
		#set up Plotwidget
		setUpRRPlot(self)

		#create legend
		creat_rr_legend(self)

		#Creat Heatmap
		self.htmap = Heatmap(data)
		self.rrPlotWidget.addItem(self.htmap)
		self.htmap.setRect(bindingRect)
		self.rrPlotWidget.setLabel('bottom', text="r1")
		self.rrPlotWidget.setLabel('left', text="r2")

	else:
		#Update heatmap and dimensions
		self.htmap.updateImage(data)
		self.htmap.setRect(bindingRect)

	#Set axis ranges
	self.rrPlotWidget.setXRange(x_min, x_max)
	self.rrPlotWidget.setYRange(y_min, y_max)

	#Adjust axis scales
	x_axis = self.rrPlotWidget.getAxis('bottom')
	y_axis = self.rrPlotWidget.getAxis('left')
	x_axis.setScale(1/self.scaleX)
	y_axis.setScale(1/self.scaleY)

	#Force update events
	QApplication.processEvents()

"***Create legend for Heatmap***"
#Hacky but it works...
def creat_rr_legend(self):
	#Create Legend
	self.rr_lgd = self.rrPlotWidget.addLegend()
	gradients =[(37,52,148), (44,127,184), (127,205,187), (199,233,180), (255,255,204)]
	styles = []
	for color in gradients:
		styles.append(self.rrPlotWidget.plot([0], pen=None, symbolPen = (0,0,0), symbolBrush=color, symbol='s', symbolSize=15))

	self.rr_lgd.addItem(styles[0], "95%")
	self.rr_lgd.addItem(styles[1], "90%")
	self.rr_lgd.addItem(styles[2], "70%")
	self.rr_lgd.addItem(styles[3], "50%")
	self.rr_lgd.addItem(styles[4], "0%")

	for style in styles:
		self.rrPlotWidget.removeItem(style)