from PyQt5 import QtWidgets, QtCore
from PyQt5.QtGui import QPainter, QBrush, QPen, QPixmap
from PyQt5.QtWidgets import QApplication
from PyQt5.QtCore import Qt 
import pyqtgraph as pg

def setup_scene(self):
	self.visualizeWindow.items().clear()
	self.scene = QtWidgets.QGraphicsScene()
	self.visualizeWindow.setScene(self.scene)
	self.polymer_animated_tracker = [0]*self.numRows

def draw_polymers(self):
	#self.scene.clear()
	if not self.polymerArray:
		return
		
	if not self.size_set:
		#find widget width and height
		width =  self.visualizeWindow.geometry().width()

		#set size of squares
		maxSize = 15
		minSize = 7

		self.size = max(min(int( width / (self.averageDP*1.2*self.percentConversion/100)), maxSize), minSize)

		self.size_set = True

	#print("size: ", size)

	pen = QPen(QtCore.Qt.black)

	for i in range(self.numRows):
		#Get the index to start drawing at
		index_to_start = self.polymer_animated_tracker[i]
		#get the polymer to draw
		polymer = self.polymerArray[i].asArray()
		#get which monomers to draw
		monomers_to_draw = polymer[index_to_start:]
		#index to keep track of x location
		j = index_to_start
		#draw the monomers
		for monomerID in monomers_to_draw:
			brush = QBrush(pg.mkColor((monomerID-1, self.numMonomers)))
			#brush = QBrush(QtCore.Qt.red)
			self.scene.addRect(j*self.size, i*self.size, self.size, self.size, pen, brush)
			j += 1

		self.polymer_animated_tracker[i] += len(monomers_to_draw)

	#Force the plot to update
	QApplication.processEvents()







