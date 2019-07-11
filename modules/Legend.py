import pyqtgraph as pg
from PyQt5 import QtGui

class Legend(pg.LegendItem):
	def __init__(self, size=None, offset=None):
		super(Legend, self).__init__()
	def addItem(self, item, name):
		label = QtGui.QLabel(name)
		sample = item
		row = self.layout.rowCount()
		self.items.append((sample, label))
		self.layout.addItem(sample, row, 0)
		self.layout.addItem(label, row, 1)
		self.updateSize()