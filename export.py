import csv
from PyQt5.QtWidgets import QFileDialog
from PyQt5.QtGui import QPainter, QPixmap

def exportPolymerArray(self):
	csvData = []
	for polymer in self.polymerArray:
		csvData.append(polymer.asArray())
	file = str(QFileDialog.getSaveFileName(self, "Select Directory", filter = 'csv(*.csv)')[0])
	if file:
		with open(file, 'w') as csvFile:
			writer = csv.writer(csvFile)
			writer.writerows(csvData)
		csvFile.close()

def exportImage(self):
	#find image dimensions
	imgWidth = self.scene.width()
	imgHeight = self.scene.height()
	#copy scene onto QPixmap obj
	pix = QPixmap(imgWidth, imgHeight)
	painter = QPainter(pix)
	self.scene.render(painter)
	painter.end()
	#Save as png
	file = str(QFileDialog.getSaveFileName(self, "Select Directory", filter = 'png(*.png)')[0])
	if file:
		pix.save(file, "PNG")