import pyqtgraph as pg
import numpy as np

class Heatmap(pg.ImageItem):
    def __init__(self, image=None, cmap=None):

        if image is not None:
            self.image = image
        else:
            self.image = np.zeros((500, 500))

        self.lookuptbl = [(37,52,148), (44,127,184), (127,205,187), (199,233,180), (255,255,204)]

        pg.ImageItem.__init__(self, self.image, lut=self.lookuptbl)

    def updateImage_(self, image):
        self.image = image
        self.render()
        self.update()