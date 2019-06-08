import math
import random
import numpy as np
import scipy.stats
import pyqtgraph as pg
from PyQt5.QtGui import *

from modules.graph import plot_rr
from modules.parse import parseRR_Inputs, rr_testAssertions
from modules.generateUI import displayMessage


def run_compute(self):

    parseRR_Inputs(self)

    inputsValid = rr_testAssertions(self)

    if not inputsValid:
        return

    self.statusbar.showMessage("Computing...")

    try:
        compute_rr(self, self.rr1_low, self.rr1_high, self.rr2_low, self.rr2_high, self.stepDensity, self.data)

    except Exception:
        displayMessage(self, "Error", "An error occured during the simulation. Please adjust the program inputs.    ")
        return

    if self.abort_rr: 
        self.statusbar.showMessage("Aborted.", 5000)
    else:
        self.statusbar.showMessage("Done.", 5000)

def compute_rr(self, rr1_low, rr1_high, rr2_low, rr2_high, stepDensity, data):
    """
    Input: lower and upper bounds for each reactivity ratio, increment, and a list of data sets.
    Each data set is a list of 3 elements: [Conversion, Monomer Fraction at Conversion, Initial Monomer Fraction]

    Output: returns a list of reactivity raios [r1, r2]
    """

    #create a matrix containing all combinations of reactivity ratios
    #count how many values of each reactivity ratio to test
    combinationsMatrix, rr1_counts, rr2_counts = create_rr_combinations(rr1_low, rr1_high, rr2_low, rr2_high, stepDensity)

    #find the greatest reactivity ratios tested (for plotting)
    rr1_max, rr2_max = combinationsMatrix[-1][-1]

    # print("rr1_low: ", rr1_low)
    # print("rr2_low: ", rr2_low)
    # print("rr1_max: ", rr1_max)
    # print("rr2_max: ", rr2_max)

    #initiate an empty matrix of zeros to hold least squares distance
    ssMatrix = [[0 for x in range(len(combinationsMatrix[0]))] for y in range(len(combinationsMatrix))] 

    #run the calculation for each set of data, updating ssMatrix every time
    num_sets = len(data)
    for i in range(num_sets):
        p = data[i][0]
        f = data[i][1]
        fZero = data[i][2]
        SSR(self, p,f,fZero, combinationsMatrix, ssMatrix)
        if self.abort_rr:
            return
        if i > 1:
            best_rr, best_ss =  get_best_rr(combinationsMatrix, ssMatrix)   
            contourMatrix = getContourData(i+1, ssMatrix, best_ss, rr1_counts, rr2_counts)
            plot_rr(self, contourMatrix, rr1_counts, rr2_counts, rr1_max, rr2_max)


    # # #find the best reactivty ratio pair based on least squares distance
    # best_rr, best_ss =  get_best_rr(combinationsMatrix, ssMatrix)

    # # #create a countour matrix of probabulity curves based on results
    # contourMatrix = getContourData(num_sets, ssMatrix, best_ss, rr1_counts, rr2_counts)

    # # #plot the contour matrix as a heatmap
    # plot_rr(self, contourMatrix, rr1_counts, rr2_counts, rr1_max, rr2_max)

    r1, r2 = best_rr
    self.r1_doubleSpinBox.setProperty("value", r1)
    self.r2_doubleSpinBox.setProperty("value", r2)

    return best_rr


def create_rr_combinations(rr1_low, rr1_high, rr2_low, rr2_high, stepDensity):
    """
    Creating a list with all combinations of reactivty ratios, leaving a spot for the sum of squares
    """
    incr = 1/stepDensity + 0.0000001
    combinationsMatrix = []
    rr1_counts = 0
    rr2_counts = 0
    rr1 = rr1_low
    rr2 = rr2_low
    index = 0
    while rr1 < rr1_high:
        combinationsMatrix.append([])
        rr2 = rr2_low
        while rr2 < rr2_high:
            if rr1_counts == 0:
                rr2_counts += 1
            rr2 += incr
            combinationsMatrix[index].append([rr1, rr2])
        index += 1

        rr1 += incr
        rr1_counts += 1

    return [combinationsMatrix, rr1_counts, rr2_counts]




"Equation 2"
#Conversion allowed by model
#Return out(which is either 0  or 1- something) (integer)
#r1, r2: 
#p: the experimental conversion	 
#f: fraction of a given monomer at conversion
#fzero intial monomer ratios
#out = calculated monomer conversion

def ModelConv(r1,r2,p,f,fZero):   
	#defines parameters for 
    alpha = r2 / (1 - r2)
    beta = r1 / (1 - r1)
    delta = (1 - r2) / (2 - r1 - r2)
    gamma = (1 - r1 * r2) / ((1 - r1) * (1 - r2))
    if f == 0:
        term1 = 0 
    else:
        term1 = alpha * math.log10(f / fZero)
    if f == 1:
        term2 = 0
    else:
        term2 = beta * math.log10((1 - f) / (1 - fZero))
    if  delta < f and delta < fZero:
        term3 = -999    
    elif  delta > f or delta > fZero:
        term3 = -999
    if  delta == f or delta == fZero:
        term3 = -999
                           
    else:
        term3 = gamma * math.log10((fZero - delta) / ((f - delta)))
    if term1 + term2 + term3 > 0: 
        out = 0
    else:
        out = 1 - math.exp(term1 + term2 + term3)
   
    return out 

#Tries to find fraction of monomers f that hits the data points best with allowed conversions -> (out)
#Composition of the monomers at these conversions as defined by modelconv
#Given r1, r2, intial monomer ratios, determine the conversion 
#returns a composition of first monomer (f) (number betweeen 0 and 1)
def ModelComp(r1,r2,p,fZero):
    MaxIter = 20
    #tolerance not used
    Tolerance = 0.01
    fAz = (1 - r2) / (2 - r1 - r2)
    if (r1 + r2 - 2) * fZero > (r2 - 1):
        fLo = 0
        if (fAz > 0) and (fAz <= fZero):
            fLo = fAz
        fHi = fZero
        NoIter = 0
        fMid = (fLo + fHi) / 2
        pMid = ModelConv(r1,r2, p, fMid, fZero)
        #Every time 
        while(NoIter < MaxIter):
            NoIter = NoIter + 1
            #defining a new upper and lower bound for the next iteration
            if pMid < p:
                fHi = fMid
            else:
                fLo = fMid
            fMid = (fLo + fHi) / 2
            pMid = ModelConv(r1, r2, pMid, fMid, fZero)

    else:
        fLo = fZero
        fHi = 1
        if (fAz < 1) and (fZero < fAz):
            fHi = fAz
        NoIter = 0
        fMid = (fLo + fHi) / 2
        pMid = ModelConv(r1, r2, p, fMid, fZero)
        while (NoIter < MaxIter):
            NoIter = NoIter + 1
            if pMid > p:
                fHi = fMid
            else:
                fLo = fMid
            fMid = (fLo + fHi) / 2
            pMid = ModelConv(r1, r2, pMid, fMid, fZero)
    return fMid

#eqn 14
def weight(p2,f2):
        return (4 * (1 - p2) ** 2) / (1 + (1 - 2 * f2) ** 2)

def Distance2(f2,p2,f1,p1):
        return ((p1 - p2) ** 2) + ((f1 - f2) ** 2) * weight(p2,f2)

#For a single r1, r2 pair, find the weighted distanced
def ModelCompConv(r1,r2,p,f,fZero):
    MaxCount = 200
    Count = 0
    Tolerance = 0.001
    p1 = p
    f1 = ModelComp(r1, r2, p1, fZero)
    f2 = fZero
    p2 = ModelConv(r1,r2, p, f2, fZero)

    def Distance2(f2,p2,f1,p1):
        return ((p1 - p2) ** 2) + ((f1 - f2) ** 2) * weight(p2,f2)

    d1 = Distance2(f1, p1, f2, p1)
    d2 = Distance2(f2, p2, f2, p1)
    while(abs(p2 - p1) < Tolerance) or Count > MaxCount:
        p0 = p1
        Count = Count + 1
        fGuess = (d2 * f1 + d1 * f2) / (d1 + d2)
        if fGuess < 0:
            fGuess = 0.01
        pGuess = ModelConv(r1, r2, p1, fGuess, fZero)
        dGuess = Distance2(fGuess, pGuess, f2, p1)
        if d1 > d2:
            p1 = pGuess
            f1 = fGuess
            d1 = dGuess
        else:
            p2 = pGuess
            f2 = fGuess
            d2 = dGuess
    
    if d1 < d2:
        ModelCompConv = p1
    else:
        ModelCompConv = p2
    return ModelCompConv



#p: experiemental conversion
#f: fraction of monomer at given conversion p
#fZero: initial fraction of monomers
def SSR(self, p,f,fZero, combinationsMatrix, ssMatrix):
    for i in range(len(combinationsMatrix)):
        for j in range(len(combinationsMatrix[i])):
            combination = combinationsMatrix[i][j]
            pHat = ModelCompConv(combination[0], combination[1], p, f, fZero)
            fHat = ModelComp(combination[0], combination[1], pHat, fZero)
            ssMatrix[i][j] += Distance2(fHat, pHat, f, p)

            #*** Processing Application Events, no impact on calculation ***
            if self.abort_rr:
                return
            QApplication.processEvents()


def get_best_rr(combinationsMatrix, ssMatrix):
    lst = np.array(ssMatrix)
    x,y = np.unravel_index(lst.argmin(), lst.shape)
    best_ss = ssMatrix[x][y]
    best_rr = combinationsMatrix[x][y]

    # best_ss, index = min((val, index) for (index, val) in enumerate(ssMatrix))
    # best_rr = combinationsMatrix[index]

    return [best_rr, best_ss]


def getContourData(num_sets, ssMatrix, best_ss, steps1, steps2):

    contourMatrix = np.array([[0 for x in range(len(ssMatrix[0]))] for y in range(len(ssMatrix))] )
    #find probability contours
    F95 = 1 + 2 / (num_sets - 2) * scipy.stats.f.ppf(0.05, 2, num_sets - 2)
    F90 = 1 + 2 / (num_sets - 2) * scipy.stats.f.ppf(0.1, 2, num_sets - 2)
    F70 = 1 + 2 / (num_sets - 2) * scipy.stats.f.ppf(0.3, 2, num_sets - 2)
    F50 = 1 + 2 / (num_sets - 2) * scipy.stats.f.ppf(0.5, 2, num_sets - 2)

    for i in range(len(ssMatrix)):
        for j in range(len(ssMatrix[0])):
            if ssMatrix[i][j] == best_ss:
                contourMatrix[i][j] += (best_ss)
            elif ssMatrix[i][j]/F95 < best_ss:
                contourMatrix[i][j] += 1
            elif ssMatrix[i][j]/F90 < best_ss:
                contourMatrix[i][j] += 2
            elif ssMatrix[i][j]/F70 < best_ss:
                contourMatrix[i][j] += 3
            elif ssMatrix[i][j]/F50 < best_ss:
                contourMatrix[i][j] += 4
            else:
                contourMatrix[i][j] += 5

    return contourMatrix
