import xlwt
import xlrd
import math
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns; sns.set()
import csv
import scipy.stats

StR1 = 0.0
EndR1 = 1.50
StR2 = 0.0
EndR2 = 1.50
DR = 0.011

"""Creating a list with all combinations of reactivty ratios, leaving a spot for the sum of squares"""
r1=StR1
r2=StR2
Steps1=0
Steps2=0
while r1 < EndR1:
    Steps1 +=1
    r1 += DR
r1=StR1
SS = []
while r2 < EndR2:
    Steps2+=1
    r1 = StR1
    while r1 < EndR1:
        r1 += DR
        SS.append([r1,r2,0])
    r2 += DR

r1=StR1
r2=StR2
Steps1=0
Steps2=0
SeaRow = []
SeaColumn=[]

while r1 < EndR1:
	#Fix the rounding: so number on plots wont overlap
    SeaColumn.append(str(round(r1,3)))
    Steps1 +=1
    r1 += DR
while r2 < EndR2:
    SeaRow.append(str(round(r2,3)))
    Steps2 +=1
    r2 += DR


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


points = 0

#p: experiemental conversion
#f: fraction of monomer at given conversion p
#fZero: initial fraction of monomers
def SSR(p,f,fZero):
    global points
    points +=1
    for x in range(len(SS)):
        pHat = ModelCompConv(SS[x][0], SS[x][1], p, f, fZero)
        fHat = ModelComp(SS[x][0], SS[x][1], pHat, fZero)
        SS[x][2] = SS[x][2]+ Distance2(fHat, pHat, f, p)
    print ("Ok")
    


#Run on 3 sets of data
SSR(0.318, 0.763,0.75) 
SSR(0.266, 0.520,0.50)
SSR(0.681, 0.273,0.25)

print(points)


RRlist =[]
for x in range(len(SS)):
    RRlist.append(SS[x][2])
Best = (SS[RRlist.index(min(RRlist))])
print (Best)

#find probability countours
F95 = 1 + 2 / (points - 2) * scipy.stats.f.ppf(0.05, 2, points - 2)
F90 = 1 + 2 / (points - 2) * scipy.stats.f.ppf(0.1, 2, points - 2)
F70 = 1 + 2 / (points - 2) * scipy.stats.f.ppf(0.3, 2, points - 2)
F50 = 1 + 2 / (points - 2) * scipy.stats.f.ppf(0.5, 2, points - 2)


#assign value to each point given prob contour
for x in range(len(SS)):
    if SS[x][2] == Best[2]:
        SS[x][2]= Best[2] 
    elif SS[x][2]/F95 < Best[2]:
        SS[x][2]=1
    elif SS[x][2]/F90 < Best[2]:
        SS[x][2]=2
    elif SS[x][2]/F70 < Best[2]:
        SS[x][2]=3
    elif SS[x][2]/F50 < Best[2]:
        SS[x][2]=4
    else:
        SS[x][2]=5


"""Setting the range for the for loop for the plotting in the spreadsheet based on the range of reactivity ratios"""
"""Plots out the list of lists with the sum of squares with the reactivity ratios as the coordinates"""
plotMatrix = [[0 for x in range(Steps2)] for y in range(Steps1)] 
for x in range (Steps2):
    for y in range (Steps1):
        plotMatrix[x][y] +=  SS[(y)+(x*Steps1)][2]


plt.figure(figsize=(20,17))
ax = sns.heatmap(plotMatrix,xticklabels=SeaColumn, yticklabels = SeaRow, cbar=True, cmap="YlGnBu_r")
ax.set(xlabel='R1', ylabel='R2')
ax.set_xlabel('R1')    
ax.xaxis.set_label_position('top') 
ax.xaxis.tick_top()

for label in ax.xaxis.get_ticklabels()[::2]:
    label.set_visible(False)
for label in ax.yaxis.get_ticklabels()[::2]:
    label.set_visible(False)
plt.xticks(rotation=90)



plt.show(ax)
print(Best)

