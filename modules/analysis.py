import numpy as np

"***Run Length Analysis***"
def getRunLengthHistData(self, monomer):
	
	histogramData = getSingleRunLengthData(self.polymerArray, monomer)

	maximum = max(histogramData)

	binwidth = 1

	y, x = np.histogram(histogramData, bins=range(1, maximum + 1, binwidth))

	monomersUsed = getMonomersUsed(self.polymerArray)

	# normalization step
	y = y / monomersUsed

	return [x,y]

#returns an array of numbers for each consecutive monomer; to be used in histogram plotting
def getSingleRunLengthData(polymerArray, monomerID):
	#counting number of monomers in polymerArray
	#flat_list = [item for sublist in polymerArray for item in sublist]
	#print('number of monomers used: ' , len(flat_list))
	#print('indexLimit: ', indexLimit)
	#initializing array to be returned
	histogramData = []
	#count through all polymers
	for polymerObj in polymerArray: 
		#get the polymer in list form
		polymer = polymerObj.asArray()
		polymerLength = len(polymer)
		#constant to keep track of polymer index, to know when to stop counting
		polymerIndex = 1
		numConsecutive = 0
		#iterate through each polymer until monomerLimit is hit and count consecutive polymers
		for monomer in polymer:
			#if end of polymer chain reached, add any consecutives to histogram, and break from for loop
			if polymerIndex >= polymerLength:
				if  monomer == monomerID:
					numConsecutive += 1
				#print(numConsecutive)
				if numConsecutive > 0:
					count = 0
					while count < numConsecutive:
						histogramData.append(numConsecutive)
						count += 1
				break
			#if monomer is not consecutive and monomer before was, add consecutive number to data and reset values
			elif monomer != monomerID and numConsecutive > 0:
				count = 0
				while count < numConsecutive:
					histogramData.append(numConsecutive)
					count += 1
				numConsecutive = 0
				#polymerIndex += 1
				#continue
			#increment consecutive counter by 1 if monomer is consecutive
			elif  monomer == monomerID:
				numConsecutive += 1
				#polymerIndex += 1
				#continue
			polymerIndex += 1
			#print(polymerIndex)
	#print("length histogramData: ", len(histogramData))
	if not histogramData:
		histogramData = [0]
	return histogramData

def getMonomersUsed(polymerArray):
	monomersUsed = sum([polymer.len() for polymer in polymerArray])
	return monomersUsed

"***Monomer composition analysis***"
#Returns a list of the percent compostion of each monomer, in monomerID order
def getComposition(self, polymerArray):
	compositionList = [0] * self.numMonomers
	for monomerID in range(1, self.numMonomers + 1):
		for polymerObj in polymerArray:
			polymer	= polymerObj.asArray()
			for monomer in polymer:
				if monomer == monomerID:
					compositionList[monomerID - 1] += 1
	totalMonomers = sum(compositionList)
	#print("precomp: ", compositionList)
	#print("totalMonomers: ", totalMonomers)
	for monomerID in range(1, self.numMonomers + 1):
		compositionList[monomerID - 1] = int(compositionList[monomerID - 1] / totalMonomers * 1000)/10
	#print("compositionList: ", compositionList)
	return compositionList

"***DP Distribution Analysis***"
def DP_Distribution(self):
	#Get the DP Distribution
	DP_Distribution = []
	for polymerObj in self.polymerArray:
		polymer = polymerObj.asArray()
		DP_Distribution.append(len(polymer))

	#Process data into a plottable form
	maximum = max(DP_Distribution)
	minimum = min(DP_Distribution)
	binwidth = 1
	hist, bins = np.histogram(DP_Distribution, bins=range(minimum, maximum + 1, binwidth))
	x, y = bins, hist
	return [x, y]

"***Analysis for Total Mass***"
def totalMass(self):
	totalMass = 0
	for polymerObj in self.polymerArray:
		polymer = polymerObj.asArray()
		totalMass += len(polymer)
	return totalMass

"***Analysis for Number Average***"
#Given a list of polymers, returns the number average DP
def numberAverageDP(self):
	mass = totalMass(self)
	numberAverageDP = mass / self.numPolymers
	return numberAverageDP

#Given a list of polymers, return the weight average DP
"***Analysis for Weight Average***"
def weightAverageDP(self, numberAverageDP):
	weightAverageDP = 0
	for polymerObj in self.polymerArray:
		polymer = polymerObj.asArray()
		M = len(polymer)
		weightAverageDP += M*M
	weightAverageDP = weightAverageDP / numberAverageDP / self.numPolymers
	return weightAverageDP

#Given number and weight average DP, return the dispersity index
"***Analysis for Dispersity Index***"
def dispersityIndex(numberAverageDP, weightAverageDP):
	return weightAverageDP / numberAverageDP 

