from parse import parseUI_Inputs, testAssertions
from Polymer import Polymer
from graph import plotData, clearGraph
from visualize import setup_scene, draw_polymers
import random

def run_simulation(self):

	"***Set Global Variables***"
	#Order of variable setting is important!
	self.simulation_running = True

	"***Parse all relevant GUI Inputs***"
	parseUI_Inputs(self)



	"***Calculate number of polymer chains by dividing total number of monomers (pool size) by average DP***"

	#initiate a Polymer Array, which will hold all of the Polymer Objects to be simulated and created
	self.polymerArray = []

	#each monomer amount is calculated based off its ratio and the pool size
	monomerAmounts = getMonomerAmounts(self)

	#preserve a copy of the original amounts for use in Monomer Usage graph
	originalMonomerAmounts = getMonomerAmounts(self)

	#Pool size is adjusted to account for rounding errors while calculating monmer amounts
	self.adjustedPoolSize = sum(monomerAmounts)

	self.numPolymers = int(self.adjustedPoolSize / self.averageDP)

	#Updating GUI input limits, no effect on simulation
	self.rowsSpinBox.setMaximum(self.numPolymers)

	"***Check if user inputs are valid***"
	inputsValid = testAssertions(self)
	if not inputsValid:
		self.simulation_running = False
		return

	"***Set Up Visualization Scene***"
	setup_scene(self)

	"***Set up variables***"
	#list for keeping track of monomer occurrences for plotting
	self.monomerOccurrenceList = []

	#list for keeping track of monomer usage for plotting
	self.monomerRemainingList = []

	for i in range(self.numMonomers):
		self.monomerOccurrenceList.append([])
		self.monomerRemainingList.append([])

	#list for keeping count of monomer occuerences in the initiation step, to be used for calculating norm monomer occurrences
	monomerOccurrence_initiation = [0]*self.numMonomers

	#Function to append a data point, calculated by monomerOccurrence_tracker, to monomerOccurrenceList
	def update_monomerOccurrenceList(monomerOccurrence_tracker):
		for i in range(self.numMonomers):
			#how many times each monomer was used/occurred in the initiation step
			occurrence = monomerOccurrence_tracker[i]
			#normalize the number of occurrences in the initiation step
			normalizedOccurrence = float(occurrence)/sum(monomerOccurrence_tracker)
			#append this data point to monomerOccurrenceList, which will eventually be plotted
			self.monomerOccurrenceList[i].append(normalizedOccurrence)

	#Calculate a data point which represents the percentage of the monomer remaining in the pool for each monomer,
	#and add that data point to monomerRemainingList, will will eventually be plotted
	def update_monomerRemainingList(monomerAmounts, originalMonomerAmounts):
		for i in range(self.numMonomers):
			currMonomerAmount = monomerAmounts[i]
			self.monomerRemainingList[i].append(currMonomerAmount / originalMonomerAmounts[i])



	"---------------------------------------------------------------------------------------------------------------------------------"
	"***INITIATION STEP***"
	"---------------------------------------------------------------------------------------------------------------------------------"

	"In this step, all polymers will be initiated with one starting monomer. The algorithm continually intiates"
	"a chain until the chain count reaches self.numPolymers"

	for i in range(self.numPolymers):  

		"In the case where the user fixes all chains to start with a certain monomer, we simply initiate all chains"
		"with the monomer of choice."
		
		#case for setting first monomer
		"Otherwise, we will choose the initial monomer randomly based on the the mole fractions 'f' of each monomer"
		"in the feed solution. For 2- and 3- monomer systems, the instantaneous form of the Mayo Lewis Equation is "
		"used to determine the probabilities of each monomer initiating the chain."

		"Initiate a variable 'choices' that keeps track of the weight probabilty assigned to each monomer."
		#Example:

		#>>>print(choices)
		#[[1, 1.5], [2, 0.5]]

		#This means that monomer 1 has a weight of 1.5 assigned to it and monomer 2 has a weight of 0.5 assigned to it.
		#Thus, monomer 1 will initiate 3x more often than monomer 2 will.

		choices = []


		if self.numMonomers == 2 and self.model == "Mayo-Lewis":
			"Case for 2-monomer system: use the Mayo Lewis Equation"
			f1 = monomerAmounts[0]
			f2 = monomerAmounts[1]
			r1 = self.rrList[0][1]
			r2 = self.rrList[1][0]
			#print("r1: ", r1)
			#print("r2: ", r2)
			weight = (r1*f1**2 + f1*f2) / (r1*f1**2 + 2*f1*f2 + r2*f2**2)
			#print("weight: ", weight)
			choices.append([1, weight])
			choices.append([2, 1 - weight])

		

		elif self.numMonomers == 3 and self.model == "Mayo-Lewis":
			"Case for 3-monomer system: use an altered Mayo Lewis Equation"
			m1 = monomerAmounts[0]
			m2 = monomerAmounts[1]
			m3 = monomerAmounts[2]
			F = m1 + m2 + m3
			f1 = m1/F
			f2 = m2/F
			f3 = m3/F
			r11 = self.rateConstantList[0][0]
			r12 = self.rateConstantList[0][1]
			r13 = self.rateConstantList[0][2]
			r21 = self.rateConstantList[1][0]
			r22 = self.rateConstantList[1][1]
			r23 = self.rateConstantList[1][2]
			r31 = self.rateConstantList[2][0]
			r32 = self.rateConstantList[2][1]
			r33 = self.rateConstantList[2][2]
			R1 = r11 + r12 + r13
			R2 = r21 + r22 + r23
			R3 = r31 + r32 + r33
			a = f1*r11*f1/(r11*f1+r12*f2+r13*f3) + f2*r21*f1/(r21*f1+r22*f2+r23*f3) + f3*r31*f1/(r31*f1+r32*f2+r33*f3)
			b = f1*r12*f2/(r11*f1+r12*f2+r13*f3) + f2*r22*f2/(r21*f1+r22*f2+r23*f3) + f3*r32*f2/(r31*f1+r32*f2+r33*f3)
			c = 1 - a - b
			#print("startingRatioList: ", [a, b, c])
			choices.append([1,a])
			choices.append([2 ,b])
			choices.append([3,c])
		else:
			"Case for any monomer system > 3: initial solely based on feed mole ratios 'f'"

			"Also Case for any Penultimate system"

			"Cycle through each monomer, finding its initial amount and using that value as the weight."

			for i in range(self.numMonomers):
				#weight chance of monomer initation: (amount of starting monomer)
				weight = monomerAmounts[i]
				#Adds a two element list to choices containing monomer and weight
				choices.append([i+1, weight])


		"Randomly choose a monomer to initiate the polymer chain by using a weighted random selector."
		"Monomer with higher relative weights will be chosen more often. The 'weighted_choices' function"
		"takes in the 'choices' variable which has relevant weights for each monomer and runs a weighted"
		"random selection on it."

		try:
			startingMonomer = weighted_choice(choices)

		#If all weights are zero due to coefficients all being zero, then sort by relative amounts of monomer instead
		except AssertionError:
			monomerID = 1
			choices = []
			while monomerID <= self.numMonomers:
				choices.append([monomerID, monomerAmounts[monomerID - 1]])
				monomerID += 1
			startingMonomer = weighted_choice(choices)

		#Starts a new polymer with startingMonomer, represented by an array, 
		#and adds that array to polymerArray
		"A polymer is represented as a Python Class, as seen in Polymer.py. Here, we inititate an instance of the Polymer object, "
		"and append that Polymer to a list containing all of the polymers"
		polymer = Polymer([startingMonomer], self.numMonomers, self.rateConstantList, self.model)

		self.polymerArray.append(polymer)

		"Remove the monomer that was used to intitate the polymer in order to accurately update the monomer pool counts."
		"If Hold Composition is checked, then the simulation will not use up any monomers, and end when the expected number of"
		"mononomers are consumed instead"

		if not self.holdComposition:
			#Uses up one monomer of the startingMonomer
			monomerAmounts[startingMonomer - 1] -= 1





		"***For Plotting Purposes Only (no affect on simulation)***"

		#Keep track of exactly how much of each monomer is used in the total initition step by incrementing counters
		#for each monomer, represented as a list the size of self.numMonomers. Once the entire initiation step is
		#completed, a data point which represents normalized monomer occurrences can be calculated
		monomerOccurrence_initiation[startingMonomer - 1] += 1



	"***For Plotting Purposes Only (no affect on simulation)***"
	#Append a data point calculated by monomerOccurrence_initiation and append to monomerOccurrenceList
	update_monomerOccurrenceList(monomerOccurrence_initiation)
	update_monomerRemainingList(monomerAmounts, originalMonomerAmounts)


	"For Penultimate Model, the second monomer in chain must also be chosen before propagation, based off of "
	"HETEROGENOUS constants (works only for 2 monomer system)"

	if self.model == "Penultimate":
		monomerOccurrence_penultimate_initiation = [0]*self.numMonomers
		#iterating through each polymer in batch
		for polymer in self.polymerArray:
			choices = []
			#calculates weight chance for each monomer
			for i in range(self.numMonomers):
				#retrieving coefficient based on previous and currrent monomer, assumes heterogoenous penultimate monomer
				if polymer.lastMonomer() == 1:
					secondToLastMonomer = 0
				else:
					secondToLastMonomer = 1
				rr = self.rateConstantList[secondToLastMonomer-1][polymer.lastMonomer()-1][i]
				#print("coeff: ", coeff)
				# weight chance calulations for monomer attaching: coefficient*(amount of monomer remaining)
				chance = monomerAmounts[i] * rr
				#adds a two element list to choices containing monomer and weight
				choices.append([i+1, chance])
			#Using weighted_choice, selects next monomer
			nextMonomer = weighted_choice(choices)

			#Attach the monomer to the polymer
			polymer.append(nextMonomer)

			if not self.holdComposition:
				#Remove monomer from pool
				monomerAmounts[nextMonomer-1] -= 1


			"***For Plotting Purposes Only (no affect on simulation)***"

			#Calculate a data point which represents the percentage of the monomer remaining in the pool for each monomer,
			#and add that data point to monomerRemainingList, which will eventually be plotted

			update_monomerRemainingList(monomerAmounts, originalMonomerAmounts)

			#Keep track of exactly how much of each monomer is used in the total initition step by incrementing counters
			#for each monomer, represented as a list the size of self.numMonomers. Once the entire initiation step is
			#completed, a data point which represents normalized monomer occurrences can be calculated
			monomerOccurrence_penultimate_initiation[nextMonomer - 1] += 1



		"***For Plotting Purposes Only (no affect on simulation)***"
		#Append a data point calculated by monomerOccurrence_penultimate_initiation and aapend to monomerOccurrenceList
		update_monomerOccurrenceList(monomerOccurrence_penultimate_initiation)



	"***For Plotting Purposes Only (no affect on simulation)***"

	#Calculate a data point which represents the percentage of the monomer remaining in the pool for each monomer,
	#and add that data point to monomerRemainingList, which will eventually be plotted
	update_monomerRemainingList(monomerAmounts, originalMonomerAmounts)


	#print([polymer.asArray() for polymer in self.polymerArray])

	"------------------------------------------------------------------------------------------------------------------------"
	"***PROPAGATION STEP***"
	"------------------------------------------------------------------------------------------------------------------------"

	#Calculating how many monomers the reaction will use, taking into account initiated monomers
	
	if self.model == "Mayo-Lewis":
		monomersUsed = 1*self.numPolymers

	elif self.model == "Penultimate":
		monomersUsed = 2*self.numPolymers

	monomers_to_consume = int((self.percentConversion/100) * self.adjustedPoolSize) - monomersUsed

	"***Propogation for standard Mayo-Lewis Case***"


	"The 'for' loop will add monomers to the chain until we reach we use up a percentage of the total starting monomers"
	"defined by user input into the 'Percent Conversion' box. Setting 'Percent Conversion' to 100% will continue to add "
	"monomers to growing chains until there are no remaining monomers, simulating a living polymerization."
	
	monomerOccurrence_propagation = [0]*self.numMonomers

	#For plotting purposes: counter to keep track of how many monomers are added before a collective data point for 
	#monomerOccurrence is calculated and added to monomerOccurrenceList
	counter = 0


	for i in range(monomers_to_consume):

		#If process is interrupted, stop the simulation
		if self.interrupt:
			self.simulation_running = False
			run_simulation(self)
			return


		"Randomly choose a polymer chain to grow."
		polymer = random.choice(self.polymerArray)

		"For the polymer selected, iterate through all possible monomers which can be added. For each monomer, calculate"
		"the weight chance of the monomer to be added to the chain , defined as the product of the relevant rate "
		"constant 'k' times the number of monomers left unreacted 'f'"

		#A variable keeping track of monomer choices and weights
		choices = []

		for monomerID in range(1, self.numMonomers + 1):
			#Retrieveing coefficient based on previous and current monomer

			#retrieve the relevant rate constant k
			k = polymer.rateConstant(monomerID)

			# weight chance calulations for monomer attaching: coefficient*(amount of monomer remaining)
			chance = monomerAmounts[monomerID - 1] * k
			#print(monomerID, " amount of monomer remaining: ", monomerAmounts[monomerID - 1]);
			#Adds a two element list to choices containing monomer and weight
			choices.append([monomerID, chance])

		#print ("currPolymerLength: ", currPolymerLength)
		#print("choices2: ", choices)

		"Using weighted_choice, select next monomer to be appended to the growing chain"
		try:
			nextMonomer = weighted_choice(choices)

		#If all weights are zero due to coefficients all being zero, then sort by relative amounts of monomer instead
		except AssertionError:
			monomerID = 1
			choices = []
			while monomerID <= self.numMonomers:
				choices.append([monomerID, monomerAmounts[monomerID - 1]])
				monomerID += 1
			nextMonomer = weighted_choice(choices)

		"Attach next monomer to polymer chain"
		polymer.append(nextMonomer)

		"Remove the selected monomer from the pool to reflect the monomer being used up in the reaction"
		"If Hold Composition is checked, then the simulation will not use up any monomers, and end when the expected number of"
		"mononomers are consumed instead"

		if not self.holdComposition:
			monomerAmounts[nextMonomer - 1] -= 1
		#print("monomerAmounts: ", monomerAmounts)


		"***For Plotting Purposes Only (no affect on simulation)***"

		#Keep track of exactly how much of each monomer is used in the total initition step by incrementing counters
		#for each monomer, represented as a list the size of self.numMonomers. Once the batch limit (numPolymers) is
		#reached, a data point which represents normalized monomer occurrences can be calculated
		monomerOccurrence_propagation[nextMonomer - 1] += 1

		counter += 1

		#When the amount of growth cycles reaches numPolymers:

		#(1) calculate the normalized monomer occurences in the batch of data from monomerOccurrence_propagation 
		# and add the data point to monomerOccurrenceList, which will eventually be plotted
		
		#(2)Calculate a data point which represents the percentage of the monomer remaining in the pool for each monomer,
		#and add that data point to monomerRemainingList, which will eventually be plotted

		if counter == self.numPolymers:
			#print("monomerAmounts: ", monomerAmounts)

			#(1)
			update_monomerOccurrenceList(monomerOccurrence_propagation) 

			#(2)
			update_monomerRemainingList(monomerAmounts, originalMonomerAmounts)

			#Reset necessary counters and lists
			monomerOccurrence_propagation = [0]*self.numMonomers
			counter = 0
			if self.animate:
				plotData(self)
				draw_polymers(self)


	"***Plotting***"
	#print("monomerOccurrenceList: ", monomerOccurrenceList)
	if sum(monomerOccurrence_propagation) != 0:
		update_monomerOccurrenceList(monomerOccurrence_propagation) 
		update_monomerRemainingList(monomerAmounts, originalMonomerAmounts)
	clearGraph(self)
	plotData(self)
	draw_polymers(self)

	"***Update Global Variables***"
	self.simulation_running = False
	self.simulated = True
		
			

	




"***Calculates the number of each monomer given poolSize and monomerRatios***"

def getMonomerAmounts(self):
	monomerAmounts = []
	sumMonomerRatios = sum(self.monomerRatios)

	for i in range(self.numMonomers):
		monomerAmounts.append(int(self.monomerRatios[i]/sumMonomerRatios*self.poolSize))

	#print("monomerAmounts: ", monomerAmounts)
	return monomerAmounts

"***Weighted choice function***"
#Takes in a list of tuple lists with item and weight, and returns a random item based on 
#weight
def weighted_choice(choices):
    total = sum(w for c, w in choices)
    r = random.uniform(0, total)
    upto = 0
    for c, w in choices:
        if upto + w >= r:
        	return c
        upto += w
    assert False, "Shouldn't get here"
