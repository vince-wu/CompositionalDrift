from modules.parse import parseUI_Inputs, testAssertions
from modules.graph import plotData, clearGraph
from modules.visualize import setup_scene, draw_polymers
from modules.structures import Polymer, Reaction
import modules.analysis as analysis
import random
import numpy as np


def run_simulation(self):

	"***Set Global Variables***"
	#Order of variable setting is important!
	self.simulation_running = True

	"***Parse all relevant GUI Inputs***"
	parseUI_Inputs(self)

	# assign variables for creating self.reaction
	model = self.model
	num_monomer_species = self.numMonomers
	avg_dispersity = self.averageDP
	conversion = self.percentConversion/100
	chain_transfer_probability = self.chainTransferPercentage/100
	hold_composition = self.holdComposition
	rate_constants = self.rateConstantList
	init_pool_size = self.poolSize
	monomer_ratios = self.monomerRatios
	normalized_monomer_ratios = np.array(monomer_ratios)/sum(monomer_ratios)
	monomer_amounts = [int(ratio*init_pool_size) for ratio in normalized_monomer_ratios]

	# create a Reaction object
	self.reaction = Reaction(num_monomer_species, model)

	# set up reaction 
	self.reaction.set_monomer_amounts(monomer_amounts)
	self.reaction.set_rate_constants(rate_constants)
	self.reaction.set_average_DP(avg_dispersity)
	self.reaction.set_conversion(conversion)
	self.reaction.set_chain_transfer_probability(chain_transfer_probability)
	self.reaction.set_hold_composition(hold_composition)

	self.numPolymers = self.reaction.num_polymers
	self.adjustedPoolSize = self.reaction.init_pool_size
	self.polymerArray = self.reaction.polymer_list

	#Updating GUI input limits, no effect on simulation
	self.rowsSpinBox.setMaximum(self.reaction.num_polymers)

	"***Check if user inputs are valid***"
	inputsValid = testAssertions(self)
	if not inputsValid:
		self.simulation_running = False
		return
	


	"***Set Up Visualization Scene***"
	setup_scene(self)

	"***Set up variables***"
	#update statusbar - NO EFFECT ON SIMULATION
	self.statusbar.showMessage("Simulating...")
	#list for keeping track of monomer occurrences for plotting
	self.monomer_species_occurrences = []

	#list for keeping track of monomer usage for plotting
	self.monomer_species_remaining = []

	for i in range(self.numMonomers):
		self.monomer_species_occurrences.append([])
		self.monomer_species_remaining.append([])

	# list to keep track of conversion indexes
	self.conversion_index = []

	#Function to append a data point, calculated by monomers_consumed_this_step, to monomer_species_occurrences
	def update_monomer_occurences(monomers_consumed_this_step):
		# normalize the number of occurences for the step. If no monomers were consumed in this step, set value occurrences
		# equal to the previous occurences value in self.monomer_species_occurrences
		if sum(monomers_consumed_this_step) == 0:
			normalized_occurences = [self.monomer_species_occurrences[i][-1] for i in range(len(self.monomer_species_occurrences))]
		else:
			normalized_occurences = [float(occurences)/sum(monomers_consumed_this_step) for occurences in monomers_consumed_this_step]
		for i in range(len(normalized_occurences)):
			#append this data point to monomerOccurrenceList, which will eventually be plotted
			self.monomer_species_occurrences[i].append(normalized_occurences[i])

	#Calculate a data point which represents the percentage of the monomer remaining in the pool for each monomer,
	#and add that data point to monomer_species_remaining, which will eventually be plotted
	def update_monomers_remaining(monomerAmounts, originalMonomerAmounts):
		for i in range(self.numMonomers):
			currMonomerAmount = monomerAmounts[i]
			self.monomer_species_remaining[i].append(self.reaction.curr_monomer_amounts[i] / self.reaction.init_monomer_amounts[i])

	def update_conversion_index():
		num_monomers_consumed = sum(self.reaction.monomer_species_consumed)
		conversion = float(num_monomers_consumed)/self.reaction.init_pool_size*100
		self.conversion_index.append(conversion)



	"---------------------------------------------------------------------------------------------------------------------------------"
	"***INITIATION STEP***"
	"---------------------------------------------------------------------------------------------------------------------------------"

	"In this step, all polymers will be initiated with one starting monomer. The algorithm continually intiates"
	"a chain until the chain count reaches self.numPolymers"

	self.reaction.run_initiatation()
	monomers_consumed_this_step = np.array(self.reaction.monomer_species_consumed).copy()
	"***For Plotting Purposes Only (no affect on simulation)***"
	#Append a data point calculated by monomerOccurrence_initiation and append to monomer_species_occurrences
	update_monomer_occurences(monomers_consumed_this_step)
	update_monomers_remaining(self.reaction.curr_monomer_amounts, self.reaction.init_monomer_amounts)
	update_conversion_index()

	"------------------------------------------------------------------------------------------------------------------------"
	"***PROPAGATION STEP***"
	"------------------------------------------------------------------------------------------------------------------------"

	#For plotting purposes: counter to keep track of how many monomers are added, used to determine when plots and visualizations are updated
	counter = 0
	monomer_occurence_tracker = [0]*self.numMonomers
	monomers_consumed_prior = np.array(self.reaction.monomer_species_consumed).copy()
	while sum(self.reaction.monomer_species_consumed) < self.reaction.init_pool_size*self.reaction.conversion:
		# print(counter)
		self.reaction.run_single_propagation()
		"***For Plotting Purposes Only (no affect on simulation)***"

		#Keep track of exactly how much of each monomer is used in the total initition step by incrementing counters
		#for each monomer, represented as a list the size of self.numMonomers. Once the batch limit (numPolymers) is
		#reached, a data point which represents normalized monomer occurrences can be calculated

		counter += 1
		# for i in range(len(monomers_consumed)):
		# 	monomer_occurence_tracker[i] += monomers_consumed[i]

		#When the amount of growth cycles reaches numPolymers:

		#(1) calculate the normalized monomer occurences in the batch of data from monomerOccurrence_propagation 
		# and add the data point to monomer_species_occurrences, which will eventually be plotted
		
		#(2)Calculate a data point which represents the percentage of the monomer remaining in the pool for each monomer,
		#and add that data point to monomer_species_remaining, which will eventually be plotted

		#(3) Update the conversion index

		if counter == self.reaction.num_polymers:

			#(1)
			# amount of monomers consumed in this step is equal to total monomer species consumed subtracted by monomers consumed
			# in the previous step. Note that "monomers_consumed_this_step" is an array keeping track of monomer consumption for
			# each monomer species
			monomers_consumed_this_step = np.array(self.reaction.monomer_species_consumed) - monomers_consumed_prior
			update_monomer_occurences(monomers_consumed_this_step) 
			#(2)
			update_monomers_remaining(self.reaction.curr_monomer_amounts, self.reaction.init_monomer_amounts)
			# (3)
			update_conversion_index()
			#Reset necessary counters and lists
			counter = 0
			monomers_consumed_prior = np.array(self.reaction.monomer_species_consumed).copy()
			# plot data and update visualizations
			if self.animate:
				plotData(self)
				draw_polymers(self)

	"***Plotting***"
	#print("monomer_species_occurrences: ", monomer_species_occurrences)
	monomers_consumed_this_step = np.array(self.reaction.monomer_species_consumed) - monomers_consumed_prior
	update_monomer_occurences(monomers_consumed_this_step)
	update_monomers_remaining(self.reaction.curr_monomer_amounts, self.reaction.init_monomer_amounts)
	update_conversion_index()
	
	# print("monomer 0: {}".format(len(self.monomer_species_remaining[0])))
	# print(self.monomer_species_remaining[0])
	# print("monomer 1: {}".format(len(self.monomer_species_remaining[0])))
	# print(self.monomer_species_remaining[1])
	clearGraph(self)
	plotData(self)
	draw_polymers(self)

	"***Update Global Variables***"
	lambdaValue = analysis.calculate_theta(self)
	self.simulation_running = False
	self.simulated = True
	self.statusbar.showMessage("Done.", 3000)
			

	




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
