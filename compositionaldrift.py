#All neccesary imports
import matplotlib
import random
import math
import json
import os.path
matplotlib.use('TkAgg')
from numpy import arange, sin, pi
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
# implement the default mpl key bindings
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import sys
if sys.version_info[0] < 3:
    import Tkinter as Tk
    import Tkinter.filedialog
else:
    import tkinter as Tk
    import tkinter.filedialog

"""
    Created by Vincent Wu on 8/17/16:
    This program uses the Mayo-Lewis equation and Monte Carlo method to simulate copolymer growth, and
    represents results both graphically and visually
    """
#Static Final variables
MONOMER_CAP = 5000000
NUM_UNIQUE_MONOMERS = 2
NUM_SIMULATIONS = 200
NUM_POLYMERS_SHOW = 8
GRAPH_OCCURENCE = 1
TOTAL_STARTING_MONOMERS = 1000
SETTING = 0
RAFT_RATIO = 0.01
LOAD_SUCCESSFUL = False
CONFIGS = [["Number of Unique Monomers", 1], ["Number of Simulations", 1],
 ["Number of Polymers to Show", 1], 
 ["Graph Monomer Occurence", 1], ["Total Starting Monomers", 1], ["RAFT to Monomers Ratio", 0], 
 ["Default Setting", 1], ["Monomer Cap", 1]]
#generates a config file if needed
def generateConfigFile():
	if not os.path.exists("config.txt"):
		file = open("config.txt", "w")
		file.write("Number of Unique Monomers = 2 \nNumber of Simulations = 200 \nNumber of Polymers to Show = 8 \nGraph Monomer Occurence = 1 \n")
		file.write("Total Starting Monomers = 1000 \nRAFT to Monomers Ratio = 0.01 \nDefault Setting = 0 \nMonomer Cap = 5000000 \n")
		file.write("Setting 1 \nMonomer 1 Ratio = 50 \nMonomer 2 Ratio = 25 \nMonomer 3 Ratio = 20 \nMonomer 4 Ratio = 5 \n") 
		file.write("1-1 = 0.89 \n1-2 = 1 \n1-3 = 1 \n1-4 = 1 \n2-1 = 1 \n2-2 = 1.1 \n2-3 = 1.1 \n2-4 = 1.1 \n3-1 = 1 \n3-2 = 1.1 \n3-3 = 1.1 \n3-4 = 1.1 \n")
		file.write("4-1 = 1 \n4-2 = 1.1 \n4-3 = 1.1 \n4-4 = 1.1 \nend")
		file.close()
#reads the config file and sets static variables based on config
def readConfigFile():
	#error checking to see if config file exists
	if not os.path.exists("config.txt"):
		errorMessage("config.txt does not exist!", 220)
		return
	file = open("config.txt", "r")
	#parses each line in the file. If it is a valid configuration line, it sets that correspondin static variable accordingly
	#variable to keep track of number of invalid config lines
	invalidLines = 0 
	#variable to keep track of whether or not line is a setting config line
	settingConfig = False
	#variable to keep track of whether to parse setting config
	parseSetting = False
	#variable to keep track of number of monomers
	global numMonomers
	numMonomers = 0
	#global variable to shwo whether or not a setting loaded successfully
	global LOAD_SUCCESSFUL
	for line in file:
		line = line.strip()
		#If line is end, read next lines as regular config line
		if line == "end":
			if settingConfig == True and parseSetting == True:
				#test to see that the settings are all valid and can be loaded properly
				try:
					#test to see that ratios for all monomers are given
					for monomerID in numMonomerArray:
						assert(numMonomerArray[monomerID -1] == monomerID)
					#test to see that all ratios are valid (> 0)
					for ratio in RATIO_ARRAY:
						assert(ratio > 0)
					for column in COEFF_ARRAY:
						for coeff in column:
							assert(coeff >= 0)
				except AssertionError:
					parseSetting = False
					LOAD_SUCCESSFUL = False
					continue
				except IndexError:
					parseSetting = False
					LOAD_SUCCESSFUL = False
				LOAD_SUCCESSFUL = True
			settingConfig = False
			parseSetting = False
			continue
		#handling for setting config line
		if settingConfig:
			#if not relevant setting number, skip
			if not parseSetting:
				continue
			lineArray = line.split(" ")
			#first line of setting config must state number of monomers
			if len(lineArray) == 6:
				try:
					assert(lineArray[0] == "Number")
					assert(lineArray[1] == "of")
					assert(lineArray[2] == "Unique")
					assert(lineArray[3] == "Monomers")
					assert(lineArray[4] == "=")
					numMonomers = int(lineArray[5])
					#global array to keep track of monomer ratios, initalized as an array of -1
					global RATIO_ARRAY
					RATIO_ARRAY = [-1] * numMonomers
					#array to keep track of all neccesary monomers
					numMonomerArray = [0] * numMonomers
					#global array to keep track monomer coefficients, initialized of an array of arrays of -1
					global COEFF_ARRAY
					COEFF_ARRAY = [-1] * numMonomers
					COEFF_ARRAY = [[-1 for i in range(numMonomers)] for j in range(numMonomers)]
				except AssertionError:
					invalidLines += 1
				except ValueError:
					invalidLines += 1
				continue
			if len(lineArray) == 5:
				#checking that line is correct syntax
				try:
					assert(lineArray[0] == "Monomer")
					assert(lineArray[2] == "Ratio")
					assert(lineArray[3] == "=")
					#adding monomer number to numMonomerArray
					numMonomerArray[int(lineArray[1]) - 1] = int(lineArray[1])
					#add monomer ratio to RATIO_ARRAY
					RATIO_ARRAY[int(lineArray[1]) - 1] = int(lineArray[4])
				except AssertionError:
					invalidLines += 1
					continue
				except ValueError:
					invalidLines += 1
					continue 
				except IndexError:
					invalidLines += 1
					continue
			if len(lineArray) == 3:
				#add the coefficient to the correct place in COEFF_ARRAY
				try:
					assert(len(lineArray[0]) == 3)
					assert(lineArray[1] == "=")
					coeffTag = lineArray[0]
					assert(coeffTag[1] == "-")
					column = int(coeffTag[0]) - 1
					row = int(coeffTag[2]) - 1
					COEFF_ARRAY[column][row] = float(lineArray[2])
				except AssertionError:
					invalidLines += 1
					continue
				except ValueError:
					invalidLines += 1
					continue
				except IndexError:
					invalidLines += 1
					continue
		#If line is a "Setting X", read next lines as setting config line
		intermediate = line.split("=")
		if len(intermediate) == 1:
			intermediate2 = intermediate[0].split(" ")
			#cehck if line is two words
			if len(intermediate2) == 2:
				#check if line starts with "Setting"
				if intermediate2[0] == "Setting":
					settingConfig = True
				#check if setting number is integer and if it matches SETTING
				try:
					if int(intermediate2[1]) == SETTING:
						parseSetting = True
				except AssertionError:
					invalidLines += 1
				continue
			invalidLines += 1
			continue
		#gets the config header
		configType = (line.split("="))[0].strip()
		#gets the config value
		configStringValue = (line.split("="))[1].strip()
		if (not setConfigVariable(configType, configStringValue)):
			invalidLines += 1
	file.close()
	print("Number of invalid config lines: ", invalidLines)
#helper method to set the static variable. Returns 1 if successful, 0 if not
def setConfigVariable(configType, configStringValue):
	for config in CONFIGS:
		if config[0] == configType:
			try:
				#handling for variable types
				if config[1]:
					setConfigVariableHelper(configType, int(configStringValue))
				else:
					setConfigVariableHelper(configType, float(configStringValue))
				print(configType, configStringValue)
				return 1
			except AssertionError:
				return 0
	return 0
#A helper function mapping name to static variables
def setConfigVariableHelper(configType, configValue):
	if configType == "Number of Unique Monomers":
		global NUM_UNIQUE_MONOMERS
		NUM_UNIQUE_MONOMERS = configValue
	elif configType == "Number of Simulations":
		global NUM_SIMULATIONS
		NUM_SIMULATIONS = configValue
	elif configType == "Number of Polymers to Show":
		global NUM_POLYMERS_SHOW
		NUM_POLYMERS_SHOW = configValue	
	elif configType == "Graph Monomer Occurence":
		global GRAPH_OCCURENCE
		GRAPH_OCCURENCE = configValue
	elif configType == "Total Starting Monomers":
		global TOTAL_STARTING_MONOMERS
		TOTAL_STARTING_MONOMERS = configValue
	elif configType == "RAFT to Monomers Ratio":
		global RAFT_RATIO
		RAFT_RATIO = configValue
	elif configType == "Default Setting":
		global SETTING
		SETTING = configValue
	elif configType == "Monomer Cap":
		global MONOMER_CAP
		MONOMER_CAP = configValue
	else:
		assert False, "shouldn't get here"	
#Main class 
class Application(Tk.Frame):
	def __init__(self, master = None):
		Tk.Frame.__init__(self, master)
		self.pack()
		self.initialize()
	#Initialization
	def initialize(self):
		#generates and reads config file
		generateConfigFile()
		readConfigFile()
		if LOAD_SUCCESSFUL:
			print("Load successful!")
		#Creates the init screen
		self.initScreen()
		#Creates Input Widgets
		self.createInputWidgets()
		#Creates dummy visualization Frame
		self.visualizationFrame = Tk.Frame(master = root)
		self.destroyCanvas = False
		self.destroyHide = False
	#Destroys unneccesary widgets
	def destroyWidgets(self):
		self.parentFrame.destroy()
		self.visualizationFrame.destroy()
		if self.destroyCanvas:
			self.canvas.get_tk_widget().destroy()
			self.toolbar.destroy()
		self.destroyCanvas = False
		self.destroyHide = False
	#Creates input widgets
	def createInputWidgets(self):
		def loadSettings(self):
			global useLoadedSettings
			useLoadedSettings = True
			self.createMoreInputs()
		#Confirms number of monomers, creates more input widgets
		def enter(self):
			#nonlocal numMonomers 
			self.numMonomers = int(self.monomerCountTkVar.get())
			#asserting that input is in correct range
			assert self.numMonomers < 8
			global useLoadedSettings
			useLoadedSettings = False
			#print(numMonomers) # debugging purposes
			self.createMoreInputs() 
		#Clears unneccesary widgets and creates more input widgets; called by enter()
		def _quit():
			root.quit()
			root.destroy()
		#PlaceHolder parent inputFrame
		self.parentFrame = Tk.Frame(master = root)
		self.parentFrame.pack(side = Tk.BOTTOM, fill = Tk.X, expand = 0, pady = 3)
		#The parent LabelFrame for all Input Widgets
		self.inputFrame = Tk.LabelFrame(master = self.parentFrame, text = "Input Parameters")
		self.inputFrame.pack(side = Tk.TOP, fill = Tk.X, expand = 0, padx = 3, pady = 5)
		#Frame for row one inputs
		self.rowFrame1 = Tk.Frame(master = self.inputFrame)
		self.rowFrame1.pack(side = Tk.TOP, fill = Tk.X, expand = 0, padx = 3, pady = 5)
		#A frame for all buttons in inputFrame
		self.buttonFrame = Tk.Frame(master = self.inputFrame)
		self.buttonFrame.pack(side = Tk.LEFT, padx = 5)
		#A frame for spinbox and label
		countFrame = Tk.Frame(master = self.buttonFrame)
		countFrame.pack(side = Tk.TOP, padx = 5, pady = 0)
		#Label for spinbox
		self.monomerCountLabel = Tk.Label(master = countFrame, text = "Number of Unique Monomers:")
		self.monomerCountLabel.pack(side = Tk.LEFT, padx = 0, pady = 0)
		#tkvar for monomercount
		self.monomerCountTkVar = Tk.IntVar()
		#MonomerCount spinbox
		self.monomerCountSpinbox = Tk.Spinbox(master = countFrame, from_ = 1, to = 7, width = 2, textvariable = self.monomerCountTkVar)
		self.monomerCountSpinbox.pack(side = Tk.LEFT, padx = 5, pady = 0)
		#sets monomerCOuntTkVar to default setting
		self.monomerCountTkVar.set(NUM_UNIQUE_MONOMERS)
		#countConfirm Button
		self.countConfirmButton = Tk.Button(master = self.inputFrame, text = "Enter",
		command = lambda:enter(self), bg = "light blue", activebackground = "light slate blue", width = 9)
		self.countConfirmButton.pack(side = Tk.LEFT, padx = 5, pady = 5)
		#frame for loaded input option
		#self.rowFrame2 = Tk.Frame(master = self.inputFrame)
		#self.rowFrame2.pack(side = Tk.TOP, fill = Tk.X, expand = 0, padx = 3, pady = 5)
		#use loaded inputs button
		self.loadButton = Tk.Button(master = self.inputFrame, text = "Load from Settings",
		command = lambda:loadSettings(self), bg = "light blue", activebackground = "light slate blue", width = 15)
		self.loadButton.pack(side = Tk.LEFT, padx = 5, pady = 5)
	#Creates more input widgets based on numMonomers
	def createMoreInputs(self):
		#Destroys or edits current widgets
		#case for loading inputs
		if LOAD_SUCCESSFUL and useLoadedSettings:
			self.numMonomers = numMonomers
		self.initFrame.destroy()
		self.monomerCountSpinbox.destroy()
		self.monomerCountLabel.config(text = "Number of Unique Monomers: " + str(self.numMonomers))
		self.countConfirmButton.destroy()
		self.loadButton.destroy()
		#Commands
		# Quit command: quits window
		def _quit():
			root.quit()
			root.destroy()
		# Back Command: goes back to numMonomers Entry
		def back(self):
			#Destroys all neccesary widgets
			self.destroyWidgets()
			#Re-Initiates
			self.initScreen()
			self.createInputWidgets()
		#Frame for numSimulations label and Entry
		self.numSimsFrame = Tk.Frame(master = self.buttonFrame)
		self.numSimsFrame.pack(side = Tk.TOP)
		#Label for numSims Entry
		self.numSimsLabel = Tk.Label(master = self.numSimsFrame, text = "Number of Simulations:")
		self.numSimsLabel.pack(side = Tk.LEFT, pady = 3)
		#Entry for numSimulations
		self.numSimsEntry = Tk.Entry(master = self.numSimsFrame, width = 5)
		self.numSimsEntry.pack(side = Tk.LEFT, padx = 3, pady = 3)
		#Setting number of simulations to 1000
		self.numSimsTkVar = Tk.IntVar()
		self.numSimsEntry["textvariable"] = self.numSimsTkVar
		self.numSimsTkVar.set(NUM_SIMULATIONS)
		#Frame for numPolyToShow SpinBox
		self.numPolyToShowFrame = Tk.Frame(master = self.buttonFrame)
		self.numPolyToShowFrame.pack(side = Tk.TOP)
		#Label for numPolyToShow spinbox
		self.numPolyToShowLabel = Tk.Label(master = self.numPolyToShowFrame, text = "Number of Polymers to Show:")
		self.numPolyToShowLabel.pack(side = Tk.LEFT, padx = 0, pady = 0)
		#numPolyToShow spinbox
		self.numPolyToShowBox = Tk.Spinbox(master = self.numPolyToShowFrame, from_ = 1, to = 12, width = 2)
		self.numPolyToShowBox.pack(side = Tk.LEFT, padx = 5, pady = 0)
		#Setting numPolyToShow to 6
		self.numPolyToShow = Tk.IntVar()
		self.numPolyToShowBox["textvariable"] = self.numPolyToShow
		self.numPolyToShow.set(NUM_POLYMERS_SHOW)
		#Frame for Back, Quit, and Simulate buttons
		self.backSimFrame = Tk.Frame(master = self.buttonFrame)
		self.backSimFrame.pack(side = Tk.TOP)
		#A simulate button to simulate polymer formation
		self.simulateButton = Tk.Button(master = self.backSimFrame, text = "Simulate", width = 7,
		 command = self.simulate, bg = "light blue", activebackground = "light slate blue")
		self.simulateButton.pack(side = Tk.LEFT, padx = 6, pady = 4)
		#A back button to enter a diff number of monomers
		self.backButton = Tk.Button(master = self.backSimFrame, text = "Back", width = 7,
		 command = lambda:back(self), bg = "light blue", activebackground = "light slate blue")
		self.backButton.pack(side = Tk.LEFT, padx = 6, pady = 4)	
		#Quit Button Widget
		quitButton = Tk.Button(master = self.backSimFrame, text = "Quit",
		 command = _quit, width = 7, bg = "light blue", activebackground = "light slate blue")
		quitButton.pack(side = Tk.LEFT, padx = 6)
		createCount = 0;
		# A list of Tk.Entry objects for Monomer ratios
		self.startingRatiosTkList = [] 
		# A 2D list of Tk.Entry objects for Coefficicients
		self.coefficientList = [] 
		#Frame for total number of monomers and RAFT to monomer ratio
		self.initialConditionsFrame = Tk.Frame(master = self.inputFrame)
		self.initialConditionsFrame.pack(side = Tk.LEFT, padx = 5)
		#variable to keep track of type of graph. 0 = percentage, 1 = occurences
		self.graphTypeTkIntVar = Tk.IntVar()
		self.graphTypeTkIntVar.set(GRAPH_OCCURENCE)
		#RadioButtons for display type
		self.percentageRadioButton = Tk.Radiobutton(master = self.initialConditionsFrame, text = "Graph Monomer Occurences",
		 variable = self.graphTypeTkIntVar, value = 1)
		self.occurenceRadioButton = Tk.Radiobutton(master = self.initialConditionsFrame, text = "Graph Percentage Monomer",
		 variable = self.graphTypeTkIntVar, value = 0)
		self.percentageRadioButton.pack(side = Tk.TOP, anchor = Tk.W)
		self.occurenceRadioButton.pack(side = Tk.TOP, anchor = Tk.W)
		#Frame for totalMonomers label and Entry
		self.totalMonomersFrame = Tk.Frame(master = self.initialConditionsFrame)
		self.totalMonomersFrame.pack(side = Tk.TOP)
		#Label for totalMonomers Entry
		self.totalMonomersLabel = Tk.Label(master = self.totalMonomersFrame, text = "Total Starting Monomers:")
		self.totalMonomersLabel.pack(side = Tk.LEFT, pady = 3)
		#Entry for totalMonomers
		self.totalMonomersEntry = Tk.Entry(master = self.totalMonomersFrame, width = 5)
		self.totalMonomersEntry.pack(side = Tk.LEFT, padx = 3, pady = 3)
		#Setting totalMonomers to 1000
		self.totalMonomersTkVar = Tk.IntVar()
		self.totalMonomersEntry["textvariable"] = self.totalMonomersTkVar
		self.totalMonomersTkVar.set(TOTAL_STARTING_MONOMERS)
		#Frame for raftRatio label and Entry
		self.raftRatioFrame = Tk.Frame(master = self.initialConditionsFrame)
		self.raftRatioFrame.pack(side = Tk.TOP)
		#Label for raftRatio Entry
		self.raftRatioLabel = Tk.Label(master = self.raftRatioFrame, text = "RAFT to Monomers Ratio:")
		self.raftRatioLabel.pack(side = Tk.LEFT, pady = 3)
		#Entry for raftRatio
		self.raftRatioEntry = Tk.Entry(master = self.raftRatioFrame, width = 5)
		self.raftRatioEntry.pack(side = Tk.LEFT, padx = 3, pady = 3)
		#Setting number of simulations to 1000
		self.raftRatioTkVar = Tk.StringVar()
		self.raftRatioEntry["textvariable"] = self.raftRatioTkVar
		self.raftRatioTkVar.set(RAFT_RATIO)
		#Frame for Monomer Amounts
		self.amountFrame = Tk.Frame(master = self.inputFrame) 
		self.amountFrame.pack(side = Tk.LEFT, padx = 5)
		#Frame for Monomer Coefficients
		self.coefficientFrame = Tk.Frame(master = self.inputFrame)
		self.coefficientFrame.pack(side = Tk.LEFT, padx = 5)
		#While loop creating number of neccesary amount Entry boxes
		while createCount < self.numMonomers:
			#Label for inputAmount
			monomerAmountFrame = Tk.Frame(master = self.amountFrame)
			monomerAmountFrame.pack(side = Tk.TOP, padx = 5, pady = 3)
			inputAmountLabel = Tk.Label(master = monomerAmountFrame, text = "     Monomer " 
				+ str(createCount + 1) + " Ratio:")
			inputAmountLabel.pack(side = Tk.LEFT)
			#Entry for inputAmount
			inputAmount = Tk.Entry(master = monomerAmountFrame, width = 5)
			inputAmount.pack(side = Tk.LEFT, padx = 5)
			#Setting Default Value to 20
			amount = Tk.IntVar()
			inputAmount["textvariable"] = amount
			if LOAD_SUCCESSFUL and useLoadedSettings:
				amount.set(RATIO_ARRAY[createCount])
			else:
				amount.set(20)
			#Add Tk.Entry object to startingAmountList
			self.startingRatiosTkList.append(inputAmount)
			createCount += 1
		#Debugging purposes
		#print("startingAmountList: ", self.startingRatiosTkList) 
		createCount2 = 0
		#While loop creating number of neccesary coefficient Entry boxes
		while createCount2 < self.numMonomers:
			combinations = 0
			#Appends to coefficient list a list containing coefficients for the polymer index
			singleMonoCoeffList = []
			self.coefficientList.append(singleMonoCoeffList)
			#Frame for Coefficients for Single Monomer
			singleCoeffFrame = Tk.Frame(master = self.coefficientFrame)
			singleCoeffFrame.pack(side = Tk.LEFT, fill = Tk.X, expand = 1)
			while combinations < self.numMonomers:				
				#Label for inputAmount
				coeffValFrame = Tk.Frame(master = singleCoeffFrame)
				coeffValFrame.pack(side = Tk.TOP, padx = 5, pady = 3)
				inputCoeffLabel = Tk.Label(master = coeffValFrame, text = str(createCount2 + 1)
				 + "-" + str(combinations + 1) + " Constant:" )
				inputCoeffLabel.pack(side = Tk.LEFT)
				#Entry for inputAmount
				inputCoeff = Tk.Entry(master = coeffValFrame, width = 4)
				inputCoeff.pack(side = Tk.LEFT, padx = 5)
				#Setting Default Coefficient to 1
				coeff = Tk.IntVar()
				inputCoeff["textvariable"] = coeff
				if LOAD_SUCCESSFUL and useLoadedSettings:
					coeff.set(COEFF_ARRAY[createCount2][combinations])
				else:
					coeff.set(1)
				#Add a Tk.Entry object to singleMonoCoeffList
				singleMonoCoeffList.append(inputCoeff)
				combinations += 1
			createCount2 += 1
		#Debugging Purposes: prints the 2D array which should store coefficients 	
		#print("coefficientList: " , self.coefficientList)
		#counter = 0
		"""for coefflist in self.coefficientList:
			print("list at index " + str(counter) + ": ", coefflist)
			counter += 1"""
	#Frame and contents for opening screen
	def initScreen(self):
		self.initFrame = Tk.Frame(master = root)
		self.initFrame.pack(side = Tk.TOP, fill = Tk.X, expand = 0, padx = 3, pady = 5)
		message1 = Tk.Message(master = self.initFrame, width = 250, font = ("Times New Roman", 10, "bold"),
		 text = "Compositional Drift Simulator v1.2")
		message1.pack(side = Tk.TOP, padx = 0, pady = 0)
		message2 = Tk.Message(master = self.initFrame, width = 200,
			text = "Author: Vincent Wu")
		message2.pack(side = Tk.TOP)
		message3 = Tk.Message(master = self.initFrame, width = 300,
			text = "Welcome to Compositional Drift Simulator! This program uses the Mayo Lewis Equation and the" 
			" Monte Carlo method to simulate the growth of a copolymer chain."
			" Please input your desired number of unique copolymers and press 'Enter' to continue.")
		message3.pack(side = Tk.TOP)
	#Hides input params
	def showInputParams(self):
		self.inputFrame.pack(side = Tk.BOTTOM, pady = 3)
		self.showButton.destroy()
	#Shows input params
	def hideInputParams(self):
		self.inputFrame.pack_forget()
		self.showButton = Tk.Button(master = self.parentFrame, text = "Show Input Parameters", width = 27,
			bg = "pale turquoise", activebackground = "light slate blue", command = self.showInputParams)
		self.showButton.pack(side = Tk.TOP)
	#Simulates polymer reaction based on input values
	def simulate(self):
		#print("Simulating!")
		#Asserts that there are no input errors. Shows errorMessage if error caught.test
		try:
			self.totalMonomers = int(self.totalMonomersTkVar.get())
			monomerAmounts = self.getMonomerAmounts()
			singleCoeffList = self.getCoefficients()
			self.numSimulations = int(self.numSimsTkVar.get())
			self.raftRatio = float(self.raftRatioTkVar.get())
			assert(self.totalMonomers > 0)
			assert(self.numSimulations > 0)
			assert(self.raftRatio > 0)
			assert(self.totalMonomers * self.numSimulations <= MONOMER_CAP)
		except ValueError:
			errorMessage("Please input valid parameters!", 220)
			return
		except notInEuropeError:
			errorMessage("You are not in Europe!", 220)
			return
		except AssertionError:
			errorMessage("Please input valid parameters!", 220)
			return
		"""Calculates Initial Conditions based on inputs"""
		#number of monomers in each polymer assuming reaction goes to completion, based on raftRatio and itotalMonomers
		print(self.totalMonomers)
		print(self.raftRatio)
		self.polymerLength = int(1 / self.raftRatio)
		if self.polymerLength == 0:
			self.polymerLength = 1
		self.numPolymers = int(self.totalMonomers / self.polymerLength)
		#Asserting valid inputs for RAFT ration and totalMonomers
		try:
			assert self.numPolymers > 0
		except AssertionError:
			errorMessage("RAFT Ratio too small!", 220)
			return
		#Debugging Purposes
		print("Polymer Length: ", self.polymerLength)
		print("Number of Polymers: ", self.numPolymers)
		#destroys hideButton if necessary
		if self.destroyHide:
			self.hideButton.destroy()
		#hideButton creation, _DEPRECATED_
		"""self.hideButton = Tk.Button(master = self.buttonFrame, bg = "pale turquoise", activebackground = "light slate blue",
			width = 27, text = "Hide Input Paramters", command = self.hideInputParams)
		self.hideButton.pack(side = Tk.TOP, pady = 3)
		self.destroyHide = True"""
		#destroys canvas if necessary
		if self.destroyCanvas:
			self.canvas.get_tk_widget().destroy()
			self.toolbar.destroy()
		self.destroyCanvas = True
		#print("monomerAmounts: ", monomerAmounts)
		#print("singleCoeffList: ", singleCoeffList)
		self.monomerCountLabel.config(text = "Polymers Per Simulation: " + str(self.numPolymers))
		#An array of polymers
		polymerArray = []
		#keeps track of number of simulations
		simCounter = 1
		self.originalMonomerAmounts = list(monomerAmounts)
		while simCounter <= self.numSimulations:
			#a local array of polymers
			localPolymerArray = []
			#sets monomerAmounts to orignalMonomerAmounts
			monomerAmounts = list(self.originalMonomerAmounts)
			#print("originalMonomerAmounts: ", originalMonomerAmounts)
			#variable keeping track of current number of polymers
			currNumPolymers = 0
			#print(self.numMonomers)
			#Builds number of polymers equal to self.numpolymers
			while currNumPolymers < self.numPolymers:
				"""Initiation: Chooses a monomer to initiate the chain based of weighted chance"""
				#A variable keeping track of current monomer
				monomerID = 1
				choices = []
				while monomerID <= self.numMonomers:
					#weight chance of monomer initation: (amount of starting monomer)
					weight = monomerAmounts[monomerID - 1]
					#Adds a two element list to choices containing monomer and weight
					choices.append([monomerID, weight])
					monomerID += 1
				#print("choices: ", choices)
				#Using weighted_choice, selects first monomer
				startingMonomer = weighted_choice(choices)
				#Starts a new polymer with startingMonomer, represented by an array, 
				#and adds that array to polymerArray
				localPolymerArray.append([startingMonomer])
				#Uses up one monomer of the startingMonomer
				monomerAmounts[startingMonomer - 1] -= 1
				#increases number of polymers by 1
				currNumPolymers += 1
			#debugging starting monomers
			#print("polymer array", polymerArray)
			#variable to keep track of polymer length
			currPolymerLength = 1
			"""Grows each polymer at the same time until they all reach desired polymer size"""
			#print("self.polymerLength: ", self.polymerLength)
			while currPolymerLength < self.polymerLength:
				"""RAFT polymerization: each polymer propagates at the same time, represented here as
				monomers being attached to each polymer chain one at a time"""
				for polymer in localPolymerArray:
					#A variable keeping track of current monomer
					monomerID = 1
					choices = []
					"""Propogation: Iterates through all monomers, calculating weight chance to bind for each, 
					and binds monomer with highest weight chance"""
					while monomerID <= self.numMonomers:
						#Retrieveing coefficient based on previous and current monomer
						coeff = singleCoeffList[polymer[-1] - 1][monomerID - 1]
						#print("coeff: ", coeff);
						# weight chance calulations for monomer attaching: coefficient*(amount of monomer remaining)
						chance = monomerAmounts[monomerID - 1] * coeff
						#print(monomerID, " amount of monomer remaining: ", monomerAmounts[monomerID - 1]);
						#Adds a two element list to choices containing monomer and weight
						choices.append([monomerID, chance])
						monomerID += 1
					#print ("currPolymerLength: ", currPolymerLength)
					#print("choices2: ", choices)
					#Using weighted_choice, selects next monomer
					try:
						nextMonomer = weighted_choice(choices)
					#If all weights are zero due to coefficients, then sort by relative amounts of monomer instead
					except AssertionError:
						monomerID = 1
						choices = []
						while monomerID <= self.numMonomers:
							choices.append([monomerID, monomerAmounts[monomerID - 1]])
							monomerID += 1
						nextMonomer = weighted_choice(choices)
					#Reduces number of nextMonomer by 1, since it is being used up in reaction
					monomerAmounts[nextMonomer - 1] -= 1
					#print("monomerAmounts: ", monomerAmounts)
					#Attaches next monomer to polymer chain
					polymer.append(nextMonomer)
				#increase current polymer length by 1
				currPolymerLength += 1
			#print("simCounter: ", simCounter)
			#adds local polymer array to global polymer array
			polymerArray = polymerArray + localPolymerArray
			#increases simCounter by 1
			simCounter += 1
		#Debugging purposes
		"""Important Debug: Prints out array of polymers"""
		#print("Array of Polymers: ", polymerArray)
		#outputs polymer onto textfile
		text_file = open("polymerArray.txt", "w")
		json.dump(polymerArray, text_file)
		text_file.close()
		self.visualizationFrame.destroy()
		self.visualizePolymers(polymerArray)
		self.plotCompositions(polymerArray)
		center(root)
		#self.inputFrame.pack_forget()
	#plots compositions given a PolymerArray
	def plotCompositions(self, polymerArray):
		self.lineColors = []
		#Plot and Figure formatting
		self.plotFigure = Figure(figsize=(5.5, 3.3), dpi=100)
		frequencyPlot = self.plotFigure.add_subplot(111)
		frequencyPlot.tick_params(labelsize = 7)
		frequencyPlot.set_xlabel("Monomer Position Index", labelpad = 0, fontsize = 9)
		#Iterates through each unique monomer and plots composition
		for monomer in range(1, self.numMonomers + 1):
			#x-axis array
			polymerIndex = list(range(1, self.polymerLength + 1))
			#y-axis array initation
			monomercounts = [0] * self.polymerLength
			#inputs counts into y-axis array
			graphType = self.graphTypeTkIntVar.get()
			#graphs Percentage of Monomer Remaining
			if graphType == 0:
				#adjust axis title
				frequencyPlot.set_ylabel("Percentage of Monomer Remaining", labelpad=5, fontsize = 9)
				#adjust y axis limiys
				frequencyPlot.set_ylim([0,1])
				#variable to keep track of average number of monomers consumed
				monomersConsumed = 0
				for index in polymerIndex:
					count = 0
					for polymer in polymerArray:
						if polymer[index - 1] == monomer:
							count += 1
					startingMonomerAmount = self.originalMonomerAmounts[monomer - 1]
					#calculates monomer consumed
					monomersConsumed += count / self.numSimulations
					#calculated percentage of monomer remaining
					percentageRemaining = (startingMonomerAmount - monomersConsumed) / startingMonomerAmount
					monomercounts[index - 1] = percentageRemaining
			else:
				#adjust axis title
				frequencyPlot.set_ylabel("Average Total Monomer Occurences", labelpad=5, fontsize = 9)
				for index in polymerIndex:
					count = 0
					for polymer in polymerArray:
						if polymer[index - 1] == monomer:
							count += 1
					monomercounts[index - 1] = count / self.numSimulations
			#debugging purposes
			#print(polymerIndex)
			#print(monomercounts)
			#plots x and y arrays
			curve = frequencyPlot.plot(polymerIndex, monomercounts, label = "Monomer " + str(monomer))
			self.lineColors.append(curve[0].get_color())
		#legend-screw matplotlib; so fucking hard to format
		handles, labels = frequencyPlot.get_legend_handles_labels()
		lgd = frequencyPlot.legend(handles, labels, prop = {'size':7})
		# A tk.DrawingArea
		#imbedding matplotlib graph onto canvas
		self.canvas = FigureCanvasTkAgg(self.plotFigure, master = root)
		self.canvas.show()
		#Imbedding matplotlib toolbar onto canvas
		self.toolbar = NavigationToolbar2TkAgg(self.canvas, root)
		self.toolbar.update()
		self.canvas._tkcanvas.pack(side = Tk.BOTTOM, fill = Tk.BOTH, expand = 1)
		self.canvas.get_tk_widget().pack(side = Tk.BOTTOM, fill = Tk.BOTH, expand = 1)
	#Visualizes the polymers with colored squares representing monomers
	def visualizePolymers(self, polymerArray):
		numRows = self.numPolyToShow.get()
		if numRows > self.numPolymers * self.numSimulations: 
			numRows = self.numPolymers * self.numSimulations
		#parameters for canvas height and width
		canvasHeight = 120
		canvasWidth = 1000
		#Maximizes size of squares
		if (canvasHeight - 50) / numRows <= (canvasWidth - 50) / self.polymerLength:
			size = (canvasHeight - 50) / numRows
			canvasWidth = self.polymerLength * size + 20
		else:
			size = (canvasWidth - 50) / self.polymerLength
			canvasHeight = numRows * size + 10
		#neccesary instance variables
		lineColors = ["blue", "green", "red", "cyan", "magenta", "yellow", "black"]
		#LabelFrame for visualizeCanvas
		self.visualizationFrame = Tk.LabelFrame(master = root, text = "Polymer Visualization", padx = 3, pady = 3)
		self.visualizationFrame.pack(side = Tk.BOTTOM, fill = Tk.BOTH, expand = 0)
		#Canvas for visualization
		visualizeCanvas = Tk.Canvas(master = self.visualizationFrame, width = canvasWidth, height = canvasHeight)
		visualizeCanvas.pack()
		#Pad Parameters
		ulx = 20
		uly = 10
		#Visualizes polymers, number of polymers visualized based on numRows
		for row in range(0, numRows):
			#iterates through an array representation of monomer and adds a square with corresponding color
			for monomer in polymerArray[row]:
				color = lineColors[monomer - 1]
				visualizeCanvas.create_rectangle(ulx, uly + size * row, ulx + size, uly + size * (row + 1), fill = color)
				ulx += size
			ulx = 20
	#Converts an array of Tk.Entrys of ratios into an int array of starting monomer amounts
	#Note: might have result in total amount fo momoners being slightly more than inital total monomers due to ceiling divide
	def getMonomerAmounts(self):
		#A list of starting monomer amounts
		monomerAmounts = []
		#A list of starting ratios (not tkvars!)
		startingWeightList = []
		#Converts a list of Tkvars to a list of floats
		for entry in self.startingRatiosTkList:
			assert(float(entry.get()) > 0)
			startingWeightList.append(float(entry.get()))
		totalWeight = sum(startingWeightList)
		for weight in startingWeightList:
			#calculates number of monomers from monomer ratio and total monomers, ceiling dividing
			numMonomers = math.ceil(self.totalMonomers * weight / totalWeight)
			monomerAmounts.append(numMonomers)
		print("monomerAmounts: ", monomerAmounts)
		return monomerAmounts
	#Converts a 2D array of Tk.Entrys for coefficients into a 2D double array
	def getCoefficients(self):
		coeffList = []
		for entry in self.coefficientList:
			if entry == 0:
				coeffList.append(0)
			else:
				singleCoeffList = []
				coeffList.append(singleCoeffList)
				for innerEntry in entry: 
					if innerEntry == 0:
						singleCoeffList.append(0)
					elif ',' in list(innerEntry.get()):
						raise notInEuropeError("You're in America!")
					else:
						assert(float(innerEntry.get()) >= 0)
						singleCoeffList.append(float(innerEntry.get()))
		return coeffList
#When called, makes a pop out error informing user of invalid inputs
def errorMessage(message, width):
	#Toplevel parameters
	top = Tk.Toplevel()
	top.wm_title("Error")
	top.geometry("%dx%d%+d%+d" % (width, 70, 250, 125))
	#Message
	msg = Tk.Message(master = top, text = message, width = 500)
	msg.pack(side = Tk.TOP, pady = 5)
	#OK button to exit
	exitButton = Tk.Button(master = top, text = "Ok", command = top.destroy, width = 7)
	exitButton.pack(side = Tk.TOP, pady = 5)
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
def center(toplevel):
    toplevel.update_idletasks()
    w = toplevel.winfo_screenwidth()
    h = toplevel.winfo_screenheight()
    size = tuple(int(_) for _ in toplevel.geometry().split('+')[0].split('x'))
    x = w/2 - size[0]/2
    y = h/2 - size[1]/2
    toplevel.geometry("%dx%d+%d+%d" % (size + (x, y - 30)))
class notInEuropeError(Exception):
	def __init__(self, value):
		self.value = value
root = Tk.Tk()
root.wm_title("Compositional Drift v1.2")
app = Application(master = root)
app.mainloop()
