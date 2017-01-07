#All neccesary imports
import matplotlib
import random
import math
import json
import os.path
matplotlib.use('TkAgg')
from numpy import arange, sin, pi
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib import style
# implement the default mpl key bindings
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import sys
if sys.version_info[0] < 3:
    import Tkinter as Tk
    from Tkinter import ttk
    import Tkinter.filedialog
    import Tkinter.messagebox
    print(sys.version_info[0], "hi")
else:
    import tkinter as Tk
    from tkinter import ttk
    import tkinter.filedialog
    import tkinter.messagebox
    print(sys.version_info[0])
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
GRAPH_TYPE = 1
TOTAL_STARTING_MONOMERS = 1000
DP = 100
SETTING = 1
LOAD_SUCCESSFUL = False
GRAPH1_TYPE = "Monomer Occurences"
GRAPH2_TYPE = "Percentage Monomer"
HISTOGRAM1_MONOMER = 1
HISTOGRAM2_MONOMER = 2
HISTOGRAM_LIMIT = 0.8
PENULTIMATE = 0
VERSION = "v1.5.1"
CONFIGS = [["Number of Unique Monomers", 1], ["Number of Simulations", 1],
 ["Number of Polymers to Show", 1], 
 ["Graph Monomer Occurence", 1], ["Total Starting Monomers", 1], ["Degree of Polymerization", 1], 
 ["Default Setting", 1], ["Monomer Cap", 1], ["Graph 1 Type", 1], ["Graph 2 Type", 1], ["Histogram 1 Monomer", 1], ["Histogram 2 Monomer", 1],
 ["Percentage to Analyze for Histogram", 0], ["Penultimate", 1]]
#generates a config file if needed
def generateConfigFile():
	if not os.path.exists("config.txt"):
		file = open("config.txt", "w")
		file.write("Number of Unique Monomers = 2 \nNumber of Simulations = 200 \nNumber of Polymers to Show = 8 \n")
		file.write("Graph 1 Type = 0 \nGraph 2 Type = 1 \n")
		file.write("Histogram 1 Monomer = 1 \nHistogram 2 Monomer = 2 \nPercentage to Analyze for Histogram = 0.8 \n")
		file.write("Total Starting Monomers = 1000 \nDegree of Polymerization = 100 \nDefault Setting = 1 \nMonomer Cap = 5000000 \nPenultimate = 0 \n")
		file.write("Setting 1 \nNumber of Unique Monomers = 4 \nMonomer 1 Ratio = 50 \nMonomer 2 Ratio = 25 \nMonomer 3 Ratio = 20 \nMonomer 4 Ratio = 5 \n") 
		file.write("1-1 = 0.89 \n1-2 = 1 \n1-3 = 1 \n1-4 = 1 \n2-1 = 1 \n2-2 = 1.1 \n2-3 = 1.1 \n2-4 = 1.1 \n3-1 = 1 \n3-2 = 1.1 \n3-3 = 1.1 \n3-4 = 1.1 \n")
		file.write("4-1 = 1 \n4-2 = 1.1 \n4-3 = 1.1 \n4-4 = 1.1 \nend")
		file.close()
#reads config before initScreen for early variables
def earlyConfigRead():
	#error checking to see if config file exists
	if not os.path.exists("config.txt"):
		errorMessage("config.txt does not exist!", 220)
		return
	file = open("config.txt", "r")
	#variable to keep track of number of invalid config lines
	validLines = 0 
	notParsed = True
	for line in file:
		line = line.strip()
		intermediate = line.split("=")
		if len(intermediate) == 2:
			global NUM_UNIQUE_MONOMERS
			global SETTING
			global PENULTIMATE
			try:
				if intermediate[0].strip() == "Number of Unique Monomers" and notParsed:
					NUM_UNIQUE_MONOMERS = int(intermediate[1].strip())
					notParsed = False
					validLines += 1
				if intermediate[0].strip() == "Default Setting":
					SETTING = int(intermediate[1].strip())
					validLines += 1
				if intermediate[0].strip() == "Penultimate":
					PENULTIMATE= int(intermediate[1].strip())
					validLines += 1
					print("read penultimate!")
			except ValueError:
				continue
	print("earlyConfigRead valid lines: %s" %validLines)
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
	#global array to keep track of monomer ratios, initalized as an array of -1
	global RATIO_ARRAY
	#global array to keep track monomer coefficients, initialized of an array of arrays of -1
	global COEFF_ARRAY
	#global variable to shwo whether or not a setting loaded successfully
	global LOAD_SUCCESSFUL
	#global variable to kep track of last setting number
	global LAST_SETTING
	LAST_SETTING = 0
	LOAD_SUCCESSFUL = False
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
					RATIO_ARRAY = [-1] * numMonomers
					#array to keep track of all neccesary monomers
					numMonomerArray = [0] * numMonomers
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
					RATIO_ARRAY[int(lineArray[1]) - 1] = float(lineArray[4])
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
			else:
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
					#check if setting is bigger than the last setting
					if int(intermediate2[1]) > LAST_SETTING:
						LAST_SETTING = int(intermediate2[1])
					if int(intermediate2[1]) == SETTING:
						print("true setting", SETTING)
						parseSetting = True
				except AssertionError:
					invalidLines += 1
				continue
			invalidLines += 1
			continue
		try:
			#gets the config header
			configType = (line.split("="))[0].strip()
			#gets the config value
			configStringValue = (line.split("="))[1].strip()
		except:
			invalidLines += 1
			continue
		if (not setConfigVariable(configType, configStringValue)):
			invalidLines += 1
			continue
	file.close()
	print("Number of invalid config lines: ", invalidLines)
#helper method to set the static variable. Returns 1 if successful, 0 if not
def setConfigVariable(configType, configStringValue):
	for config in CONFIGS:
		if config[0] == configType:
			try:
				#handling for variable types
				if config[1] == 0:
					setConfigVariableHelper(configType, float(configStringValue))
				elif config[1] == 1:
					setConfigVariableHelper(configType, int(configStringValue))
				elif config[1] == 2:
					#print("configstring ", configStringValue)
					setConfigVariableHelper(configType, configStringValue)
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
	elif configType == "Total Starting Monomers":
		global TOTAL_STARTING_MONOMERS
		TOTAL_STARTING_MONOMERS = configValue
	elif configType == "Degree of Polymerization":
		global DP
		DP = configValue
	elif configType == "Default Setting":
		pass
	elif configType == "Monomer Cap":
		global MONOMER_CAP
		MONOMER_CAP = configValue
	elif configType == "Graph 1 Type":
		global GRAPH1_TYPE
		if configValue == 0:
			GRAPH1_TYPE = "Monomer Occurences"
		elif configValue == 1:
			GRAPH1_TYPE = "Percentage Monomer"
		elif configValue == 2:
			GRAPH1_TYPE = "Monomer Separation"
		elif configValue == 3:
			GRAPH1_TYPE = "None"
		else:
			assert False
	elif configType == "Graph 2 Type":
		global GRAPH2_TYPE
		if configValue == 0:
			GRAPH2_TYPE = "Monomer Occurences"
		elif configValue == 1:
			GRAPH2_TYPE = "Percentage Monomer"
		elif configValue == 2:
			GRAPH2_TYPE = "Monomer Separation"
		elif configValue == 3:
			GRAPH2_TYPE = "None"
		else:
			assert False
	elif configType == "Histogram 1 Monomer":
		global HISTOGRAM1_MONOMER
		HISTOGRAM1_MONOMER = configValue
	elif configType == "Histogram 2 Monomer":
		global HISTOGRAM2_MONOMER
		HISTOGRAM2_MONOMER = configValue
	elif configType == "Percentage to Analyze for Histogram":
		global HISTOGRAM_LIMIT
		HISTOGRAM_LIMIT = configValue
	elif configType == "Penultimate":
		global PENULTIMATE
		PENULTIMATE = configValue
	else:
		assert False, "shouldn't get here"	
#Main class 
class Application(ttk.Frame):
	def __init__(self, master = None):
		ttk.Frame.__init__(self, master)
		self.pack()
		self.initialize()
	#Initialization
	def initialize(self):
		try:
			#create config file
			generateConfigFile()
		except:
			errorMessage("Unable to generate config file!", 220)
		try:
			#early read on config file
			earlyConfigRead()
		except:
			errorMessage("Unable to read config file!", 220)
		#Creates the init screen
		self.initScreen()
		#Creates Input Widgets
		self.createInputWidgets()
		#Creates dummy visualization Frame
		self.visualizationFrame = ttk.Frame(master = root)
		self.destroyCanvas = False
		self.destroyCanvas2 = False
		self.destroyHide = False
	#Destroys unneccesary widgets
	def destroyWidgets(self):
		self.inputFrame.destroy()
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
			global SETTING
			try: 
				SETTING = int(self.settingTkVar.get())
				print("Setting", SETTING)
				assert int(SETTING) >= 0
			except:
				errorMessage("Please input valid parameters!", 220)
				return
			#reads config file
			try:
				readConfigFile()
			except:
				pass
			if LOAD_SUCCESSFUL:
				print("Load successful!")
				self.createMoreInputs()
			else:
				errorMessage("Unable to load config settings! Please fix config file.", 300)
		#Confirms number of monomers, creates more input widgets
		def enter(self):
			#read config file
			try:
				readConfigFile()
			except:
				pass
			if LOAD_SUCCESSFUL:
				print("Load Successful!")
			else:
				errorMessage("Unable to load settings! Please fix config file.", 300)
			#nonlocal numMonomers 
			self.numMonomers = int(self.monomerCountTkVar.get())
			try:
				#asserting that input is in correct range
				assert self.numMonomers < 8
				assert self.monomerCountTkVar.get() > 0
			except:
				errorMessage("Please input valid paramters!", 220)
			global useLoadedSettings
			useLoadedSettings = False
			#print(numMonomers) # debugging purposes
			self.createMoreInputs() 
		#Clears unneccesary widgets and creates more input widgets; called by enter()
		def _quit():
			root.quit()
			root.destroy()
		#The parent LabelFrame for all Input Widgets
		self.inputFrame = ttk.Frame(master = root, borderwidth = 1, relief = "ridge")
		self.inputFrame.pack(side = Tk.BOTTOM, fill = Tk.X, expand = 0, padx = 7, pady = 7)
		#top seperator
		#self.topSep = ttk.Separator(master = self.inputFrame, orient = Tk.HORIZONTAL)
		#self.topSep.pack(side = Tk.TOP, fill = Tk.BOTH)
		#leftmost separator
		#self.sideSep1 = ttk.Separator(master = self.inputFrame, orient = Tk.VERTICAL)
		#self.sideSep1.pack(side = Tk.LEFT, fill = Tk.BOTH)
		#A frame for column 1 of inputs
		self.columnFrame = ttk.Frame(master = self.inputFrame)
		self.columnFrame.pack(side = Tk.LEFT, padx = 5, pady = 5)
		#A frame for row 1 of inputs
		self.rowFrame1 = ttk.Frame(master = self.columnFrame)
		self.rowFrame1.pack(side = Tk.TOP, padx = 5, pady = 0)
		#A frame for row 2 of inputs
		self.rowFrame2 = ttk.Frame(master = self.columnFrame)
		self.rowFrame2.pack(side = Tk.TOP, padx = 5, pady = 5)
		#Label for monomerCount spinbox
		self.monomerCountLabel = ttk.Label(master = self.rowFrame1, text = "Number of Unique Monomers:")
		self.monomerCountLabel.pack(side = Tk.LEFT, padx = 0, pady = 0)
		#tkvar for monomercount
		self.monomerCountTkVar = Tk.IntVar()
		#MonomerCount spinbox
		self.monomerCountSpinbox = Tk.Spinbox(master = self.rowFrame1, from_ = 2, to = 7, width = 2, textvariable = self.monomerCountTkVar)
		self.monomerCountSpinbox.pack(side = Tk.LEFT, padx = 5, pady = 0)
		#sets monomerCOuntTkVar to default setting
		self.monomerCountTkVar.set(NUM_UNIQUE_MONOMERS)
		#countConfirm Button
		self.countConfirmButton = ttk.Button(master = self.rowFrame1, text = "Enter",
		command = lambda:enter(self), width = 9)
		self.countConfirmButton.pack(side = Tk.LEFT, padx = 5, pady = 5)
		#Label for msetting spinbox
		self.monomerCountLabel = ttk.Label(master = self.rowFrame2, text = "Setting Number to Use:")
		self.monomerCountLabel.pack(side = Tk.LEFT, padx = 0, pady = 0)
		#tkvar for setting
		self.settingTkVar = Tk.IntVar()
		#settingCount spinbox
		self.settingSpinbox = Tk.Spinbox(master = self.rowFrame2, from_ = 1, to = 1000, width = 2, textvariable = self.settingTkVar)
		self.settingSpinbox.pack(side = Tk.LEFT, padx = 5, pady = 0)
		#sets settingTkVar to 1
		global SETTING
		self.settingTkVar.set(SETTING)
		#use loaded inputs button
		self.loadButton = ttk.Button(master = self.rowFrame2, text = "Load from Settings",
		command = lambda:loadSettings(self), width = 18)
		self.loadButton.pack(side = Tk.LEFT, padx = 5, pady = 5)
		#frame for penultimate model choice, 3 row of inputs
		self.penultimateFrame = ttk.Frame(master = self.columnFrame)
		self.penultimateFrame.pack(side = Tk.TOP, padx = 5, pady = 0)
		#Label for penultimate checkbox
		self.penultimateLabel = ttk.Label(master = self.penultimateFrame, text = "Use Penultimate Model:")
		self.penultimateLabel.pack(side = Tk.LEFT, padx = 0, pady = 0)
		self.penultimateCheckButton = ttk.Checkbutton(master = self.penultimateFrame, text = None)
		self.penultimateCheckButton.pack(side = Tk.LEFT, padx = 5)
		#Setting penultimateTkVar to PENULTIMATE global variable
		self.penultimateTkVar = Tk.IntVar()
		self.penultimateCheckButton["variable"] = self.penultimateTkVar
		self.penultimateTkVar.set(PENULTIMATE)
		#rightmost separator
		#self.sideSep2 = ttk.Separator(master = self.inputFrame, orient = Tk.VERTICAL)
		#self.sideSep2.pack(side = Tk.LEFT, fill = Tk.BOTH)
		#bottom separator
		#self.bottomSep = ttk.Separator(master = self.inputFrame, orient = Tk.HORIZONTAL)
		#self.bottomSep.pack(side = Tk.BOTTOM, fill = Tk.BOTH)
	#Creates more input widgets based on numMonomers
	def createMoreInputs(self):
		root.maxsize(width = 10000, height = 10000)
		root.resizable(width = True, height = True)
		root.sizefrom(who = "program")
		#case for loading inputs
		if LOAD_SUCCESSFUL and useLoadedSettings:
			self.numMonomers = numMonomers
		else: 
			try:
				assert(self.monomerCountTkVar.get() > 0)
				global PENULTIMATE
				PENULTIMATE = self.penultimateTkVar.get()
			except:
				errorMessage("Please input valid parameters!", 220)
				return
		#Destroys or edits current widgets
		self.initFrame.destroy()
		self.loadButton.destroy()
		self.rowFrame1.destroy()
		self.rowFrame2.destroy()
		self.penultimateFrame.destroy()
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

		#Frame for totalMonomers label and Entry
		self.totalMonomersFrame = ttk.Frame(master = self.columnFrame)
		self.totalMonomersFrame.pack(side = Tk.TOP)
		#Label for totalMonomers Entry
		self.totalMonomersLabel = ttk.Label(master = self.totalMonomersFrame, text = "Total Starting Monomers:")
		self.totalMonomersLabel.pack(side = Tk.LEFT, pady = 3)
		#Entry for totalMonomers
		self.totalMonomersEntry = ttk.Entry(master = self.totalMonomersFrame, width = 5)
		self.totalMonomersEntry.pack(side = Tk.LEFT, padx = 3, pady = 3)
		#Setting totalMonomers to 1000
		self.totalMonomersTkVar = Tk.IntVar()
		self.totalMonomersEntry["textvariable"] = self.totalMonomersTkVar
		self.totalMonomersTkVar.set(TOTAL_STARTING_MONOMERS)
		#Frame for raftRatio label and Entry
		self.raftRatioFrame = ttk.Frame(master = self.columnFrame)
		self.raftRatioFrame.pack(side = Tk.TOP)
		#Label for raftRatio Entry
		self.raftRatioLabel = ttk.Label(master = self.raftRatioFrame, text = "Degree of Polymerization:")
		self.raftRatioLabel.pack(side = Tk.LEFT, pady = 3)
		#Entry for raftRatio
		self.raftRatioEntry = ttk.Entry(master = self.raftRatioFrame, width = 5)
		self.raftRatioEntry.pack(side = Tk.LEFT, padx = 3, pady = 3)
		self.raftRatioTkVar = Tk.StringVar()
		self.raftRatioEntry["textvariable"] = self.raftRatioTkVar
		self.raftRatioTkVar.set(DP)
		#Frame for numSimulations label and Entry
		self.numSimsFrame = ttk.Frame(master = self.columnFrame)
		self.numSimsFrame.pack(side = Tk.TOP)
		#Label for numSims Entry
		self.numSimsLabel = ttk.Label(master = self.numSimsFrame, text = "Number of Simulations:   ")
		self.numSimsLabel.pack(side = Tk.LEFT, pady = 3)
		#Entry for numSimulations
		self.numSimsEntry = ttk.Entry(master = self.numSimsFrame, width = 5)
		self.numSimsEntry.pack(side = Tk.LEFT, padx = 3, pady = 3)
		#Setting number of simulations to 1000
		self.numSimsTkVar = Tk.IntVar()
		self.numSimsEntry["textvariable"] = self.numSimsTkVar
		self.numSimsTkVar.set(NUM_SIMULATIONS)
		#Frame for numPolyToShow SpinBox
		self.numPolyToShowFrame = ttk.Frame(master = self.columnFrame)
		self.numPolyToShowFrame.pack(side = Tk.TOP)
		#Label for numPolyToShow spinbox
		self.numPolyToShowLabel = ttk.Label(master = self.numPolyToShowFrame, text = "Number of Polymers to Show:")
		self.numPolyToShowLabel.pack(side = Tk.LEFT, padx = 0, pady = 0)
		#numPolyToShow spinbox
		self.numPolyToShowBox = Tk.Spinbox(master = self.numPolyToShowFrame, from_ = 1, to = 12, width = 2)
		self.numPolyToShowBox.pack(side = Tk.LEFT, padx = 5, pady = 0)
		#Setting numPolyToShow to 6
		self.numPolyToShow = Tk.IntVar()
		self.numPolyToShowBox["textvariable"] = self.numPolyToShow
		self.numPolyToShow.set(NUM_POLYMERS_SHOW)
		#Frame for Back, Quit, and Simulate buttons
		self.backSimFrame = ttk.Frame(master = self.columnFrame)
		self.backSimFrame.pack(side = Tk.TOP)
		#A simulate button to simulate polymer formation
		self.simulateButton = ttk.Button(master = self.backSimFrame, text = "Simulate", width = 9,
		 command = self.simulate)
		self.simulateButton.pack(side = Tk.LEFT, padx = 6, pady = 4)
		#A back button to enter a diff number of monomers
		self.backButton = ttk.Button(master = self.backSimFrame, text = "Back", width = 7,
		 command = lambda:back(self))
		self.backButton.pack(side = Tk.LEFT, padx = 6, pady = 4)	
		#An Update Button Widget
		updateButton = ttk.Button(master = self.backSimFrame, text = "Update",
		 command = lambda:self.plotCompositions(True), width = 7)
		updateButton.pack(side = Tk.LEFT, padx = 6)
		#A save state button
		self.saveButton = ttk.Button(master = self.backSimFrame, text = "Save", width = 7,
		 command = self.saveState)
		self.saveButton.pack(side = Tk.LEFT, padx = 6, pady = 4)
		#seperator for column
		self.col1Sep = ttk.Separator(master = self.inputFrame, orient = Tk.VERTICAL)
		self.col1Sep.pack(side = Tk.LEFT, expand = True, fill = Tk.BOTH, pady = 1)
		#Frame for graph options
		self.column2Frame = ttk.Frame(master = self.inputFrame)
		self.column2Frame.pack(side = Tk.LEFT, padx = 8)
		#tkvar for graphType1
		self.graphType1TkVar = Tk.StringVar()
		#tkVar for graphType2
		self.graphType2TkVar = Tk.StringVar()
		#frame for graphType1
		self.graphFrame1 = ttk.Frame(master = self.column2Frame)
		self.graphFrame1.pack(side = Tk.TOP, pady = 3)
		#frame for graphType2
		self.graphFrame2 = ttk.Frame(master = self.column2Frame)
		self.graphFrame2.pack(side = Tk.TOP, pady = 3)
		#Label for graphType1 spinbox
		self.graphType1Label = ttk.Label(master = self.graphFrame1, text = "Graph 1 Type:")
		self.graphType1Label.pack(side = Tk.LEFT)
		#Label for graphType2 spinbox
		self.graphType2Label = ttk.Label(master = self.graphFrame2, text = "Graph 2 Type:")
		self.graphType2Label.pack(side = Tk.LEFT)
		#combobox
		self.graphType1ComboBox = ttk.Combobox(master = self.graphFrame1, values = ("Monomer Occurences", "Percentage Monomer", 
			"Monomer Separation", "None"), textvariable = self.graphType1TkVar, state = "readonly")
		self.graphType1ComboBox.pack(side = Tk.LEFT)

		#Spinbox for graphType1
		#self.graphType1Spinbox = Tk.Spinbox(master = self.graphFrame1, values = ("Monomer Occurences", "Percentage Monomer", 
		#	"Monomer Separation", "None"), width = 20, textvariable = self.graphType1TkVar, state = "readonly")
		#self.graphType1Spinbox.pack(side = Tk.LEFT)
		#setting default variable for graphType1
		self.graphType1TkVar.set(GRAPH1_TYPE)
		#combobox
		self.graphType2ComboBox = ttk.Combobox(master = self.graphFrame2, values = ("Monomer Occurences", "Percentage Monomer", 
			"Monomer Separation", "None"), textvariable = self.graphType2TkVar, state = "readonly")
		self.graphType2ComboBox.pack(side = Tk.LEFT)
		#Spinbox for graphType2
		#self.graphType2Spinbox = Tk.Spinbox(master = self.graphFrame2, values = ("Monomer Occurences", "Percentage Monomer", 
		#"Monomer Separation", "None"), width = 20, textvariable = self.graphType2TkVar, state = "readonly")
		#self.graphType2Spinbox.pack(side = Tk.LEFT)
		#setting default variable for graphType2
		self.graphType2TkVar.set(GRAPH2_TYPE)
		#Frame for histogramMonomer1
		self.histogramMonomer1Frame = ttk.Frame(master = self.column2Frame)
		self.histogramMonomer1Frame.pack(side = Tk.TOP)
		#Label for histogramMonomer1
		self.histogramMonomer1Label = ttk.Label(master = self.histogramMonomer1Frame, text = "Histogram 1 Monomer:")
		self.histogramMonomer1Label.pack(side = Tk.LEFT)
		#tkvar for histogramMonomer1
		self.histogramMonomer1TkVar = Tk.IntVar()
		#spinbox for histogramMonomer1
		self.histogramMonomer1Spinbox  = Tk.Spinbox(master = self.histogramMonomer1Frame, from_ = 1, to = self.numMonomers, width = 2, state = "readonly")
		self.histogramMonomer1Spinbox.pack(side = Tk.LEFT)
		#setting spinbox to default value
		self.histogramMonomer1Spinbox["textvariable"] = self.histogramMonomer1TkVar
		self.histogramMonomer1TkVar.set(HISTOGRAM1_MONOMER)
		#Frame for histogramMonomer2
		self.histogramMonomer2Frame = ttk.Frame(master = self.column2Frame)
		self.histogramMonomer2Frame.pack(side = Tk.TOP)
		#Label for histogramMonomer2
		self.histogramMonomer2Label = ttk.Label(master = self.histogramMonomer2Frame, text = "Histogram 2 Monomer:")
		self.histogramMonomer2Label.pack(side = Tk.LEFT)
		#tkvar for histogramMonomer2
		self.histogramMonomer2TkVar = Tk.IntVar()
		#spinbox for histogramMonomer2
		self.histogramMonomer2Spinbox  = Tk.Spinbox(master = self.histogramMonomer2Frame, from_ = 1, to = self.numMonomers, width = 2, state = "readonly")
		self.histogramMonomer2Spinbox.pack(side = Tk.LEFT)
		#setting spinbox to default value
		self.histogramMonomer2Spinbox["textvariable"] = self.histogramMonomer2TkVar
		self.histogramMonomer2TkVar.set(HISTOGRAM2_MONOMER)
		#Frame for histogramLimit
		self.histogramLimitFrame = ttk.Frame(master = self.column2Frame)
		self.histogramLimitFrame.pack(side = Tk.TOP)
		#Label for histogramLimit
		self.histogramLimitLabel = ttk.Label(master = self.histogramLimitFrame, text = "Percentage to Analyze for Histogram:")
		self.histogramLimitLabel.pack(side = Tk.LEFT)
		#entry for histogramLimit
		self.histogramLimitEntry = ttk.Entry(master = self.histogramLimitFrame, width = 4)
		self.histogramLimitEntry.pack(side = Tk.LEFT)
		#tkvar for histogramLimit
		self.histogramLimitTkVar = Tk.StringVar()
		self.histogramLimitEntry["textvariable"] = self.histogramLimitTkVar
		self.histogramLimitTkVar.set(HISTOGRAM_LIMIT)
		#seperator for column2
		self.col2Sep = ttk.Separator(master = self.inputFrame, orient = Tk.VERTICAL)
		self.col2Sep.pack(side = Tk.LEFT, expand = True, fill = Tk.BOTH, padx = 1, pady = 1)
		#Frame for Monomer Amounts
		self.amountFrame = ttk.Frame(master = self.inputFrame) 
		self.amountFrame.pack(side = Tk.LEFT, padx = 0)
		#Frame for Monomer Coefficients
		self.coefficientFrame = ttk.Frame(master = self.inputFrame)
		self.coefficientFrame.pack(side = Tk.LEFT, padx = 0)
		# A list of ttk.Entry objects for Monomer ratios
		self.startingRatiosTkList = [] 
		# A 2D list of ttk.Entry objects for Coefficicients
		self.coefficientList = [] 
		#a list of the ttk.Entry objects for monomer input
		self.monomerAmountTkVarArray = []
		#variable for keeping track of monomer ratio loop
		createCount = 0;
		#While loop creating number of neccesary amount Entry boxes
		while createCount < self.numMonomers:
			#Label for inputAmount
			monomerAmountFrame = ttk.Frame(master = self.amountFrame)
			monomerAmountFrame.pack(side = Tk.TOP, padx = 5, pady = 3)
			inputAmountLabel = ttk.Label(master = monomerAmountFrame, text = "     Monomer " 
				+ str(createCount + 1) + " Ratio:")
			inputAmountLabel.pack(side = Tk.LEFT)
			#Entry for inputAmount
			amount = Tk.IntVar()
			inputAmount = ttk.Entry(master = monomerAmountFrame, width = 5, textvariable = amount, text = 20)
			inputAmount.pack(side = Tk.LEFT, padx = 2)
			#Setting Default Value to 20
			inputAmount["textvariable"] = amount
			if LOAD_SUCCESSFUL and useLoadedSettings:
				amount.set(RATIO_ARRAY[createCount])
			else:
				amount.set(20)
			self.monomerAmountTkVarArray.append(amount)
			#Add ttk.Entry object to startingAmountList
			self.startingRatiosTkList.append(inputAmount)
			createCount += 1
		setDefaultCount = 0
		"""while setDefaultCount < self.numMonomers:
			self.startingRatiosTkList[setDefaultCount].configure(text = 20)
			setDefaultCount += 1"""
		#Debugging purposes
		#print("startingAmountList: ", self.startingRatiosTkList) 
		if not PENULTIMATE:
			createCount2 = 0
			#IMPORTANT: list so that inputs can be udated and not trash collected
			self.coeffTkVarArray = []
			#While loop creating number of neccesary coefficient Entry boxes
			while createCount2 < self.numMonomers:
				combinations = 0
				#Appends to coefficient list a list containing coefficients for the polymer index
				singleMonoCoeffList = []
				self.coefficientList.append(singleMonoCoeffList)
				#Frame for Coefficients for Single Monomer
				singleCoeffFrame = ttk.Frame(master = self.coefficientFrame)
				singleCoeffFrame.pack(side = Tk.LEFT, fill = Tk.X, expand = 1)
				while combinations < self.numMonomers:				
					#Label for inputAmount
					coeffValFrame = ttk.Frame(master = singleCoeffFrame)
					coeffValFrame.pack(side = Tk.TOP, padx = 5, pady = 3)
					inputCoeffLabel = ttk.Label(master = coeffValFrame, text = str(createCount2 + 1)
					 + "-" + str(combinations + 1) + " Constant:" )
					inputCoeffLabel.pack(side = Tk.LEFT)
					#Entry for inputAmount
					inputCoeff = ttk.Entry(master = coeffValFrame, width = 4)
					inputCoeff.pack(side = Tk.LEFT, padx = 5)
					#Setting Default Coefficient to 1
					coeff = Tk.IntVar()
					inputCoeff["textvariable"] = coeff
					if LOAD_SUCCESSFUL and useLoadedSettings:
						coeff.set(COEFF_ARRAY[createCount2][combinations])
					else:
						coeff.set(1)
					#IMPORTANT: saves input so won't be garbage collected. DO NOT DELETE THIS LINE< EVEN IF IT SEEMS USELESS
					self.coeffTkVarArray.append(coeff)
					#Add a ttk.Entry object to singleMonoCoeffList
					singleMonoCoeffList.append(inputCoeff)
					combinations += 1
				createCount2 += 1
		#if PENULTIMATE is true, create input entries for penultimate model
		if PENULTIMATE:
			self.coeffTkVarArray = []
			#for loop creating neccesary number of penultimate coefficient entry boxes
			for nextMonomer in range(0, self.numMonomers):
				nextMonomerList = []
				self.coefficientList.append(nextMonomerList)
				#Frame for Coefficients for Single Monomer
				singleCoeffFrame = ttk.Frame(master = self.coefficientFrame)
				singleCoeffFrame.pack(side = Tk.LEFT, fill = Tk.X, expand = 1)
				for ultimate in range(0, self.numMonomers):
					CoeffList = []
					nextMonomerList.append(CoeffList)
					for penultimate in range(0, self.numMonomers):
						#Label for inputAmount
						coeffValFrame = ttk.Frame(master = singleCoeffFrame)
						coeffValFrame.pack(side = Tk.TOP, padx = 5, pady = 3)
						inputCoeffLabel = ttk.Label(master = coeffValFrame, text = str(penultimate + 1) + "-" + str(ultimate + 1)
						 + "-" + str(nextMonomer + 1) + " Constant:" )
						inputCoeffLabel.pack(side = Tk.LEFT)
						#Entry for inputAmount
						inputCoeff = ttk.Entry(master = coeffValFrame, width = 4)
						inputCoeff.pack(side = Tk.LEFT, padx = 5)
						#Setting Default Coefficient to 1
						coeff = Tk.IntVar()
						inputCoeff["textvariable"] = coeff
						if LOAD_SUCCESSFUL and useLoadedSettings:
							coeff.set(COEFF_ARRAY[createCount2][combinations])
						else:
							coeff.set(1)
						self.coeffTkVarArray.append(coeff)
						#Add a ttk.Entry object to singleMonoCoeffList
						CoeffList.append(inputCoeff)

		# Syntax": for i in range(0,x)
		# fpr i in iterable
		#Debugging Purposes: prints the 2D array which should store coefficients 	
		#print("coefficientList: " , self.coefficientList)
		#counter = 0
		"""for coefflist in self.coefficientList:
			print("list at index " + str(counter) + ": ", coefflist)
			counter += 1"""
	#Frame and contents for opening screen
	def initScreen(self):
		root.maxsize(width = 330, height = 270)
		self.initFrame = ttk.Frame(master = root)
		self.initFrame.pack(side = Tk.TOP, fill = Tk.X, expand = 0, padx = 3, pady = 5)
		message1 = Tk.Message(master = self.initFrame, width = 250, font = ("Times New Roman", 10, "bold"),
		 text = "Compositional Drift Simulator %s" % VERSION)
		message1.pack(side = Tk.TOP, padx = 0, pady = 0)
		message2 = Tk.Message(master = self.initFrame, width = 200,
			text = "Author: Vincent Wu")
		message2.pack(side = Tk.TOP)
		message3 = Tk.Message(master = self.initFrame, width = 300,
			text = "Welcome to Compositional Drift Simulator! This program uses the Mayo Lewis Equation and the" 
			" Monte Carlo method to simulate the growth of a copolymer chain."
			" Please input your desired number of unique monomers or choose a default setting.")
		message3.pack(side = Tk.TOP)
		root.resizable(width=True, height=True)
		root.sizefrom(who = "program")
	#Hides input params
	def showInputParams(self):
		self.inputFrame.pack(side = Tk.BOTTOM, pady = 3)
		self.showButton.destroy()
	#Shows input params
	def hideInputParams(self):
		self.inputFrame.pack_forget()
		self.showButton = ttk.Button(master = self.parentFrame, text = "Show Input Parameters", width = 27,
			bg = "pale turquoise", activebackground = "light slate blue", command = self.showInputParams)
		self.showButton.pack(side = Tk.TOP)
	#Simulates polymer reaction based on input values
	def simulate(self):
		#print("Simulating!")
		#Asserts that there are no input errors. Shows errorMessage if error caught.test
		try:
			self.totalMonomers = int(self.totalMonomersTkVar.get())
			monomerAmounts = self.getMonomerAmounts()
			if not PENULTIMATE:
				singleCoeffList = self.getCoefficients()
			else:
				singleCoeffList = self.getPenultimateCoeff()
			self.numSimulations = int(self.numSimsTkVar.get())
			self.raftRatio = float(self.raftRatioTkVar.get())
			self.histogramLimit = float(self.histogramLimitTkVar.get())
			assert(self.histogramLimit <= 1)
			assert(self.histogramLimit > 0)
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
		#progress bar widget setup
		simulationProgressTkVar = Tk.IntVar()
		self.simulationProgressBar = ttk.Progressbar(master = self.column2Frame, mode = "determinate", 
			variable = simulationProgressTkVar, maximum = self.numSimulations, length = 200)
		self.simulationProgressBar.pack(side = Tk.TOP, padx = 5, pady = 2)
		self.simulationProgressBar.update()
		simulationProgressTkVar.set(0)
		currProgress = 0
		self.simulationProgressBar.start()
		"""Calculates Initial Conditions based on inputs"""
		#number of monomers in each polymer assuming reaction goes to completion, based on raftRatio and itotalMonomers
		#print(self.totalMonomers)
		#print(self.raftRatio)
		self.polymerLength = int(self.raftRatio)
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
		#print("Polymer Length: ", self.polymerLength)
		#print("Number of Polymers: ", self.numPolymers)
		#destroys hideButton if necessary
		if self.destroyHide:
			self.hideButton.destroy()
		#hideButton creation, _DEPRECATED_
		"""self.hideButton = ttk.Button(master = self.buttonFrame, bg = "pale turquoise", activebackground = "light slate blue",
			width = 27, text = "Hide Input Paramters", command = self.hideInputParams)
		self.hideButton.pack(side = Tk.TOP, pady = 3)
		self.destroyHide = True"""
		#print("monomerAmounts: ", monomerAmounts)
		#print("singleCoeffList: ", singleCoeffList)
		#An array of polymers
		self.polymerArray = []
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
			#case for penultimate model: choosing the second monomer in chain, based off of HETEROGENOUS constants (works only for 2 monomer system)
			if PENULTIMATE:
				#iterating through each polymer in batch
				for polymer in localPolymerArray:
					choices = []
					#calculates weight chance for each monomer
					for monomerID in range(1, self.numMonomers + 1):
						#retrieving coefficient based on previous and currrent monomer, assumes hetergoenous penultimate monomer
						coeff = singleCoeffList[monomerID - 1][polymer[-1] - 1][1]
						#print("coeff: ", coeff)
						# weight chance calulations for monomer attaching: coefficient*(amount of monomer remaining)
						chance = monomerAmounts[monomerID - 1] * coeff
						#adds a two element list to choices containing monomer and weight
						choices.append([monomerID, chance])
					#Using weighted_choice, selects next monomer
					try:
						#print(choices)
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
			#case for penultimate model: propogating polymer chain once first 2 monomers have been decided
			if PENULTIMATE:
				for currPolymer in range(2, self.polymerLength):	
					for polymer in localPolymerArray:
						choices = []
						#calculates weight chance for each monomer
						for monomerID in range(1, self.numMonomers + 1):
							#retrieving coefficient based on previous and currrent monomer, assumes hetergoenous penultimate monomer
							#assigning variables for clarity
							nextMonomer = monomerID - 1
							ultimateMonomer = polymer[-1] - 1
							penultimateMonomer = polymer[-2] - 1
							coeff = singleCoeffList[nextMonomer][ultimateMonomer][penultimateMonomer]
							# weight chance calulations for monomer attaching: coefficient*(amount of monomer remaining)
							chance = monomerAmounts[monomerID - 1] * coeff
							#adds a two element list to choices containing monomer and weight
							choices.append([monomerID, chance])
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
			#case for non-penultimate simulation
			else:
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
			self.polymerArray = self.polymerArray + localPolymerArray
			#increases simCounter by 1
			simCounter += 1
			#increase progress
			currProgress += 1
			simulationProgressTkVar.set(currProgress)
			self.simulationProgressBar.update()
			#print(currProgress)
		#Debugging purposes
		"""Important Debug: Prints out array of polymers"""
		#print("Array of Polymers: ", polymerArray)
		#stops progress bar
		self.simulationProgressBar.stop()
		#outputs polymer onto textfile
		text_file = open("polymerArray.txt", "w")
		json.dump(self.polymerArray, text_file)
		text_file.close()
		#destroys canvas if necessary
		if self.destroyCanvas:
			self.canvas.get_tk_widget().destroy()
			self.toolbar.destroy()
		self.destroyCanvas = True
		self.visualizationFrame.destroy()
		self.visualizePolymers(self.polymerArray)
		self.plotCompositions(False)
		#destroy progress bar
		self.simulationProgressBar.destroy()
		center(root)
		#self.inputFrame.pack_forget()
	#plots compositions given a PolymerArray
	def plotCompositions(self, update):
		try:
			polymerArray = self.polymerArray
		except:
			errorMessage("Please simulate first!", 300)
			return
		#destroys canvas if necessary
		if update:
			self.canvas.get_tk_widget().destroy()
			self.toolbar.destroy()
		#style to use
		style.use("bmh")
		#retrieving graphType variables
		self.graph1Type = self.graphType1TkVar.get()
		self.graph2Type = self.graphType2TkVar.get()
		#Plot and Figure formatting
		self.plotFigure = Figure(figsize=(5.5, 3.3), dpi=100)
		if self.graph1Type == "None" and self.graph2Type == "None":
			return
		if self.graph1Type == "None":
			self.subplot1 = self.plotFigure.add_subplot(111)
			self.subplot2 = self.plotFigure.add_subplot(122)
			self.subplot1.set_color_cycle(self.colorArray)
			self.subplot1.tick_params(labelsize = 7)
			self.graphSubPlot(polymerArray, self.graph2Type, self.subplot1, 1)
		elif self.graph2Type == "None":
			self.subplot1 = self.plotFigure.add_subplot(111)
			self.subplot1.set_color_cycle(self.colorArray)
			self.subplot1.tick_params(labelsize = 7)
			self.graphSubPlot(polymerArray, self.graph1Type, self.subplot1, 1)
		else:
			self.subplot1 = self.plotFigure.add_subplot(121)
			self.subplot2 = self.plotFigure.add_subplot(122)
			self.subplot1.set_color_cycle(self.colorArray)
			self.subplot2.set_color_cycle(self.colorArray)
			self.subplot1.tick_params(labelsize = 7)
			self.subplot2.tick_params(labelsize = 7)
			self.graphSubPlot(polymerArray, self.graph1Type, self.subplot1, 1)
			self.graphSubPlot(polymerArray, self.graph2Type, self.subplot2, 2)	
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
		#LabelFrame for visualizeCanvas
		self.visualizationFrame = ttk.LabelFrame(master = root, text = "Polymer Visualization")
		self.visualizationFrame.pack(side = Tk.BOTTOM, fill = Tk.BOTH, expand = 0, padx = 7, pady = 0)
		#update visuals to get correct sizing
		self.visualizationFrame.update()
		#variable to keep track of frame width
		self.visualFrameWidth = self.visualizationFrame.winfo_width()
		#print("Width", self.visualFrameWidth)
		numRows = self.numPolyToShow.get()
		if numRows > self.numPolymers * self.numSimulations: 
			numRows = self.numPolymers * self.numSimulations
		#parameters for canvas height and width
		canvasHeight = 130
		canvasWidth = self.visualFrameWidth
		#Maximizes size of squares
		if (canvasHeight - 25) / numRows <= (canvasWidth - 40) / self.polymerLength:
			size = (canvasHeight - 25) / numRows
			canvasWidth = self.polymerLength * size + 20
		else:
			size = (canvasWidth - 40) / self.polymerLength
			canvasHeight = numRows * size + 10
		self.visualFrameWidth = self.visualizationFrame.winfo_height()
		#Canvas for visualization
		visualizeCanvas = Tk.Canvas(master = self.visualizationFrame, width = canvasWidth, height = canvasHeight)
		visualizeCanvas.pack()
		#colors
		#line colors to use
		self.colorArray = ['#4D4D4D','#5DA5DA', '#F15854', '#DECF3F', '#60BD68', '#F17CB0', '#B276B2', '#FAA43A']
		#Pad Parameters
		ulx = 20
		uly = 10
		#Visualizes polymers, number of polymers visualized based on numRows
		for row in range(0, numRows):
			#iterates through an array representation of monomer and adds a square with corresponding color
			for monomer in polymerArray[row]:
				color = self.colorArray[monomer - 1]
				visualizeCanvas.create_rectangle(ulx, uly + size * row, ulx + size, uly + size * (row + 1), fill = color)
				ulx += size
			ulx = 20
	#Converts an array of ttk.Entrys of ratios into an int array of starting monomer amounts
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
		#print("monomerAmounts: ", monomerAmounts)
		return monomerAmounts
	#Converts a 2D array of ttk.Entrys for coefficients into a 2D double array
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
	def getPenultimateCoeff(self):
		coeffList = []
		for nextMonomerList in self.coefficientList:
			ultimateList = []
			coeffList.append(ultimateList)
			for ultimateMonomerList in nextMonomerList:
				penultimateCoeffList = []
				ultimateList.append(penultimateCoeffList)
				for penultimateCoeffEntry in ultimateMonomerList:
					if ',' in list(penultimateCoeffEntry.get()):
						raise notInEuropeError("You're in America!")
					else:
						assert(float(penultimateCoeffEntry.get()) >= 0)
						penultimateCoeffList.append(float(penultimateCoeffEntry.get()))
		print("coeffList: ", coeffList)
		return coeffList


	#returns an array of numbers for each consecutive monomer; to be used in histogram plotting
	def getHistogramData(self, polymerArray, monomerID, indexLimit):
		#initializing array to be returned
		histogramData = []
		#count through all polymers
		for polymer in polymerArray: 
			#constant to keep track of polymer index, to know when to stop counting
			polymerIndex = 1
			numConsecutive = 0
			#iterate through each polymer until monomerLimit is hit and count consecutive polymers
			for monomer in polymer:
				#if index limit is reached, add any consecutives to histogram, and break from for loop
				if polymerIndex > indexLimit:
					if numConsecutive > 0:
						count = 0
						while count < numConsecutive:
							histogramData.append(numConsecutive)
							count += 1
					break
				#if monomer is not consecutive and monomer before was, add consecutive number to data and reset values
				if monomer != monomerID and numConsecutive > 0:
					count = 0
					while count < numConsecutive:
						histogramData.append(numConsecutive)
						count += 1
					numConsecutive = 0
					polymerIndex += 1
					continue
				#increment consecutive counter by 1 if monomer is consecutive
				if  monomer == monomerID:
					numConsecutive += 1
					polymerIndex += 1
					continue
		return histogramData
	def graphSubPlot(self, polymerArray, graphType, subplot, number):
		if graphType == "Percentage Monomer" or graphType == "Monomer Occurences":
			#Iterates through each unique monomer and plots composition
			for monomer in range(1, self.numMonomers + 1):
				#x-axis array
				polymerIndex = list(range(1, self.polymerLength + 1))
				#y-axis array initation
				monomercounts = [0] * self.polymerLength
				#inputs counts into y-axis array
				#graphs Percentage of Monomer Remaining
				if graphType == "Percentage Monomer":
					#adjust axis title
					subplot.set_ylabel("Percentage of Monomer Remaining", labelpad=5, fontsize = 9)
					#adjust y axis limiys
					subplot.set_ylim([0,1])
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
				elif graphType == "Monomer Occurences":
					#adjust axis title
					subplot.set_ylabel("Average Total Monomer Occurences", labelpad=5, fontsize = 9)
					for index in polymerIndex:
						count = 0
						for polymer in polymerArray:
							if polymer[index - 1] == monomer:
								count += 1
						monomercounts[index - 1] = float(float(count) / float(self.numSimulations))
				#debugging purposes
				#print(polymerIndex)
				#print(monomercounts)
				#plots x and y arrays
				curve = subplot.plot(polymerIndex, monomercounts, label = "Monomer " + str(monomer))
			#legend-screw matplotlib; so fucking hard to format
			handles, labels = subplot.get_legend_handles_labels()
			lgd = subplot.legend(handles, labels, prop = {'size':7}, loc = "best")
			subplot.set_xlabel("Monomer Position Index", labelpad = 0, fontsize = 9)
		elif graphType == "Monomer Separation":
			#obtain histogram limit
			try:
				self.histogramLimit = float(self.histogramLimitTkVar.get())
				assert(self.histogramLimit <= 1)
				assert(self.histogramLimit > 0)
			except AssertionError:
				errorMessage("Percent to Analyze Value Invalid!", 220)
				return
			#retrieving all needed variables from inputs
			if number == 1:
				histogramMonomer = int(self.histogramMonomer1TkVar.get())
			if number == 2:
				histogramMonomer = int(self.histogramMonomer2TkVar.get())
			histogramNumberLimit = self.histogramLimit * self.polymerLength
			histogramData = self.getHistogramData(polymerArray, histogramMonomer, histogramNumberLimit)
			#print(histogramData)
			binwidth = 1
			subplot.hist(histogramData, bins=range(min(histogramData), max(histogramData) + binwidth, binwidth),
			 color = self.colorArray[histogramMonomer - 1], normed = True)
			subplot.set_ylabel("Normalized Separation", labelpad=5, fontsize = 9)
			subplot.set_xlabel("Monomer %i Block Size" %histogramMonomer, labelpad = 0, fontsize = 9)
			subplot.set_xticks(arange(min(histogramData), max(histogramData) + 1, 1))
			#print(min(histogramData))
			#print(max(histogramData))
	def saveState(self):
		ratiosList = []
		#Asserts that there are no input errors. Shows errorMessage if error caught.test
		try:
			self.totalMonomers = int(self.totalMonomersTkVar.get())
			monomerAmounts = self.getMonomerAmounts()
			singleCoeffList = self.getCoefficients()
			self.numSimulations = int(self.numSimsTkVar.get())
			self.raftRatio = float(self.raftRatioTkVar.get())
			self.histogramLimit = float(self.histogramLimitTkVar.get())
			#Converts a list of Tkvars to a list of floats
			for entry in self.startingRatiosTkList:
				assert(float(entry.get()) > 0)
				ratiosList.append(float(entry.get()))
			assert(self.histogramLimit <= 1)
			assert(self.histogramLimit > 0)
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
		startingWeightList = []
		#error checking to see if config file exists
		if not os.path.exists("config.txt"):
			errorMessage("config.txt does not exist!", 220)
			return
		nextSetting = LAST_SETTING + 1
		file = open("config.txt", "a")
		file.write("\nSetting %i \nNumber of Unique Monomers = %i " %(nextSetting, self.numMonomers))
		monomerIndex = 1
		while monomerIndex <= self.numMonomers:
			file.write("\nMonomer %i Ratio = %i" %(monomerIndex, ratiosList[monomerIndex - 1]))
			monomerIndex += 1
		monomerIndex = 1
		print("single:", singleCoeffList)
		while monomerIndex <= self.numMonomers:
			innerIndex = 1
			while innerIndex <= self.numMonomers:
				coefflist = singleCoeffList[monomerIndex - 1]
				file.write("\n%i-%i = %f" %(monomerIndex, innerIndex, (coefflist[innerIndex - 1])))
				innerIndex += 1
			monomerIndex += 1
		file.write("\nend")
		file.close()
		global LAST_SETTING
		LAST_SETTING += 1
		infoMessage("Save Successful", "State successfully saved into setting %i!" %nextSetting, 300)
		#self.saveButton.update()
		return

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
	exitButton = ttk.Button(master = top, text = "Ok", command = top.destroy, width = 7)
	exitButton.pack(side = Tk.TOP, pady = 5)
def infoMessage(title, message, width):
	#Toplevel parameters
	top = Tk.Toplevel()
	top.wm_title(title)
	top.geometry("%dx%d%+d%+d" % (width, 70, 250, 125))
	#Message
	msg = Tk.Message(master = top, text = message, width = 500)
	msg.pack(side = Tk.TOP, pady = 5)
	#OK button to exit
	exitButton = ttk.Button(master = top, text = "Ok", command = top.destroy, width = 7)
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
root.wm_title("Compositional Drift %s" % VERSION)
app = Application(master = root)
app.mainloop()
