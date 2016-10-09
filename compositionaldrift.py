#All neccesary imports
import matplotlib
import random
import math
import json
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
#Main class 
class Application(Tk.Frame):
	def __init__(self, master = None):
		Tk.Frame.__init__(self, master)
		self.pack()
		self.initialize()
	#Initialization
	def initialize(self):
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
		#Confirms number of monomers, creates more input widgets
		def enter(self):
			#nonlocal numMonomers 
			self.numMonomers = int(self.monomerCount.get())
			#asserting that input is in correct range
			assert self.numMonomers < 8
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
		#A frame for all buttons in inputFrame
		self.buttonFrame = Tk.Frame(master = self.inputFrame)
		self.buttonFrame.pack(side = Tk.LEFT, padx = 5)
		#A frame for spinbox and label
		countFrame = Tk.Frame(master = self.buttonFrame)
		countFrame.pack(side = Tk.TOP, padx = 5, pady = 0)
		#Label for spinbox
		self.monomerCountLabel = Tk.Label(master = countFrame, text = "Number of Unique Monomers:")
		self.monomerCountLabel.pack(side = Tk.LEFT, padx = 0, pady = 0)
		#MonomerCount spinbox
		self.monomerCount = Tk.Spinbox(master = countFrame, from_ = 1, to = 7, width = 2)
		self.monomerCount.pack(side = Tk.LEFT, padx = 5, pady = 0)
		#countConfirm Button
		self.countConfirm = Tk.Button(master = self.inputFrame, text = "Enter",
		 command = lambda:enter(self), bg = "light blue", activebackground = "light slate blue", width = 9)
		self.countConfirm.pack(side = Tk.LEFT, padx = 5, pady = 5)
	#Creates more input widgets based on numMonomers
	def createMoreInputs(self):
		#Destroys or edits current widgets
		self.initFrame.destroy()
		self.monomerCount.destroy()
		self.monomerCountLabel.config(text = "Number of Unique Monomers: " + str(self.numMonomers))
		self.countConfirm.destroy()
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
		self.numSimsTkVar.set(200)
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
		self.numPolyToShow.set(8)
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
		self.graphTypeTkIntVar.set(0)
		#RadioButtons for display type
		self.percentageRadioButton = Tk.Radiobutton(master = self.initialConditionsFrame, text = "Graph Monomer Occurences",
		 variable = self.graphTypeTkIntVar, value = 0)
		self.occurenceRadioButton = Tk.Radiobutton(master = self.initialConditionsFrame, text = "Graph Percentage Monomer",
		 variable = self.graphTypeTkIntVar, value = 1)
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
		self.totalMonomersTkVar.set(1000)
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
		self.raftRatioTkVar.set(0.01)
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
						assert(float(innerEntry.get()) > 0)
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
