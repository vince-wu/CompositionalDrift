#All neccesary imports
import matplotlib
import random
matplotlib.use('TkAgg')
from numpy import arange, sin, pi
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
# implement the default mpl key bindings
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import sys
if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk
    """"
    Created by Vincent Wu on 8/17/16:
    This program uses the Mayo-Lewis equation and Monte Carlo method to simulate copolymer growth, and
    represents results both graphically and visually
    """
#Main class 
class Application(Tk.Frame):
	def __init__(self, master = None):
		Tk.Frame.__init__(self, master)
		self.pack()
		#GUI Widget Initiation
		#Creates the init screen
		self.initScreen()
		#Creates Input Widgets
		self.createInputWidgets()
		#Creates dummy visualization Frame
		self.visualizationFrame = Tk.Frame(master = root)
		self.destroyCanvas = False
		self.destroyHide = False
	#Destroys unneccesary widgets
	def destroy(self):
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
			#print(numMonomers) # debugging purposes
			self.createMoreInputs(self.numMonomers) 
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
	def createMoreInputs(self, numMonomers):
		#Destroys or edits current widgets
		self.initFrame.destroy()
		self.monomerCount.destroy()
		self.monomerCountLabel.config(text = "Number of Unique Monomers: " + str(numMonomers))
		self.countConfirm.destroy()
		#sets a class instance varaible for numMonomers
		self.numMonomers = numMonomers
		#Commands
		# Quit command: quits window
		def _quit():
			root.quit()
			root.destroy()
		# Back Command: goes back to numMonomers Entry
		def back(self):
			#Destroys all neccesary widgets
			self.destroy()
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
		self.numSimulations = Tk.Entry(master = self.numSimsFrame, width = 5)
		self.numSimulations.pack(side = Tk.LEFT, padx = 3, pady = 3)
		#Setting number of simulations to 1000
		self.numSims = Tk.IntVar()
		self.numSimulations["textvariable"] = self.numSims
		self.numSims.set(1000)
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
		# A list of Tk.Entry objects for Monomer Amount
		self.startingAmountList = [] 
		# A 2D list of Tk.Entry objects for Coefficicients
		self.coefficientList = [] 
		#Frame for Monomer Amounts
		self.amountFrame = Tk.Frame(master = self.inputFrame) 
		self.amountFrame.pack(side = Tk.LEFT, padx = 5)
		#Frame for Monomer Coefficients
		self.coefficientFrame = Tk.Frame(master = self.inputFrame)
		self.coefficientFrame.pack(side = Tk.LEFT, padx = 5)
		#While loop creating number of neccesary amount Entry boxes
		while createCount < numMonomers:
			#Label for inputAmount
			monomerAmountFrame = Tk.Frame(master = self.amountFrame)
			monomerAmountFrame.pack(side = Tk.TOP, padx = 5, pady = 3)
			inputAmountLabel = Tk.Label(master = monomerAmountFrame, text = "     Monomer " 
				+ str(createCount + 1) + " Amount:")
			inputAmountLabel.pack(side = Tk.LEFT)
			#Entry for inputAmount
			inputAmount = Tk.Entry(master = monomerAmountFrame, width = 5)
			inputAmount.pack(side = Tk.LEFT, padx = 5)
			#Setting Default Value to 20
			amount = Tk.IntVar()
			inputAmount["textvariable"] = amount
			amount.set(20)
			#Add Tk.Entry object to startingAmountList
			self.startingAmountList.append(inputAmount)
			createCount += 1
		#Debugging purposes
		#print("startingAmountList: ", self.startingAmountList) 
		createCount2 = 0
		#While loop creating number of neccesary coefficient Entry boxes
		while createCount2 < numMonomers:
			combinations = 0
			#Appends to coefficient list a list containing coefficients for the polymer index
			singleMonoCoeffList = []
			self.coefficientList.append(singleMonoCoeffList)
			#Frame for Coefficients for Single Monomer
			singleCoeffFrame = Tk.Frame(master = self.coefficientFrame)
			singleCoeffFrame.pack(side = Tk.LEFT, fill = Tk.X, expand = 1)
			while combinations < numMonomers:				
				#Label for inputAmount
				coeffValFrame = Tk.Frame(master = singleCoeffFrame)
				coeffValFrame.pack(side = Tk.TOP, padx = 5, pady = 3)
				inputCoeffLabel = Tk.Label(master = coeffValFrame, text = str(createCount2 + 1)
				 + "-" + str(combinations + 1) + " Constant:" )
				inputCoeffLabel.pack(side = Tk.LEFT)
				#Entry for inputAmount
				inputCoeff = Tk.Entry(master = coeffValFrame, width = 3)
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
		 text = "Compositional Drift Simulator v1.1")
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
		#Asserts that there are no input errors. Shows errorMessage if error caught.
		try:
			monomerAmounts = self.getMonomerAmounts()
			singleCoeffList = self.getCoefficients()
			simTest = int(self.numSims.get())
		except ValueError:
			errorMessage("Please input valid parameters!")
			return
		except notInEuropeError:
			errorMessage("You are not in Europe!")
			return
		#destroys hideButton if necessary
		if self.destroyHide:
			self.hideButton.destroy()
		self.destroyHide = True
		#hideButton creation
		self.hideButton = Tk.Button(master = self.buttonFrame, bg = "pale turquoise", activebackground = "light slate blue",
			width = 27, text = "Hide Input Paramters", command = self.hideInputParams)
		self.hideButton.pack(side = Tk.TOP, pady = 3)
		#destroys canvas if necessary
		if self.destroyCanvas:
			self.canvas.get_tk_widget().destroy()
			self.toolbar.destroy()
		self.destroyCanvas = True
		print("monomerAmounts: ", monomerAmounts)
		print("singleCoeffList: ", singleCoeffList)
		#An array of polymers
		polymerArray = []
		counter = 0
		#print(self.numMonomers)
		#Builds number of polymers equal to self.numSims
		while counter < self.numSims.get():
			#Saving original monomerAmounts so resets are possible by copying monomerAmounts
			originalMonomerAmounts = monomerAmounts[:]
			#An array representing a polymer chain
			polymer = []
			#initiating polymer chain
			"""Initiation: Iterates through all monomers, calculating weight chance to initiate for each,
			 and initiates monomer with highest weight chance"""
			count = 1
			largestChance = 0
			total = sum(monomerAmounts)
			choices = []
			while count <= self.numMonomers:
				#weight chance of monomer initation: (amount of starting monomer)
				weight = monomerAmounts[count - 1]
				#Adds a two element list to choices containing monomer and weight
				choices.append([count, weight])
				count += 1
			#Using weighted_choice, selects first monomer
			startingMonomer = weighted_choice(choices)
			#Uses up one monomer of the startingMonomer
			monomerAmounts[startingMonomer - 1] -= 1
			#startingPolymer becomes first monomer in polymer chain
			polymer.append(startingMonomer)
			#building polymer chain
			while sum(monomerAmounts) > 0:
				polyCounter = 1
				largestChance = 0
				choices = []
				"""Propogation: Iterates through all monomers, calculating weight chance to bind for each, 
				and binds monomer with highest weight chance"""
				while polyCounter <= self.numMonomers:
					#Retrieveing coefficient based on previous and current monomer
					coeff = singleCoeffList[polymer[-1] - 1][polyCounter - 1]
					# weight chance calulations for monomer attaching: coefficient*(amount of monomer remaining)
					chance = monomerAmounts[polyCounter - 1] * coeff
					#Adds a two element list to choices containing monomer and weight
					choices.append([polyCounter, chance])
					polyCounter += 1
				#Using weighted_choice, selects next monomer
				nextMonomer = weighted_choice(choices)
				#Reduces number of nextMonomer by 1, since it is being used up in reaction
				monomerAmounts[nextMonomer - 1] -= 1
				#Attaches next monomer to polymer chain
				polymer.append(nextMonomer)
			#Adds finished polymer to list of polymers
			polymerArray.append(polymer)
			#Resets the monomerAmounts to originalMonomerAmounts
			monomerAmounts = originalMonomerAmounts
			counter += 1
		#Debugging purposes
		print("Array of Polymers: ", polymerArray)
		self.visualizationFrame.destroy()
		self.visualizePolymers(polymerArray)
		self.plotCompositions(polymerArray)
		#self.inputFrame.pack_forget()
	#plots compositions given a PolymerArray
	def plotCompositions(self, polymerArray):
		self.lineColors = []
		#Plot and Figure formatting
		self.plotFigure = Figure(figsize=(5.5, 3.3), dpi=100)
		frequencyPlot = self.plotFigure.add_subplot(111)
		frequencyPlot.tick_params(labelsize = 7)
		frequencyPlot.set_ylabel("Total Monomer Occurences", labelpad=5, fontsize = 9)
		frequencyPlot.set_xlabel("Monomer Position Index", labelpad = 0, fontsize = 9)
		#Iterates through each unique monomer and plots composition
		for monomer in range(1, self.numMonomers + 1):
			#x-axis array
			polymerIndex = list(range(1, self.totalNumMonomers + 1))
			#y-axis array initation
			monomercounts = [0] * self.totalNumMonomers
			#inputs counts into y-axis array
			for index in polymerIndex:
				count = 0
				for polymer in polymerArray:
					if polymer[index - 1] == monomer:
						count += 1
				monomercounts[index - 1] = count
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
		self.totalNumMonomers = sum(self.getMonomerAmounts())
		numRows = self.numPolyToShow.get()
		#parameters for canvas height and width
		canvasHeight = 120
		canvasWidth = 1000
		#Maximizes size of squares
		if (canvasHeight - 50) / numRows <= (canvasWidth - 50) / self.totalNumMonomers:
			size = (canvasHeight - 50) / numRows
			canvasWidth = self.totalNumMonomers * size + 20
		else:
			size = (canvasWidth - 50) / self.totalNumMonomers
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


	#Converts an array of Tk.Entrys for numMonomers into an int array
	def getMonomerAmounts(self):
		numStartingAmtList = []
		for entry in self.startingAmountList:
			if entry == 0:
				numStartingAmtList.append(entry)
			else:
				numStartingAmtList.append(int(entry.get()))
		return numStartingAmtList

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
						singleCoeffList.append(float(innerEntry.get()))
		return coeffList
#When called, makes a pop out error informing user of invalid inputs
def errorMessage(message):
	top = Tk.Toplevel()
	top.wm_title("Error")
	top.geometry("%dx%d%+d%+d" % (220, 70, 250, 125))
	msg = Tk.Message(master = top, text = message, width = 400)
	msg.pack(side = Tk.TOP, pady = 5)
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
class notInEuropeError(Exception):
	def __init__(self, value):
		self.value = value
root = Tk.Tk()
root.wm_title("Compositional Drift")
app = Application(master = root)
app.mainloop()
