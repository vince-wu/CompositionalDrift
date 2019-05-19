"Polymer Object"

class Polymer:

	def __init__(self, polymerList, numMonomers, rateConstantList, model):
		self.polymerList = polymerList
		self.numMonomers = numMonomers
		self.rateConstantList = rateConstantList
		self.model = model

	def lastMonomer(self):
		return self.polymerList[-1]

	def secondToLastMonomer(self):
		return self.polymerList[-2]

	def asArray(self):
		return self.polymerList

	def len(self):
		return len(self.polymerList)

	def append(self, monomer):
		self.polymerList.append(monomer)

	def rateConstant(self, monomer):

		if self.model == "Mayo-Lewis":
			return self.rateConstantList[self.lastMonomer()-1][monomer-1]

		elif self.model == "Penultimate":
			return self.rateConstantList[self.secondToLastMonomer()-1][self.lastMonomer()-1][monomer-1]











