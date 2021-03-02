import numpy as np
import random

class Polymer:

    def __init__(self, polymerList, model, rateConstantList):
        self.polymerList = polymerList
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
            return self.rateConstantList[self.lastMonomer()][monomer]

        elif self.model == "Penultimate":
            return self.rateConstantList[self.secondToLastMonomer()][self.lastMonomer()][monomer]
    def __str__(self):
        return str(self.polymerList)
    def __repr__(self):
        return self.polymerList

class Reaction():
    def __init__(self, num_monomer_species):
        self.num_monomer_species = num_monomer_species
        self.model = None
        self.average_DP = None
        self.conversion = None
        self.init_monomer_amounts = None
        self.curr_monomer_amounts = None
        self.reactivity_ratios = None
        self.rate_constants = None
        self.init_pool_size = None
        self.curr_pool_size = None
        self.num_polymers = None
        self.chain_transfer_probability = None
        self.hold_composition = False
        # list to keep track of how much of each monomer species is consumed
        self.monomer_species_consumed = [0 for i in range(num_monomer_species)]
        # a list of Polymer objects to keep track of all the polymers in the reaction 
        self.polymer_list = []
    def set_model(self, model): 
        self.model = model
    def set_average_DP(self, average_DP):
        self.average_DP = average_DP
        if self.init_pool_size:
            self.num_polymers = int(self.init_pool_size/self.average_DP)
    def set_conversion(self, conversion):
        self.conversion = conversion
    def set_monomer_amounts(self, monomer_amounts):
        self.init_monomer_amounts = monomer_amounts
        self.curr_monomer_amounts = monomer_amounts
        self.init_pool_size = sum(monomer_amounts)
        self.curr_pool_size = sum(monomer_amounts)
        if self.average_DP:
            self.num_polymers = int(self.init_pool_size/self.average_DP)
    def set_reactivity_ratios(self, reactivity_ratios):
        self.reactivity_ratios = reactivity_ratios
    def set_chain_transfer_probability(self, chain_transfer_probability):
        self.chain_transfer_probability = chain_transfer_probability
    def monomer_distribution(self):

        dist = np.array([sum(polymer.count(monomer) for polymer in self.polymer_list) for monomer in range(self.num_monomer_species)])
        dist = dist/sum(dist)
        
        return dist
    
    def add_monomer(self, monomer, polymer):
        """
        DESCRIPTION: Adds a single monomer species to the specified polymer. Will update attributes of the Reaction object
                     such as curr_monomer_amounts and curr_pool_size accordingly
        PARAMETERS:
            monomer: integer
                An integer representing the monomer species to add to the polymer. Note that monomer indexes start at 0
            polymer: Polymer object
                The polymer to add the monomer to. The polymer object should be part of self.polymer_list, else errors may occur
        RETURNS: N/A
        """
        polymer.append(monomer)
        if not self.hold_composition:
            self.curr_monomer_amounts[monomer] -= 1
            self.curr_pool_size -= 1
        self.monomer_species_consumed[monomer] -= 1
    
    def initiation_distribution(self):
        """
        DESCRIPTION: returns an array representing the probability distribution of which monomer species intiates the reaction
                     depending on the model and number of monomer species, a different initiation calculation will be done
        PARAMETERS: N/A
        RETURNS: 1D array of floats 
                 [P(1), P(2), .... P(N)], where P(x) is the probability of monomer species x initiating the reaction,
                 and N is the number of monomer species. This is a normalized probability distribution
        """

        if self.num_monomer_species == 2 and self.model == "Mayo-Lewis":
            "Case for 2-monomer system: use the Mayo Lewis Equation"
            f1 = self.init_monomer_amounts[0]
            f2 = self.init_monomer_amounts[1]
            total_f = f1+f2
            f1 = f1/total_f
            f2 = f2/total_f
            r1 = self.reactivity_ratios[0][1]
            r2 = self.reactivity_ratios[1][0]
            weight = (r1*f1**2 + f1*f2) / (r1*f1**2 + 2*f1*f2 + r2*f2**2)
            #print("weight: ", weight)pyt
            init_probability_dist = [weight, 1-weight]

        elif self.num_monomer_species == 3 and self.model == "Mayo-Lewis":
            "Case for 3-monomer system: use an altered Mayo Lewis Equation"
            m1 = self.monomer_amounts[0]
            m2 = self.monomer_amounts[1]
            m3 = self.monomer_amounts[2]
            F = m1 + m2 + m3
            f1 = m1/F
            f2 = m2/F
            f3 = m3/F
            r11 = self.rate_constants[0][0]
            r12 = self.rate_constants[0][1]
            r13 = self.rate_constants[0][2]
            r21 = self.rate_constants[1][0]
            r22 = self.rate_constants[1][1]
            r23 = self.rate_constants[1][2]
            r31 = self.rate_constants[2][0]
            r32 = self.rate_constants[2][1]
            r33 = self.rate_constants[2][2]
            R1 = r11 + r12 + r13
            R2 = r21 + r22 + r23
            R3 = r31 + r32 + r33
            a = f1*r11*f1/(r11*f1+r12*f2+r13*f3) + f2*r21*f1/(r21*f1+r22*f2+r23*f3) + f3*r31*f1/(r31*f1+r32*f2+r33*f3)
            b = f1*r12*f2/(r11*f1+r12*f2+r13*f3) + f2*r22*f2/(r21*f1+r22*f2+r23*f3) + f3*r32*f2/(r31*f1+r32*f2+r33*f3)
            c = 1 - a - b
            init_probability_dist = [a, b, c]
        else:
            "Case for any monomer system > 3: initial solely based on feed mole ratios 'f'"

            "Also Case for any Penultimate system"

            "Cycle through each monomer, finding its initial amount and using that value as the weight."
            init_probability_dist = [self.monomer_amounts[i]/self.init_pool_size for i in range(self.num_monomer_species)]
        return init_probability_dist

    def run_initiatation(self):
        """
        DESCRIPTION: Initiates the polymerization. For each polymer, randomly select the starting monomer(s) from a probability
                     distribution that depends on the ratios of each starting monomer species, and (for 2- or 3- monomer Mayo-Lewis models),
                     on the reactivity ratios. By the end of this step, all polymers will have either one (model=Mayo Lewis) or
                     two (model=Penultimate) monomers in their chain
        PARAMETERS: N/A
        RETURNS: N/A
        """
        # For both Mayo-Lewis and Penultimate models, initiate the first monomer 
        for i in range(self.num_polymers):
            # Get the probability distribution of how likely each monomer species is initiate the chain. 
            probability_dist = self.initiation_distribution()
            # print(probability_dist)
            # Randomly choose a monomer using the probability distribution
            starting_monomer = np.random.choice(range(self.num_monomer_species), p=probability_dist)
            # Initiate a polymer object
            polymer = Polymer([], self.model, self.rate_constants)
            # Add the selected starting monomer to the polymer
            self.add_monomer(starting_monomer, polymer)
            # Add the polymer to the polymer list
            self.polymer_list.append(polymer)

        # For Penultimate Model, the second monomer in chain must also be chosen before propagation, based off of 
	    # HETEROGENOUS constants (works only for 2 monomer system)
        if self.model == "Penultimate":
            #iterating through each polymer in batch
            for polymer in self.polymer_list:
                probability_dist = np.array([])
                #calculate weight chance for each monomer species
                for i in range(self.num_monomer_species):
                    #retrieving coefficient based on previous and currrent monomer, assumes heterogoenous penultimate monomer
                    if polymer.lastMonomer() == 0:
                        penultimate_monomer = 1
                    else:
                        penultimate_monomer = 0
                    rr = self.rateConstantList[penultimate_monomer][polymer.lastMonomer()-1][i]
                    # weight chance calulations for monomer attaching: coefficient*(amount of monomer remaining)
                    weight = monomerAmounts[i] * rr
                    #adds a two element list to choices containing monomer and weight
                    np.append(probability_dist, weight)
                # normalize probability distribution
                probability_dist = probability_dist/sum(probability_dist)
                #Using weighted_choice, selects next monomer
                second_monomer = np.random.choice(range(self.num_monomer_species), p=probability_dist)
                #Attach the monomer to the polymer
                self.add_monomer(second_monomer, polymer)    
    
    def run_single_propagation(self):
        """
        DESCRIPTION: Executes a single propagation step by growing the length on one randomly chosen polymer chain.
                     The number of monomers is a random variable whose distribution depends on the chain transfer percentage
                     The monomer species chosen to be added to the chain is randomly selected from a distribution that is directly proportional
                     to both the amount of the monomer species remaining in the pool, and the rate constant associated to that species
                     This method also changes the state of the Reaction object, by updating attributes such as curr_monomer_amounts and 
                     curr_polymer_array
        PARAMETERS: N/A
        OUTPUT: Polymer object 
                The polymer chain that was grown. 
        """
        # Randomly choose a polymer chain to grow.
        polymer = random.choice(self.polymerArray)

        # For the polymer selected, iterate through all possible monomers which can be added. For each monomer, calculate
        # the weight chance of the monomer to be added to the chain , defined as the product of the relevant rate 
        # constant 'k' times the number of monomers left unreacted 'f'

        # A variable to keep track of number of monomers to add. This is usually 1; however if the Chain Transfer % is less than 100, there
        # is a chance that more than one monomer is added.
        num_monomers_to_add = 1

        # change that an extra monomer will be added in this single propagation step
        fudge_factor = 1 - self.chain_transfer_probability/100

        # calculate how many monomers to add (usually just 1, but with nonzero chain transfer percentage this can be more than one)
        while True:
            # If number of monomers added exceeds total monomers left to consume do not increase number of monomers to add
            # curr_pool_size is not used here because if composition of pool is maintained, curr_pool_size would never decrease
            if sum(self.monomer_species_consumed) + num_monomers_to_add >= self.init_pool_size*self.conversion:
                break
            # Randomly determine if another monomer will be added
            elif random.random() >= fudge_factor:
                break
            # keep adding another monomer, until the random selector stops
            else:
                num_monomers_to_add += 1

        # the probability for each monomer to be selected as the next monomer to be added
        growth_probability_dist = np.array([])
        for i in range(num_monomers_to_add):
            # A list representing the probability distribution of each monomer species being added to the polymer chain
            for monomer in range(self.num_monomer_species):
                # Retrieve relevant rate constant k based on previous and current monomer
                k = polymer.rateConstant(monomer)
                # weight chance calulations for monomer attaching: coefficient*(amount of monomer remaining)
                weight = self.monomer_amounts[monomer] * k
                np.append(growth_probability_dist, weight)

        # If all weights are zero due to coefficients all being zero, then sort by relative amounts of monomer instead
        if sum(growth_probability_dist) == 0:
            growth_probability_dist = np.array([self.monomer_amounts[i] for i in range(self.num_monomer_species)])
        
        # normalize the probabilty distribution
        growth_probability_dist = growth_probability_dist/sum(growth_probability_dist)

        # select next monomer to be appended to the growing chain based on the probability distribution
        next_monomer = np.random.choice(range(self.num_monomer_species), p=probability_dist)

        # Attach next monomer to polymer chain
        self.add_monomer(next_monomer, polymer)

        # return the new polymer
        return polymer
    
    def run_complete_propagation(self):
        """
        DESCRIPTION: Runs the entire propagation reaction until termination, which occurs when the number of monomers consumed
                     by the reaction equals the number of monomers in the initial pool. 
        PARAMETERS: N/A
        RETURNS: N/A
        """
        # repeat a single propagation until the number of monomers consumed equal the initial number of monomers
        # curr_pool_size is not used here because if composition of pool is maintained, curr_pool_size would never decrease
        while sum(self.monomer_species_consumed) < self.init_pool_size*self.conversion:
            self.run_single_propagation()
        return

    def __repr__(self):
        repr = """
Model: {}
Number of monomer species: {}
Dispersity: {}
Reaction conversion: {}
Number of polymers to synthesize: {}
Reactivity ratios: {}
Chain Transfer Probability: {}
Initial monomer species amounts: {}
Current monomer species amounts: {}
Initial monomer pool size: {}
Current monomer pool size: {}
Monomers consumed: {}
        """.format(self.model, self.num_monomer_species, self.average_DP,  self.conversion, self.num_polymers,
        self.reactivity_ratios, self.chain_transfer_probability, self.init_monomer_amounts, self.curr_monomer_amounts,
        self.init_pool_size, self.curr_pool_size, self.monomer_species_consumed)
        return repr

