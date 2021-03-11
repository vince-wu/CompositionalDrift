import numpy as np
import random
import math
import errors
import numbers

class Polymer(list):
    def __init__(self, polymer_list, model, rate_constants):
        list.__init__(self, polymer_list)
        self.rate_constants = rate_constants
        self.model = model

    def last_monomer(self):
        return self[-1]

    def second_to_last_monomer(self):
        return self[-2]

    def rate_constant(self, monomer):
        if self.model == "Mayo-Lewis":
            return self.rate_constants[self.last_monomer()][monomer]

        elif self.model == "Penultimate":
            return self.rate_constants[self.second_to_last_monomer()][self.last_monomer()][monomer]

class Reaction():
    def __init__(self, num_monomer_species, model):
        self.num_monomer_species = num_monomer_species
        self.set_model(model)
        self.average_DP = None
        self.conversion = 1
        self.init_monomer_amounts = None
        self.curr_monomer_amounts = None
        self.reactivity_ratios = None
        self.rate_constants = None
        self.init_pool_size = None
        self.curr_pool_size = None
        self.num_polymers = None
        self.chain_transfer_probability = 1
        self.hold_composition = False
        # list to keep track of how much of each monomer species is consumed
        self.monomer_species_consumed = [0 for i in range(num_monomer_species)]
        # a list of Polymer objects to keep track of all the polymers in the reaction 
        self.polymer_list = []

    def set_model(self, model): 
        """
        DESCRIPTION: Sets the reaction model
        PARAMETERS: 
            average_DP: string
        """
        if not isinstance(model, str):
            raise ValueError("Model should be a string")
        if model != "Mayo-Lewis" and model != "Penultimate":
            raise ValueError("Did not recognize model. Model should either be 'Mayo-Lewis' or 'Penultimate'")
        self.model = model

    def set_average_DP(self, average_DP):
        """
        DESCRIPTION: Sets the reaction's average dispersity (polymer length)
        PARAMETERS: 
            average_DP: int
        """
        if not isinstance(average_DP, numbers.Integral):
            raise ValueError("Dispersity should be an integer")
        self.average_DP = average_DP
        if self.init_pool_size:
            self.num_polymers = int(self.init_pool_size/self.average_DP)

    def set_conversion(self, conversion):
        """
        DESCRIPTION: Sets the reaction conversion
        PARAMETERS: 
            average_DP: float betwee 0 and 1
        """
        if not isinstance(conversion, numbers.Number):
            raise ValueError("Conversion should be a number between (0, 1]")
        if conversion <= 0 or conversion > 1:
            raise ValueError("Conversion should be a number between (0, 1]")

        self.conversion = conversion

    def set_monomer_amounts(self, monomer_amounts):
        """
        DESCRIPTION: Sets the initial and current amount of monomers in the pool for each monomer species
        PARAMETERS: 
            average_DP: array of size Reaction.num_monomer_species
        """
        if not hasattr(monomer_amounts, '__iter__'): 
            raise ValueError("Monomer amounts should be a list of size Reaction.num_monomer_species")
        if len(monomer_amounts) != self.num_monomer_species:
            raise ValueError("Monomer amounts should be a list of size Reaction.num_monomer_species")
        if not all([isinstance(val, numbers.Integral) and val >= 0 for val in monomer_amounts]):
            raise ValueError("Each monomer amount needs to be a nonegative integer")

        self.init_monomer_amounts = monomer_amounts.copy()
        self.curr_monomer_amounts = monomer_amounts
        self.init_pool_size = sum(self.init_monomer_amounts)
        self.curr_pool_size = sum(monomer_amounts)
        if self.average_DP:
            self.num_polymers = int(self.init_pool_size/self.average_DP)

    def set_reactivity_ratios(self, reactivity_ratios):
        """
        DESCRIPTION: Sets the reactivity ratios for the reaction. Only used for the Mayo-Lewis model. The reactivity ratio 
                     is a matrix of shape (n, n-1), where n = number of monomer species. The value in the ith column and jth row 
                     is the reactivity ratio corresponding to monomer j appending onto a polymer terminated by a monomer i. 
                     If the inputted reactivity ratio is of shpae (n-1, n), it is still valid, and this function will automatically
                     transpose the matrix into the correct shape (n, n-1). The rate constants, which 
        PARAMETERS: 2D array of shape (n, n-1) OR (n-1, n), where n = number of monomer species
        RETURNS: N/A
        """
        if self.model == "Penultimate":
            raise ValueError("Can not set reactivity ratios for Penultimate model. Use Reaction.set_rate_constants() instead")
        if not hasattr(reactivity_ratios, '__iter__'): 
            raise ValueError("Reactivity ratios should be array-like")
        shape = np.shape(reactivity_ratios)
        n = self.num_monomer_species
        if self.num_monomer_species != 2:
            if shape != (n,n-1) and shape != (n-1, n):
                raise(ValueError("Reactivity ratios should be of shape (n, n-1) or (n, n-1) for > 2 monomer system"))
        if self.num_monomer_species == 2:
            if shape == (2,):
                reactivity_ratios = [[reactivity_ratios[0]], [reactivity_ratios[1]]]
            elif shape != (1,2) and shape != (2,1):
                raise ValueError("For a 2-monomer system, reactivity ratios should be a 2-element list, or a matrix of shape (2,1) or (1,2)")
        if not all([isinstance(val, numbers.Number) and val >= 0 for rows in reactivity_ratios for val in rows]):
            raise ValueError("All reactivity ratios should be nonnegative numbers")
        dim1, dim2 = np.shape(reactivity_ratios)
        if dim1 < dim2:
            reactivity_ratios = np.array(reactivity_ratios).T
        self.reactivity_ratios = reactivity_ratios
        self.set_rate_constants(self.rr_to_rate_constants(reactivity_ratios))
    
    def set_rate_constants(self, rate_constants):
        """
        DESCRIPTION: Sets rate constants of the reaction
        PARAMETERS: 
            rate_constants: 
                Model = Mayo Lewis: a square nxn matrix (rank 2 tensor), where n = Reaction.num_monomer_species
                Model = Penultimate: a rank 3 tensor of size n
        """
        if not hasattr(rate_constants, '__iter__'): 
            raise ValueError("Rate constants should be array-like")
        shape = np.shape(rate_constants)
        n = self.num_monomer_species
        if self.model=="Mayo-Lewis" and shape != (n, n):
            raise ValueError("Rate constants should either be a square matrix of size n, where n = Reaction.num_monomer_species")
        elif self.model=="Penultimate" and shape != (n, n, n):
            raise ValueError("Rate constants should be a tensor of rank n, where n = Reaction.num_monomer_species")
        self.rate_constants = rate_constants
        return

    def set_chain_transfer_probability(self, chain_transfer_probability):
        """
        DESCRIPTION: Sets chain transfer probability
        PARAMETERS: 
            chain_transfer_probability: int, between [0,1]
        """
        if not isinstance(chain_transfer_probability, numbers.Number):
            raise ValueError("Chain transfer probability should be a number between (0, 1]")
        if chain_transfer_probability <= 0 or chain_transfer_probability > 1:
            raise ValueError("Chain transfer probability should be a number between [0, 1], inclusive")
        self.chain_transfer_probability = chain_transfer_probability

    def set_hold_composition(self, hold_composition):
        """
        DESCRIPTION: Sets whether or not to hold composition (do not take monomers away from pool during polymerization)
        PARAMETERS: 
            hold_composition: boolean
        """
        if not isinstance(hold_composition, bool):
            raise ValueError("Hold composition should be a boolean value")
        self.hold_composition = hold_composition

    def monomer_distribution(self):
        """
        DESCRIPTION: Calculates the normalized distribution of each monomer species in all monomers across current polymers in the reaction
        PARAMETERS: n/a
        RETURNS: list
            A list of size Reaction.num_monomer_species, where the ith value in the list is the decimal percentage of monomers in all polymers in 
            Reaction.polymer_list that is of species 'i'. The distribution is normalized; i.e the list sums to 1
        """
        dist = np.array([sum(polymer.count(monomer) for polymer in self.polymer_list) for monomer in range(self.num_monomer_species)])
        dist = dist/sum(dist)   
        return dist

    def rr_to_rate_constants(self, reactivity_ratios):
        """
        DESCRIPTION: Converts reactivity ratios into rate constants. Only used for Mayo-Lewis model. The rate constants are represented
                     as a square matrix of size n, where n is the number of monomer species. The rate constant is inversly related to its 
                     corresponding reactivity ratio, r12 = k11/k12. Since homopolymerization rate constants are assumed to
                     all be equal to 1; that is, k11 = k22 = k33 = ... = 1, the rate constant is simply the inverse of the corresponding
                     reactivity ratio. The methods used by the Reaction object use rate constants rather than reactivity ratios
                     to calculate probablities of a monomer appending to a growing polymer chain, as it is more convienient. 

        PARAMETERS: 2D array of shape (n, n-1), where n = number of monomer species
                    A 2D array representing the reactivity ratios for the monomer species. Since homopolymerization constants are assumed
                    to be unity, they are not included in the 2D array

        RETURNS:  Square matrix (rank 2 tensor) of size n, where n = number of monomer species.
                  A tensor representing the reactivity ratios of monomers appending to a polymer
        """
        n = self.num_monomer_species
        if self.model == "Mayo-Lewis":
            rate_constants = np.empty([n, n])
            for i in range(n):
                j_index = 0
                for j in range(n):
                    # for homopolymerization terms, always set rate constant equal to 1
                    if i == j:
                        rate_constants[i][j] = 1
                    # for heteronuclear terms, set rate constant to the inverse
                    else:
                        rr = reactivity_ratios[i][j_index]
                        if rr == 0:
                            rate_constants[i][j] = math.inf
                        else:
                            rate_constants[i][j] = 1/rr
                        j_index += 1

        return rate_constants
    
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
        assert self.curr_monomer_amounts[monomer] > 0, "There are no remaining monomer {} species left in the pool".format(monomer)
        polymer.append(monomer)
        if not self.hold_composition:
            self.curr_monomer_amounts[monomer] -= 1
            self.curr_pool_size -= 1
        self.monomer_species_consumed[monomer] += 1

    def validate_parameters(self):
        """
        DESCRIPTION: Ensures that all reaction parameters need are present and are valid. Throws an exception if this is not the case 
        """
        # Check for missing parameters
        if not self.model:
            raise errors.MissingInputError("Model is a required Reaction parameter")
        if not self.num_monomer_species:
            raise errors.MissingInputError("Number of monomer species is a required Reaction parameter")
        if self.average_DP is None:
            raise errors.MissingInputError("Average dispersity is a required Reaction parameter")
        if self.conversion is None:
            raise errors.MissingInputError("Conversion is a required Reaction parameter")
        if self.init_monomer_amounts is None:
            raise errors.MissingInputError("Initial monomer amounts is a required Reaction parameter")
        if self.reactivity_ratios is None and self.rate_constants is None:
            raise errors.MissingInputError("Reactivity ratios or rate constants is a required Reaction parameter")
        if self.chain_transfer_probability is None:
            raise errors.MissingInputError("Chain transfer probabilty is a required Reaction parameter")
        if self.hold_composition is None:
            raise errors.MissingInputError("Hold composition is a required Reaction parameter")

        # Check for valid parameters


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
            k1 = self.rate_constants[0][1]
            k2 = self.rate_constants[1][0]
            f1 = self.curr_monomer_amounts[0]
            f2 = self.curr_monomer_amounts[1]
            total_f = f1+f2
            f1 = f1/total_f
            f2 = f2/total_f
            # edge cases for zero-valued rate constants
            # if both rate constants are zero, use ratio of starting monomer species as probability distribution
            if k1 == 0 and k2 == 0:
                init_probability_dist = [f1, f2]
            # if only one rate constant is zero, set that corresponding monomer to have 100% chance of initiating
            elif k1 == 0 or k2 == 0:
                if k1 == 0:
                    init_probability_dist = [1, 0]
                else:
                    init_probability_dist = [0, 1]
            # standard case, use instantaneous Mayo Lewis equation
            else:
                r1 = 1/k1
                r2 = 1/k2
                weight = (r1*f1**2 + f1*f2) / (r1*f1**2 + 2*f1*f2 + r2*f2**2)
                init_probability_dist = [weight, 1-weight]

        elif self.num_monomer_species == 3 and self.model == "Mayo-Lewis":
            "Case for 3-monomer system: use an altered Mayo Lewis Equation"
            m1 = self.curr_monomer_amounts[0]
            m2 = self.curr_monomer_amounts[1]
            m3 = self.curr_monomer_amounts[2]
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
            # case for infinite reactivity ratio or no monomer species remaining, just use starting monomer ratios as weights
            if any([r == math.inf for r in [r11, r12, r13, r21, r22, r23, r31, r32, r33]]) or any([m == 0 for m in [m1, m2, m3]]):
                init_probability_dist = [f1, f2 ,f3]
            else:
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
            init_probability_dist = [self.init_monomer_amounts[i]/self.init_pool_size for i in range(self.num_monomer_species)]
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
        self.validate_parameters()
        # For both Mayo-Lewis and Penultimate models, initiate the first monomer 
        for i in range(self.num_polymers):
            # if there are no monomers left in the pool, stop the reaction
            if self.curr_pool_size <= 0:
                break
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
                # if there are no monomers left in the pool, stop the reaction
                if self.curr_pool_size <= 0:
                    break
                probability_dist = np.array([])
                #calculate weight chance for each monomer species
                for i in range(self.num_monomer_species):
                    #retrieving coefficient based on previous and currrent monomer, assumes heterogoenous penultimate monomer
                    if polymer.last_monomer() == 0:
                        penultimate_monomer = 1
                    else:
                        penultimate_monomer = 0
                    rr = self.rate_constants[penultimate_monomer][polymer.last_monomer()-1][i]
                    # weight chance calulations for monomer attaching: coefficient*(amount of monomer remaining)
                    weight = self.curr_monomer_amounts[i] * rr
                    #adds a two element list to choices containing monomer and weight
                    probability_dist = np.append(probability_dist, weight)
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
        RETURNS: Polymer object 
                The polymer chain that was grown. 
        """
        # Randomly choose a polymer chain to grow.
        polymer = random.choice(self.polymer_list)

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
                k = polymer.rate_constant(monomer)
                # weight chance calulations for monomer attaching: coefficient*(amount of monomer remaining)
                weight = self.curr_monomer_amounts[monomer] * k
                np.append(growth_probability_dist, weight)

        # If all weights are zero due to coefficients all being zero, then sort by relative amounts of monomer instead
        if sum(growth_probability_dist) == 0:
            growth_probability_dist = np.array([self.curr_monomer_amounts[i] for i in range(self.num_monomer_species)])
        
        # normalize the probabilty distribution
        growth_probability_dist = growth_probability_dist/sum(growth_probability_dist)

        # select next monomer to be appended to the growing chain based on the probability distribution
        next_monomer = np.random.choice(range(self.num_monomer_species), p=growth_probability_dist)

        # Attach next monomer to polymer chain
        self.add_monomer(next_monomer, polymer)

        # return the new polymer
        return polymer
    
    def run_complete_propagation(self):
        """
        DESCRIPTION: Runs the entire propagation reaction until termination, which occurs when the number of monomers consumed
                     by the reaction equals the number of monomers in the initial pool multiplied by the conversion. 
        PARAMETERS: N/A
        RETURNS: N/A
        """
        # repeat a single propagation until the number of monomers consumed equal the initial number of monomers
        # curr_pool_size is not used here because if composition of pool is maintained, curr_pool_size would never decrease
        while sum(self.monomer_species_consumed) < self.init_pool_size*self.conversion:
            self.run_single_propagation()
        return
    
    def run_reaction(self):
        """
        DESCRIPTION: Runs the entire polymerization reaction until termination, which occurs when the number of monomers consumed
                     by the reaction equals the number of monomers in the initial pool.  multiplied by the conversion
        PARAMETERS: N/A
        RETURNS: N/A
        """
        self.run_initiatation()
        self.run_complete_propagation()

    def __repr__(self):
        repr = """
Model: {}
Number of monomer species: {}
Dispersity: {}
Reaction conversion: {}
Number of polymers to synthesize: {}
Reactivity ratios: \n{}
Rate constants: \n{}
Chain Transfer Probability: {}
Initial monomer species amounts: {}
Current monomer species amounts: {}
Initial monomer pool size: {}
Current monomer pool size: {}
Monomers consumed: {}
        """.format(self.model, self.num_monomer_species, self.average_DP,  self.conversion, self.num_polymers,
        self.reactivity_ratios, self.rate_constants, self.chain_transfer_probability, self.init_monomer_amounts, 
        self.curr_monomer_amounts, self.init_pool_size, self.curr_pool_size, self.monomer_species_consumed)
        return repr

if __name__ == '__main__':
    polymer = Polymer([1,1,1,0,0,1,1,0], "Mayo-Lewis", [0.1,2])
    num_species = 2
    model = "Mayo-Lewis"
    avg_DP = np.random.randint(1, 500)
    rr1 = np.random.uniform(0, 10)
    rr2 = np.random.uniform(0, 10)
    reactivity_ratios = np.random.rand(num_species, num_species-1)
    reactivity_ratios = [math.inf,0]
    conversion = np.random.rand()
    monomer_amounts = np.random.randint(low=100, high=100000, size=num_species)
    chain_transfer_probability = np.random.rand()
    reaction = Reaction(num_species, model)
    reaction.set_monomer_amounts(monomer_amounts)
    reaction.set_reactivity_ratios(reactivity_ratios)
    reaction.set_average_DP(avg_DP)
    reaction.set_conversion(conversion)
    reaction.set_chain_transfer_probability(chain_transfer_probability)
    reaction.run_initiatation()
    reaction.run_single_propagation()
    print(reaction)
