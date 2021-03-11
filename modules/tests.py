import unittest
from structures import Reaction, Polymer
import traceback
import numpy as np
import test_cases

def generate_random_reaction(num_monomer_species=None, model="Mayo-Lewis", max_monomer_amounts=100000):
    """
    DESCRIPTION: Generates a reaction with randomly chosen inputs
    RETURNS: Reaction object 
    """
    if not num_monomer_species:
        num_monomer_species = np.random.randint(2, 4)
    avg_DP = np.random.randint(1, 500)
    rr1 = np.random.uniform(0, 10)
    rr2 = np.random.uniform(0, 10)
    if model == "Mayo-Lewis":
        reactivity_ratios = np.random.rand(num_monomer_species, num_monomer_species-1)
    elif model == "Penultimate":
        rate_constants = np.random.rand(num_monomer_species, num_monomer_species, num_monomer_species)
    else:
        raise ValueError("did not recognize model")
    conversion = np.random.rand()
    monomer_amounts = np.random.randint(low=100, high=max_monomer_amounts, size=num_monomer_species)
    chain_transfer_probability = np.random.rand()
    reaction = Reaction(num_monomer_species, model)
    reaction.set_monomer_amounts(monomer_amounts)
    if model == "Mayo-Lewis":
        reaction.set_reactivity_ratios(reactivity_ratios)
    elif model == "Penultimate":
        reaction.set_rate_constants(rate_constants)
    reaction.set_average_DP(avg_DP)
    reaction.set_monomer_amounts(monomer_amounts)
    reaction.set_conversion(conversion)
    reaction.set_chain_transfer_probability(chain_transfer_probability)
    return reaction

def test_reaction_self_consistency(TestCase, reaction):
    """
    DESCRIPTION: Given a reaction object, do static checks on the state of the reaction to see if parameters are self-consistent
    """
    # Test that the correct number of polymers have been created and/or maintained
    TestCase.assertEqual(len(reaction.polymer_list), reaction.num_polymers,\
        "\nNumber of polymers initiated is wrong. Expected: {}, Actual: {}\n System: {}".format(reaction.num_polymers, \
            len(reaction.polymer_list), str(reaction)))
    # Make sure that monomer species consumed and current monomer amounts are consistent with one another
    for i in range(reaction.num_monomer_species):
        TestCase.assertEqual(reaction.init_monomer_amounts[i] - reaction.curr_monomer_amounts[i], reaction.monomer_species_consumed[i],\
            "\nDifference in inital and final monomer amounts for monomer {} does not match monomer species consumed.\nReaction: {}"\
                .format(i, reaction))

    # Makes sure that pool size and monomer amounts are consistent with one another
    TestCase.assertEqual(sum(reaction.curr_monomer_amounts), reaction.curr_pool_size,\
        "\nSum of current monomer amounts is not equal to current monomer pool size.\nReaction: {}".format(reaction))
    TestCase.assertEqual(sum(reaction.init_monomer_amounts), reaction.init_pool_size,\
        "\nSum of initial monomer amounts is not equal to initial monomer pool size.\nReaction: {}".format(reaction))

def test_initiation(TestCase, reaction):
    """
    DESCRIPTION: Runs and tests the initiation step
    """
    # Run the initiation, make sure no errors occur
    try:
        reaction.run_initiatation()
    except Exception as e:
        tb = traceback.format_exc()
        TestCase.fail("\nInitiation step failed unexpectedly!\nReaction: {}\n{}".format(reaction, tb))
    # Test that the reaction parameters are self consistent
    test_reaction_self_consistency(TestCase, reaction)

    # Test that all polymers have exactly one monomer
    for polymer in reaction.polymer_list:
        if reaction.model == "Mayo-Lewis":
            TestCase.assertEqual(len(polymer), 1, "\nNot all polymer lengths after initiation are equal to 1.\n \
                Reaction: {}\n Polymer: {}".format(reaction, polymer))
        elif reaction.model == "Penultimate":
            TestCase.assertGreaterEqual(len(polymer), 1, "\nNot all polymer lengths after initiation are equal to 1.\n \
                Reaction: {}\n Polymer: {}".format(reaction, polymer))

    # # Test that the distribution of initial monomers is close to expected
    # pdist = reaction.initiation_distribution()
    # actual_dist = reaction.monomer_distribution()
    # for prob, actual_prob in zip(pdist, actual_dist):
    #     TestCase.assertAlmostEqual(prob, actual_prob, 1, "\nSimulated monomer distribution is too far off. \nExpected: {}\nActual:{}\
    #         \nReaction: {}".format(pdist, actual_dist, reaction))




def test_single_propagation(TestCase, reaction):
    """
    DESCRIPTION: Runs and tests the first propagation step after initiation
    """
    # Run the initiation, make sure no errors occur
    try:
        reaction.run_initiatation()
        reaction.run_single_propagation()
    except Exception as e:
        tb = traceback.format_exc()
        TestCase.fail("\First propagation failed unexpectedly!\nReaction: {}\n{}".format(reaction, tb))
    # Test that the reaction parameters are self consistent
    test_reaction_self_consistency(TestCase, reaction)

def test_complete_reaction(TestCase, reaction):
    """
    DESCRIPTION: Runs and tests the first propagation step after initiation
    """
    # Run the initiation, make sure no errors occur
    try:
        reaction.run_reaction()
    except Exception as e:
        tb = traceback.format_exc()
        TestCase.fail("\First propagation failed unexpectedly!\nReaction: {}\n{}".format(reaction, tb))
    # Test that the reaction parameters are self consistent
    test_reaction_self_consistency(TestCase, reaction)

class Test2MonomerMayoLewis(unittest.TestCase):
    "Test the Mayo-Lewis model for two monomer species"
    def setUp(self):
        self.num_short_tests = 10
        self.num_long_tests = 10
        self.num_edge_tests = 10

    def test_initiation_general_cases(self):
        "Blanket test for errors in the initiation step"
        for i in range(self.num_short_tests):
            with self.subTest(i=i):
                reaction = generate_random_reaction()
                test_initiation(self, reaction)

    def test_single_propagation_general_cases(self):
        "Blanket test for errors in a single propagation step"
        for i in range(self.num_short_tests):
            with self.subTest(i=i):
                reaction = generate_random_reaction()
                test_single_propagation(self, reaction)
        return

    def test_complete_reaction_general_cases(self):
        for i in range(self.num_long_tests):
            with self.subTest(i=i):
                reaction = generate_random_reaction(max_monomer_amounts=1000)
                test_complete_reaction(self, reaction)
        return
    
    def test_penultimate_initiation(self):
        for i in range(self.num_short_tests):
            with self.subTest(i=i):
                reaction = generate_random_reaction(model="Penultimate")
                test_initiation(self, reaction)

    def test_penultimate_single_propagation(self):
        for i in range(self.num_short_tests):
            with self.subTest(i=i):
                reaction = generate_random_reaction(model="Penultimate")
                test_single_propagation(self, reaction)
        return
    def test_penultimate_complete_reaction(self):
        for i in range(self.num_long_tests):
            with self.subTest(i=i):
                reaction = generate_random_reaction(model="Penultimate", max_monomer_amounts=1000)
                test_complete_reaction(self, reaction)
        return
        
    def test_initiation_edge_cases(self):
        for i in range(self.num_edge_tests):
            with self.subTest(i=i):
                for generator in test_cases.edge_case_reaction_generators:
                    reaction = generator()
                    test_initiation(self, reaction)

if __name__ == '__main__':
    unittest.main()