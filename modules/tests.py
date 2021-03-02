import unittest
from structures import Reaction, Polymer
import traceback
import numpy as np

class Test2MonomerMayoLewis(unittest.TestCase):
    def setUp(self):
        self.num_monomer_species = 2
        self.model = "Mayo-Lewis"
        self.avg_DP = 100
        self.reactivity_ratios = np.random.rand(2,2)
        self.monomer_amounts = np.random.rand(2)
        self.conversion = np.random.rand()
        self.chain_transfer_probability = np.random.rand()
        self.reaction =  Reaction(2)
        self.reaction.set_model("Mayo-Lewis")
        self.reaction.set_average_DP(100)
        self.reaction.set_reactivity_ratios([[1,0.1],[2,1]])
        self.reaction.set_monomer_amounts([100, 200])
        self.reaction.set_conversion(100)
        self.reaction.set_chain_transfer_probability(0)
    def test_general_cases(self):
        "Tests for errors in initiation step"
        num_tests = 1
        for i in range(num_tests):
            with self.subTest(i=i):
                avg_DP = np.random.randint(1, 500)
                rr1 = np.random.uniform(0, 10)
                rr2 = np.random.uniform(0, 10)
                reactivity_ratios = [[1,rr1],[rr2,1]]
                conversion = np.random.rand()
                monomer_amounts = np.random.randint(100000, size=2)
                chain_transfer_probability = np.random.rand()
                reaction = Reaction(2)
                reaction.set_model("Mayo-Lewis")
                reaction.set_monomer_amounts(monomer_amounts)
                reaction.set_reactivity_ratios(reactivity_ratios)
                reaction.set_average_DP(avg_DP)
                reaction.set_monomer_amounts(monomer_amounts)
                reaction.set_conversion(conversion)
                reaction.set_chain_transfer_probability(chain_transfer_probability)
                try:
                    reaction.run_initiatation()
                except Exception as e:
                    tb = traceback.format_exc()
                    self.fail("Initiation step failed unexpectedly!\n Reaction Parameters: {}\n {}".format(reaction, tb))
                self.assertEqual(len(reaction.polymer_list), reaction.num_polymers,\
                    "Number of polymers initiated is wrong. Expected: {}, Actual: {}".format(reaction.num_polymers, \
                        len(reaction.polymer_list)))
                pdist = reaction.initiation_distribution()
                actual_dist = reaction.monomer_distribution()
                self.assertAlmostEqual(pdist, actual_dist, "Simulatated monomer distribution is too far off. Expected: {}, Actual: {}"\
                    .format(pdist, actual_dist))
                for polymer in reaction.polymer_list:
                    self.assertEqual(polymer.len(), 1, "Not all polymer lengths after inititation are equal to 1!\n \
                        system: {}\n Polymer: {}".format(str(reaction), str(polymer)))
                





if __name__ == '__main__':
    unittest.main()