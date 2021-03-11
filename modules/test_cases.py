import math
import numpy as np
from structures import Reaction, Polymer


def some_zero_valued_reactivities():
    """
    DESCRIPTION: Generates a reaction with zeros in the reactivity ratio list
    RETURNS: Reaction object 
    """
    n = np.random.randint(2, 4)
    model = "Mayo-Lewis"
    max_monomer_amounts = 5000
    avg_DP = np.random.randint(1, 500)
    reactivity_ratios = np.random.rand(n, n-1)
    # choose a random number of attempts put zeros in reactivity ratios
    num_zeros = np.random.randint(n*(n-1))
    # put zeros into the reactivity ratio matrix
    for i in range(num_zeros):
        row = np.random.randint(0,n)
        col = np.random.randint(0,n-1)
        reactivity_ratios[row][col] = 0

    conversion = np.random.rand()
    monomer_amounts = np.random.randint(low=100, high=max_monomer_amounts, size=n)
    chain_transfer_probability = np.random.rand()
    reaction = Reaction(n, model)
    reaction.set_monomer_amounts(monomer_amounts)
    reaction.set_reactivity_ratios(reactivity_ratios)
    reaction.set_average_DP(avg_DP)
    reaction.set_monomer_amounts(monomer_amounts)
    reaction.set_conversion(conversion)
    reaction.set_chain_transfer_probability(chain_transfer_probability)
    return reaction

def some_inf_valued_reactivities():
    """
    DESCRIPTION: Generates a reaction with infinities in the reactivity ratio list
    RETURNS: Reaction object 
    """
    n = np.random.randint(2, 4)
    model = "Mayo-Lewis"
    max_monomer_amounts = 5000
    avg_DP = np.random.randint(1, 500)
    reactivity_ratios = np.random.rand(n, n-1)
    num_infs = np.random.randint(n*(n-1))
    for i in range(num_infs):
        row = np.random.randint(0,n)
        col = np.random.randint(0,n-1)
        reactivity_ratios[row][col] = math.inf

    conversion = np.random.rand()
    monomer_amounts = np.random.randint(low=100, high=max_monomer_amounts, size=n)
    chain_transfer_probability = np.random.rand()
    reaction = Reaction(n, model)
    reaction.set_monomer_amounts(monomer_amounts)
    reaction.set_reactivity_ratios(reactivity_ratios)
    reaction.set_average_DP(avg_DP)
    reaction.set_monomer_amounts(monomer_amounts)
    reaction.set_conversion(conversion)
    reaction.set_chain_transfer_probability(chain_transfer_probability)
    return reaction

def all_zero_valued_reactivities():
    """
    DESCRIPTION: Generates a reaction with zeros in the entire reactivity ratio list
    RETURNS: Reaction object 
    """
    n = np.random.randint(2, 4)
    model = "Mayo-Lewis"
    max_monomer_amounts = 5000
    avg_DP = np.random.randint(1, 500)
    reactivity_ratios = np.zeros((n, n-1))
    conversion = np.random.rand()
    monomer_amounts = np.random.randint(low=100, high=max_monomer_amounts, size=n)
    chain_transfer_probability = np.random.rand()
    reaction = Reaction(n, model)
    reaction.set_monomer_amounts(monomer_amounts)
    reaction.set_reactivity_ratios(reactivity_ratios)
    reaction.set_average_DP(avg_DP)
    reaction.set_monomer_amounts(monomer_amounts)
    reaction.set_conversion(conversion)
    reaction.set_chain_transfer_probability(chain_transfer_probability)
    return reaction

def all_inf_valued_reactivities():
    """
    DESCRIPTION: Generates a reaction with zeros in the entire reactivity ratio list
    RETURNS: Reaction object 
    """
    n = np.random.randint(2, 4)
    model = "Mayo-Lewis"
    max_monomer_amounts = 5000
    avg_DP = np.random.randint(1, 500)
    reactivity_ratios = np.full((n, n-1), math.inf)
    conversion = np.random.rand()
    monomer_amounts = np.random.randint(low=100, high=max_monomer_amounts, size=n)
    chain_transfer_probability = np.random.rand()
    reaction = Reaction(n, model)
    reaction.set_monomer_amounts(monomer_amounts)
    reaction.set_reactivity_ratios(reactivity_ratios)
    reaction.set_average_DP(avg_DP)
    reaction.set_monomer_amounts(monomer_amounts)
    reaction.set_conversion(conversion)
    reaction.set_chain_transfer_probability(chain_transfer_probability)
    return reaction

def dispersity_is_one():
    """
    DESCRIPTION: Generates a reaction with zeros in the entire reactivity ratio list
    RETURNS: Reaction object 
    """
    n = np.random.randint(2, 4)
    model = "Mayo-Lewis"
    max_monomer_amounts = 5000
    avg_DP = 1
    reactivity_ratios = np.random.rand(n, n-1)
    conversion = np.random.rand()
    monomer_amounts = np.random.randint(low=100, high=max_monomer_amounts, size=n)
    chain_transfer_probability = np.random.rand()
    reaction = Reaction(n, model)
    reaction.set_monomer_amounts(monomer_amounts)
    reaction.set_reactivity_ratios(reactivity_ratios)
    reaction.set_average_DP(avg_DP)
    reaction.set_monomer_amounts(monomer_amounts)
    reaction.set_conversion(conversion)
    reaction.set_chain_transfer_probability(chain_transfer_probability)
    return reaction

edge_case_reaction_generators = [
    some_zero_valued_reactivities,
    some_inf_valued_reactivities,
    all_zero_valued_reactivities,
    all_inf_valued_reactivities,
    dispersity_is_one,
    ]