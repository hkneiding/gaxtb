import pygad
import matplotlib.pyplot as plt
import matplotlib.image as img

import numpy as np

from gaxtb.molecule_handler import MoleculeHandler
from gaxtb.ga_functions import GaFunctions
from gaxtb.xtb_output_parser import XtbOutputParser
from gaxtb.xtb_runner import XtbRunner


def main():

    # problem definition / definition of chemical space to explore
    smiles_base_string = 'c1(/N=N/c2c(R)c(R)c(R)c(R)c2(R))c(R)c(R)c(R)c(R)c1(R)'
    smiles_gene_list = ['', 'C', 'O']
    # smiles_gene_list = ['', 'C', 'Cl', 'F', 'CO', 'C(F)(F)F', 'O']
    # definition of targets and weighing
    targets = ['dipole_moment']
    weights = [1]

    # start point
    # initial_population = [[0, 0, 0, 0, 0, 0], [1, 1, 1, 1, 1, 1], [2, 2, 2, 2, 2, 2], [3, 3, 3, 3, 3, 3], [4, 4, 4, 4, 4, 4], [5, 5, 5, 5, 5, 5], [6, 6, 6, 6, 6, 6], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]

    # hyper parameters
    population_size = 8
    num_generations = 10
    num_parents_mating = 2

    # fixed parameters (based on problem definition)
    n_genes = smiles_base_string.count('R')
    init_range_low = 0
    init_range_high = len(smiles_gene_list)

    # get functions for GA
    fitness_function = GaFunctions.get_fitness_function(targets, smiles_base_string, smiles_gene_list, weights=weights)
    on_start, on_stop, on_generation = GaFunctions.get_callback_functions_vis(smiles_base_string, smiles_gene_list, population_size)

    ga_instance = pygad.GA(num_generations=num_generations,
                            sol_per_pop=population_size,
                            num_parents_mating=num_parents_mating,
                            num_genes=n_genes,
                            fitness_func=fitness_function,
                            parent_selection_type="sss",
                            mutation_percent_genes=50,
                            crossover_type="single_point",
                            crossover_probability=0.5,
                            gene_type=int,
                            gene_space=[i for i in range(len(smiles_gene_list))],
                            init_range_low=init_range_low,
                            init_range_high=init_range_high,
                            on_generation=on_generation,
                            on_start=on_start,
                            on_stop=on_stop
                    )

    ga_instance.run()


if __name__ == "__main__":
    main()
