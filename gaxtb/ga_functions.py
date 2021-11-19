import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as img


from .molecule_handler import MoleculeHandler
from .xtb_output_parser import XtbOutputParser
from .xtb_runner import XtbRunner

class GaFunctions:

    @staticmethod
    def get_fitness_function(xtb_targets, smiles_base_string, smiles_gene_list, weights=None):

        xtbr = XtbRunner()

        xtb_targets = xtb_targets
        smiles_base_string = smiles_base_string
        smiles_gene_list = smiles_gene_list

        if weights == None:
            weights = [1 for x in xtb_targets]
        else:
            # check that the length of xtb_targets and weights are equal
            assert len(xtb_targets) == len(weights)
            weights = weights

        def fitness_function(solution, solution_id):

            smiles_genes = MoleculeHandler.get_smiles_genes(solution, smiles_gene_list)
            full_smiles_string = MoleculeHandler.get_full_smiles_string(smiles_base_string, smiles_genes)

            molecule_file_content = MoleculeHandler.convert_smiles_to_xyz(full_smiles_string)

            # run xtb
            xtb_output = xtbr.run_xtb_full(molecule_file_content, 'xyz')
            xtb_output_data = XtbOutputParser(xtb_output).parse()

            fitness = 0
            for i in range(len(xtb_targets)):
                fitness += weights[i] * xtb_output_data[xtb_targets[i]]

            return fitness

        return fitness_function

    @staticmethod
    def get_callback_functions_vis(smiles_base_string, smiles_gene_list, population_size):

        # set up plotting window
        fig, ax = plt.subplots(nrows=int(np.ceil(population_size / 3)), ncols=3)
        plt.ion()
        plt.show()
        fig.canvas.manager.full_screen_toggle()
        plt.draw()

        def on_start(ga_instance):

            plt.suptitle('Generation ' + str(ga_instance.generations_completed))
            plt.draw()
            plt.pause(0.01)

        def on_stop(ga_instance, fitness_values):

            input('Press Enter to exit')

        def on_generation(ga_instance):

            # get image data
            images = []
            for i in range(len(ga_instance.population)):

                # get full smiles string
                smiles_genes = MoleculeHandler.get_smiles_genes(ga_instance.population[i], smiles_gene_list)
                full_smiles_string = MoleculeHandler.get_full_smiles_string(smiles_base_string, smiles_genes)

                # get image data
                MoleculeHandler.save_image_from_smiles(full_smiles_string, 'mol.jpeg')
                images.append(img.imread('mol.jpeg'))

            # get descending order of fitnesses
            order_list = np.argsort(-ga_instance.last_generation_fitness)
            
            # plot individual subplots
            i = 0
            for row in ax:
                for col in row:

                    # skip if all individuals have been plotted already
                    if i < len(order_list):
                        col.imshow(images[i])
                        col.title.set_text('Fitness: ' + str(np.round(ga_instance.last_generation_fitness[order_list[i]], decimals=3)))
                    col.axis('off')
                    i += 1

            # add title
            plt.suptitle('Generation ' + str(ga_instance.generations_completed))
            # update plot
            plt.draw()
            plt.pause(0.01)

        return on_start, on_stop, on_generation