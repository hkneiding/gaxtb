import itertools
from openbabel import openbabel
from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt


class MoleculeHandler:

    @staticmethod
    def convert_smiles_to_xyz(smiles: str):
        return MoleculeHandler.convert_smiles(smiles, 'xyz')

    @staticmethod
    def convert_smiles(smiles: str, format: str):

        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("smi", format)

        molecule = openbabel.OBMol()
        obConversion.ReadString(molecule, smiles)
        molecule.AddHydrogens()

        gen3d = openbabel.OBOp.FindType("Gen3D")
        gen3d.Do(molecule, '--best')

        return obConversion.WriteString(molecule)

    @staticmethod
    def get_full_smiles_string(smiles_base_string: str, smiles_genes: list[str]):

            smiles_split = smiles_base_string.split('R')

            # check for correct dimensions
            assert len(smiles_split) == len(smiles_genes) + 1

            # get full smiles string
            full_smiles_string = ''.join([x for x in itertools.chain.from_iterable(itertools.zip_longest(smiles_split, smiles_genes)) if x])

            return full_smiles_string.replace('()', '')

    @staticmethod
    def get_smiles_genes(chromosome: list[int] , smiles_gene_list: list[str]):
        return [smiles_gene_list[x] for x in chromosome]

    @staticmethod
    def save_image_from_smiles(smiles_string: str, file_path: str):

        molecule = Chem.MolFromSmiles(smiles_string)
        figure = Draw.MolToMPL(molecule, size=(120,120))
        plt.axis('off')
        figure.savefig(file_path, bbox_inches='tight')
        figure.clear()
        plt.close(figure)
