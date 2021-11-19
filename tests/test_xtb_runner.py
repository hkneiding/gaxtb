import unittest
from parameterized import parameterized
from .utils.helpers import assertDeepAlmostEqual

from gaxtb.xtb_runner import XtbRunner
from gaxtb.xtb_output_parser import XtbOutputParser
from gaxtb.molecule_handler import MoleculeHandler


class TestXtbRunner(unittest.TestCase):

    @parameterized.expand([

        [
            'C=C',
            { 'homo_energy': -11.3190, 'lumo_energy': -5.6407, 'dipole_moment': 0.000, 'heat_capacity': 10.2312, 'entropy': 52.3017,
            'zpve': 0.049895220647, 'energy': -6.271260107426, 'enthalpy': -6.217359643534, 'free_energy':  -6.242209880088, 'homo_lumo_gap': 5.678272522007}
        ]

    ])
    def test_xtb_run_full(self, smiles_string, expected):

        xtbr = XtbRunner()

        molecule_data = MoleculeHandler.convert_smiles_to_xyz(smiles_string)
        xtb_output = xtbr.run_xtb_full(molecule_data, 'xyz')

        result = XtbOutputParser(xtb_output).parse()

        assertDeepAlmostEqual(self, result, expected, places=3)
