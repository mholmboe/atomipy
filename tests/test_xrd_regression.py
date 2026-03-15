import warnings
import unittest
from pathlib import Path

import atomipy as ap
import atomipy.diffraction as diffraction


REPO_ROOT = Path(__file__).resolve().parents[1]
STRUCTURE_PDB = REPO_ROOT / "Kaolinite_GII_0.0487.pdb"


class TestXrdRegression(unittest.TestCase):
    def setUp(self):
        atoms, cell = ap.import_pdb(str(STRUCTURE_PDB))
        self.atoms = ap.element(atoms)
        self.cell = cell

    def test_xrd_no_invalid_arcsin_runtime_warning(self):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            two_theta, intensity = ap.xrd(
                atoms=self.atoms,
                Box=self.cell,
                two_theta_range=(2.0, 90.0),
                angle_step=0.05,
                plot=False,
                save_output=False,
            )

        runtime_messages = [
            str(w.message)
            for w in caught
            if issubclass(w.category, RuntimeWarning)
        ]
        self.assertFalse(
            any("invalid value encountered in arcsin" in msg for msg in runtime_messages),
            runtime_messages,
        )
        self.assertEqual(len(two_theta), len(intensity))
        self.assertGreater(len(two_theta), 0)

    def test_xrd_plot_false_works_without_optional_plotting_dependencies(self):
        old_plt = diffraction.plt
        old_find_peaks = diffraction.find_peaks
        diffraction.plt = None
        diffraction.find_peaks = None
        try:
            two_theta, intensity = ap.xrd(
                atoms=self.atoms,
                Box=self.cell,
                two_theta_range=(5.0, 20.0),
                angle_step=0.1,
                plot=False,
                save_output=False,
            )
        finally:
            diffraction.plt = old_plt
            diffraction.find_peaks = old_find_peaks

        self.assertEqual(len(two_theta), len(intensity))
        self.assertGreater(len(two_theta), 0)

    def test_xrd_plot_true_requires_matplotlib(self):
        old_plt = diffraction.plt
        diffraction.plt = None
        try:
            with self.assertRaises(ImportError):
                ap.xrd(
                    atoms=self.atoms,
                    Box=self.cell,
                    two_theta_range=(5.0, 20.0),
                    angle_step=0.1,
                    plot=True,
                    save_output=False,
                )
        finally:
            diffraction.plt = old_plt


if __name__ == "__main__":
    unittest.main()
