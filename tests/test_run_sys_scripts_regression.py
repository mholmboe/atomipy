import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
STRUCTURE_PDB = REPO_ROOT / "Kaolinite_GII_0.0487.pdb"


def gro_atom_count(gro_path: Path) -> int:
    with gro_path.open("r", encoding="utf-8") as handle:
        handle.readline()
        return int(handle.readline().strip())


def psf_atom_count(psf_path: Path) -> int:
    with psf_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if "!NATOM" in line:
                return int(line.split()[0])
    raise ValueError(f"Could not find !NATOM in {psf_path}")


class TestRunSysScriptsRegression(unittest.TestCase):
    def run_script_in_tempdir(self, script_name: str, checker):
        script_path = REPO_ROOT / "scripts" / script_name
        self.assertTrue(script_path.exists(), f"Missing script: {script_path}")

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            env = os.environ.copy()
            pythonpath = env.get("PYTHONPATH")
            env["PYTHONPATH"] = (
                f"{REPO_ROOT}{os.pathsep}{pythonpath}" if pythonpath else str(REPO_ROOT)
            )

            cmd = [sys.executable, str(script_path), str(STRUCTURE_PDB)]
            proc = subprocess.run(
                cmd,
                cwd=tmp_path,
                env=env,
                capture_output=True,
                text=True,
                timeout=180,
            )
            if proc.returncode != 0:
                self.fail(
                    f"{script_name} failed with code {proc.returncode}\n"
                    f"stdout:\n{proc.stdout}\n\nstderr:\n{proc.stderr}"
                )
            checker(tmp_path)

    def test_run_sys2minff_psf_matches_full_system_gro(self):
        def checker(tmp_path):
            preem_gro = tmp_path / "preem.gro"
            minff_psf = tmp_path / "minff.psf"
            self.assertTrue(preem_gro.exists())
            self.assertTrue(minff_psf.exists())
            self.assertTrue((tmp_path / "minff.itp").exists())
            self.assertTrue((tmp_path / "minff.data").exists())
            self.assertTrue((tmp_path / "minff_no150angles.data").exists())
            self.assertEqual(psf_atom_count(minff_psf), gro_atom_count(preem_gro))

        self.run_script_in_tempdir("run_sys2minff.py", checker)

    def test_run_sys2clayff_psf_matches_full_system_gro(self):
        def checker(tmp_path):
            preem_gro = tmp_path / "preem.gro"
            clayff_psf = tmp_path / "clayff.psf"
            self.assertTrue(preem_gro.exists())
            self.assertTrue(clayff_psf.exists())
            self.assertTrue((tmp_path / "clayff.itp").exists())
            self.assertTrue((tmp_path / "clayff.data").exists())
            self.assertEqual(psf_atom_count(clayff_psf), gro_atom_count(preem_gro))

        self.run_script_in_tempdir("run_sys2clayff.py", checker)


if __name__ == "__main__":
    unittest.main()
