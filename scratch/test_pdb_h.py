import sys
import os
import io

# Add the local directory to path to find atomipy
sys.path.append(os.getcwd())
import atomipy as ap

def test_h_writing():
    atoms = [
        {'index': 1, 'type': 'O', 'x': 0.0, 'y': 0.0, 'z': 0.0, 'resname': 'SOL'},
        {'index': 2, 'type': 'H', 'x': 0.0, 'y': 1.0, 'z': 0.0, 'resname': 'SOL'},
        {'index': 3, 'type': 'HW1', 'x': 1.0, 'y': 0.0, 'z': 0.0, 'resname': 'SOL'},
    ]
    box = [10, 10, 10]
    
    buf = io.StringIO()
    ap.write_pdb(atoms, box, buf)
    pdb_content = buf.getvalue()
    print("--- PDB CONTENT ---")
    print(pdb_content)
    print("-------------------")

if __name__ == "__main__":
    test_h_writing()
