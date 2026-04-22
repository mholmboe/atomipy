import atomipy as ap
import numpy as np
import sys

def main():
    try:
        # Load Kaolinite
        atoms, Box = ap.import_gro("Kaolinite_GII_0.0487.gro")
        print(f"Loaded {len(atoms)} atoms")

        # Find an Oxygen that is bonded to a Hydrogen
        h_idx = -1
        o_idx = -1
        for i, a in enumerate(atoms):
            if a.get('element') == 'H':
                h_pos = np.array([a['x'], a['y'], a['z']])
                for j, o in enumerate(atoms):
                    if o.get('element') == 'O':
                        o_pos = np.array([o['x'], o['y'], o['z']])
                        if np.linalg.norm(o_pos - h_pos) < 1.2:
                            h_idx = i
                            o_idx = j
                            break
                if h_idx != -1: break
        
        if h_idx == -1:
            print("Could not find an O-H pair to test with.")
            sys.exit(1)
            
        target_o_type = atoms[o_idx]['type']
        print(f"Removing Hydrogen (index {h_idx}) bonded to Oxygen (index {o_idx}, type {target_o_type})")
        
        # Create a copy without that H
        atoms_test = [a for i, a in enumerate(atoms) if i != h_idx]
        
        # Test add_hydrogens_bvs
        print("Running BVS analysis to check deficit...")
        report = ap.analyze_bvs(atoms_test, Box)
        # Find the oxygen again in the new list
        # It's type should be unique enough or we find by position
        target_pos = np.array([atoms[o_idx]['x'], atoms[o_idx]['y'], atoms[o_idx]['z']])
        target_row = None
        for r in report['results']:
            idx = r['index'] - 1
            a = atoms_test[idx]
            a_pos = np.array([a['x'], a['y'], a['z']])
            if np.linalg.norm(a_pos - target_pos) < 0.1:
                target_row = r
                break
        
        if target_row:
            print(f"Target Oxygen BVS: {target_row['bvs']:.3f}, Delta: {target_row['delta']:.3f}, CN: {target_row['cn']}")
        else:
            print("Could not find target Oxygen in BVS report.")

        print("Running add_hydrogens_bvs with coordination=3...")
        new_atoms = ap.add_hydrogens_bvs(atoms_test, Box, delta_threshold=-0.2, max_additions=1, coordination=3)
        
        print(f"Original test count: {len(atoms_test)}")
        print(f"New atom count: {len(new_atoms)}")
        
        if len(new_atoms) > len(atoms_test):
            print("✓ Successfully added a hydrogen atom based on BVS deficit!")
            # Check position
            new_h = new_atoms[-1]
            new_h_pos = np.array([new_h['x'], new_h['y'], new_h['z']])
            print(f"New H position: {new_h_pos}, Dist to O: {np.linalg.norm(new_h_pos - target_pos):.3f} Å")
        else:
            print("✗ Failed to add a hydrogen atom. Likely coordination check issue in add_H_atom.")

    except Exception as e:
        print(f"Error during verification: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
