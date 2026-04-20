import atomipy as ap
import os

def test_stats_auto_bonds_and_charges():
    print("Testing get_structure_stats auto-bond and auto-charge calculation...")
    
    # Create a simple water molecule manually (no bonds/neigh/charges)
    atoms = [
        {'type': 'Ow', 'x': 0.0, 'y': 0.0, 'z': 0.0, 'element': 'O', 'resname': 'SOL'},
        {'type': 'Hw', 'x': 0.9572, 'y': 0.0, 'z': 0.0, 'element': 'H', 'resname': 'SOL'},
        {'type': 'Hw', 'x': -0.24, 'y': 0.927, 'z': 0.0, 'element': 'H', 'resname': 'SOL'}
    ]
    Box = [10.0, 10.0, 10.0]
    
    log_file = "test_native_stats_charges.log"
    if os.path.exists(log_file):
        os.remove(log_file)
        
    print("  Calling get_structure_stats (should auto-calculate bonds and charges)...")
    ap.get_structure_stats(atoms, Box=Box, log_file=log_file, ffname='Native')
    
    if os.path.exists(log_file):
        with open(log_file, 'r') as f:
            content = f.read()
            print("\n--- Log File Content ---")
            print(content)
            print("------------------------")
            
            # Verify that bond info is present
            if "Bonds count: 2" in content or "Bonds (Å)" in content or "Coordination Numbers" in content or "Angle Statistics" in content:
                 print("\nSUCCESS: Coordination/Bond data found in log!")
                 
                 # Check for charges in content (Ow charge is -0.89517 in OPC3)
                 if "-0.89517" in content or "0.44759" in content:
                     print("SUCCESS: Auto-assigned charges found in log!")
                 else:
                     print("FAILURE: Auto-assigned charges MISSING from log.")
            else:
                 print("\nFAILURE: Structural data MISSING from log.")
    else:
        print("\nFAILURE: Log file was not created.")

if __name__ == "__main__":
    test_stats_auto_bonds_and_charges()
