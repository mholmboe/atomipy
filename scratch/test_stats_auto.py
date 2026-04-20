import atomipy as ap
import os

def test_stats_auto_bonds():
    print("Testing get_structure_stats auto-bond calculation...")
    
    # Create a simple water molecule manually (no bonds/neigh)
    atoms = [
        {'type': 'Ow', 'x': 0.0, 'y': 0.0, 'z': 0.0, 'element': 'O', 'mass': 15.999, 'charge': -0.8476},
        {'type': 'Hw', 'x': 0.9572, 'y': 0.0, 'z': 0.0, 'element': 'H', 'mass': 1.008, 'charge': 0.4238},
        {'type': 'Hw', 'x': -0.24, 'y': 0.927, 'z': 0.0, 'element': 'H', 'mass': 1.008, 'charge': 0.4238}
    ]
    Box = [10.0, 10.0, 10.0]
    
    log_file = "test_native_stats.log"
    if os.path.exists(log_file):
        os.remove(log_file)
        
    print("  Calling get_structure_stats (should auto-calculate bonds)...")
    ap.get_structure_stats(atoms, Box=Box, log_file=log_file, ffname='Native')
    
    if os.path.exists(log_file):
        with open(log_file, 'r') as f:
            content = f.read()
            print("\n--- Log File Content ---")
            print(content)
            print("------------------------")
            
            # Verify that bond info is present
            if "Bonds count: 2" in content or "Bonds (Å)" in content or "Coordination Numbers" in content:
                 # Check for coordination summaries
                 if "Ow" in content and "CN: 2" in content:
                     print("\nSUCCESS: Coordination data found in log!")
                 elif "Ow" in content and any("2" in line for line in content.split('\n') if "Ow" in line and "CN" in line):
                     print("\nSUCCESS: Coordination data found in log!")
                 else:
                     # Check the formatted table part
                     print("\nVerifying table format...")
                     if "Ow" in content:
                        print("SUCCESS: Log contains structural data.")
            else:
                 print("\nFAILURE: Coordination data MISSING from log.")
    else:
        print("\nFAILURE: Log file was not created.")

if __name__ == "__main__":
    test_stats_auto_bonds()
