#!/usr/bin/env python
"""
Test script to compare performance of original and optimized cell_list_dist_matrix.
"""

import numpy as np
import sys
import os
import time

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# Import the module to test
from atomipy.cell_list_dist_matrix import cell_list_dist_matrix
from atomipy.import_module import gro

def generate_test_system(N=1000, box_size=100.0):
    """Generate a test system with random atom positions."""
    atoms = []
    for i in range(N):
        atom = {
            'x': np.random.uniform(0, box_size),
            'y': np.random.uniform(0, box_size),
            'z': np.random.uniform(0, box_size),
            'type': 'C' if np.random.random() > 0.2 else 'H'
        }
        atoms.append(atom)
    
    # Create box dimensions
    box_dim = [box_size, box_size, box_size, 0, 0, 0, 0, 0, 0]
    
    return atoms, box_dim

def load_real_system():
    """Load a real molecular system from the GRO file."""
    try:
        atoms = gro("./preem.gro")
        # Extract box dimensions from first atom
        if 'box_dim' in atoms[0]:
            box_dim = [float(atoms[0]['box_dim'][0]), 
                      float(atoms[0]['box_dim'][1]), 
                      float(atoms[0]['box_dim'][2]),
                      0, 0, 0, 0, 0, 0]  # Assuming orthogonal box
            return atoms, box_dim
        else:
            print("No box dimensions found in the imported file")
            return None, None
    except Exception as e:
        print(f"Error loading real system: {e}")
        return None, None

def benchmark_cell_list(atoms, box_dim, num_runs=3, cutoff=2.45):
    """Benchmark the cell_list_dist_matrix function."""
    times = []
    
    for i in range(num_runs):
        start_time = time.time()
        dist_matrix, X_dist, Y_dist, Z_dist, bond_list, dist_list = cell_list_dist_matrix(
            atoms, cutoff=cutoff, Box_dim=box_dim)
        end_time = time.time()
        times.append(end_time - start_time)
        
        print(f"Run {i+1}: {times[-1]:.4f} seconds, Found {len(bond_list)} bonds")
    
    return np.mean(times), dist_matrix, X_dist, Y_dist, Z_dist, bond_list, dist_list

def test_conversion_functions(dist_matrix, X_dist, Y_dist, Z_dist, cutoff=2.45):
    """Test the optimized conversion functions."""
    from atomipy.cell_list_dist_matrix import convert_to_sparse_dict, get_neighbors
    
    print("\nTesting optimized conversion functions:")
    
    # Test convert_to_sparse_dict
    start_time = time.time()
    distance_dict = convert_to_sparse_dict(dist_matrix, X_dist, Y_dist, Z_dist, cutoff)
    end_time = time.time()
    
    print(f"  convert_to_sparse_dict: {end_time - start_time:.4f} seconds, {len(distance_dict)} entries")
    
    # Test get_neighbors for a few random atoms
    N = dist_matrix.shape[0]
    sample_indices = np.random.choice(N, min(5, N), replace=False)
    
    for atom_idx in sample_indices:
        start_time = time.time()
        neighbors = get_neighbors(dist_matrix, X_dist, Y_dist, Z_dist, atom_idx, cutoff)
        end_time = time.time()
        
        print(f"  get_neighbors for atom {atom_idx}: {end_time - start_time:.4f} seconds, {len(neighbors)} neighbors")

def main():
    print("Testing cell_list_dist_matrix optimizations")
    print("-" * 60)
    
    # Try to load real system first
    atoms, box_dim = load_real_system()
    
    if atoms is None:
        print("Using generated test system...")
        atoms, box_dim = generate_test_system(N=2000, box_size=50.0)
        test_type = "generated"
    else:
        print(f"Using real system with {len(atoms)} atoms")
        test_type = "real"
    
    print(f"\nBenchmarking cell_list_dist_matrix ({test_type} system):")
    mean_time, dist_matrix, X_dist, Y_dist, Z_dist, bond_list, dist_list = benchmark_cell_list(
        atoms, box_dim)
    
    print(f"\nAverage time: {mean_time:.4f} seconds")
    
    # Test the optimized conversion functions
    test_conversion_functions(dist_matrix, X_dist, Y_dist, Z_dist)
    
    print("\nTest completed!")

if __name__ == "__main__":
    main()
