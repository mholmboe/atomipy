#!/usr/bin/env python3
"""
Test script to compare the Bond_index and Angle_index output variables and performance
of the bond_angle function when using dist_matrix vs. cell_list_dist_matrix.
"""

import os
import time
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
import sys

# Add the parent directory to sys.path to ensure imports work correctly
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# Import necessary modules
import atomipy as ap
from atomipy.dist_matrix import dist_matrix
from atomipy.cell_list_dist_matrix import cell_list_dist_matrix

# Create versions of bond_angle with each distance matrix implementation
def bond_angle_with_dist_matrix(atoms, Box_dim, rmaxH=1.2, rmaxM=2.45, same_element_bonds=False, same_molecule_only=True):
    """Implementation using dist_matrix"""
    # Get the number of atoms
    N = len(atoms)

    # Get all atom pairs with periodic boundary correction using dist_matrix
    dmat, dx, dy, dz = dist_matrix(atoms, Box_dim=Box_dim)
    
    # Create bond lists based on the distance matrix and appropriate cutoffs
    precalc_bond_list = []
    dist_list = []
    
    # Create bond lists based on the distance matrix and appropriate cutoffs
    for i in range(N):
        for j in range(i+1, N):  # Only consider each pair once
            if dmat[i, j] > 0:  # Skip diagonal and zero distances
                # Determine if either atom is hydrogen
                el_i = atoms[i].get('element', 'X')
                el_j = atoms[j].get('element', 'X')
                isH_i = el_i == 'H'
                isH_j = el_j == 'H'
                
                # Apply appropriate cutoff based on atom types
                if isH_i or isH_j:
                    cutoff = rmaxH  # Use hydrogen cutoff
                else:
                    cutoff = rmaxM  # Use non-hydrogen cutoff
                    
                # If within cutoff, add to bond list
                if dmat[i, j] <= cutoff:
                    precalc_bond_list.append([i, j])
                    dist_list.append(dmat[i, j])

    # Initialize lists for all atoms
    for i in range(N):
        atoms[i]['neigh'] = []
        atoms[i]['bonds'] = []
        atoms[i]['angles'] = []
    
    # Filter bonds based on element types and molecule IDs
    bond_pairs = []  # Store bonds as (atom1_idx, atom2_idx, distance)
    
    # At this point, precalc_bond_list and dist_list should be populated
    # from the loop that checked distance matrix against cutoffs
    if len(precalc_bond_list) > 0:
        for k in range(len(precalc_bond_list)):
            i, j = precalc_bond_list[k]
            distance = dist_list[k]
            
            # Ensure i < j for consistency - smaller index always in first column
            if i > j:
                i, j = j, i
                
            el_i = atoms[i].get('element','X')
            el_j = atoms[j].get('element','X')
            
            # Get molecule IDs if available, otherwise use None
            molid_i = atoms[i].get('molid', None)
            molid_j = atoms[j].get('molid', None)
            
            # Apply element check and molecule check if needed
            molecule_condition = True if not same_molecule_only else (molid_i == molid_j)
            element_condition = same_element_bonds or el_i != el_j
            
            if element_condition and molecule_condition:
                # Add to both atoms' neighbor and bond lists
                atoms[i]['neigh'].append(j)
                atoms[i]['bonds'].append((j, distance))
                
                atoms[j]['neigh'].append(i)
                atoms[j]['bonds'].append((i, distance))
                
                # Store bond information as tuple (low_idx, high_idx, distance)
                bond_pairs.append((i, j, distance))
    
    # Calculate angles for atoms with bonds
    angle_data = []  # Store angle data
    
    for i in range(N):
        # Skip if atom has less than 2 bonds
        if len(atoms[i]['neigh']) < 2:
            continue
            
        # Compute angles for each pair of bonded neighbors
        for m in range(len(atoms[i]['neigh'])):
            for n in range(m+1, len(atoms[i]['neigh'])):
                j = atoms[i]['neigh'][m]
                k = atoms[i]['neigh'][n]
                
                # Get vectors from atom i to atoms j and k with PBC correction
                rij = np.array([dx[i, j], dy[i, j], dz[i, j]])
                rik = np.array([dx[i, k], dy[i, k], dz[i, k]])
                
                # Normalize vectors
                rij_norm = np.linalg.norm(rij)
                rik_norm = np.linalg.norm(rik)
                
                # Calculate angle using dot product
                cos_angle = np.dot(rij, rik) / (rij_norm * rik_norm)
                
                # Clamp to valid range to prevent numerical errors
                cos_angle = max(min(cos_angle, 1.0), -1.0)
                angle = np.degrees(np.arccos(cos_angle))
                
                # Add angle to atom's data
                atoms[i]['angles'].append(((j, k), angle))
                
                # Store angle data with proper ordering for Angle_index
                # Ensure first atom has lower index than third atom
                if j < k:
                    atom1, atom3 = j, k
                    # Vector from middle atom (i) to lowest index atom (j)
                    dx12, dy12, dz12 = dx[i, j], dy[i, j], dz[i, j]
                    # Vector from middle atom (i) to highest index atom (k)
                    dx23, dy23, dz23 = dx[i, k], dy[i, k], dz[i, k]
                else:
                    atom1, atom3 = k, j
                    # Vector from middle atom (i) to lowest index atom (k)
                    dx12, dy12, dz12 = dx[i, k], dy[i, k], dz[i, k]
                    # Vector from middle atom (i) to highest index atom (j)
                    dx23, dy23, dz23 = dx[i, j], dy[i, j], dz[i, j]
                
                # Store the angle data in consistent format
                angle_data.append((atom1, i, atom3, angle, 
                                 dx12, dy12, dz12, 
                                 dx23, dy23, dz23))
    
    # Convert bond_pairs list to Nx3 numpy array
    Bond_index = np.array(bond_pairs)
    
    # For each bond, ensure the smaller atom index is in the first column
    if len(Bond_index) > 0:
        # This should already be taken care of when creating bond_pairs,
        # but let's make sure by doing a final check
        for i in range(len(Bond_index)):
            if Bond_index[i, 0] > Bond_index[i, 1]:
                # Swap indices to put smaller first
                Bond_index[i, 0], Bond_index[i, 1] = Bond_index[i, 1], Bond_index[i, 0]
        
        # Now sort rows based on first column (atom1_idx) and then second column (atom2_idx)
        sorted_indices = np.lexsort((Bond_index[:, 1], Bond_index[:, 0]))
        Bond_index = Bond_index[sorted_indices]
    
    # Convert angle_data list to Mx10 numpy array
    Angle_index = np.array(angle_data)
    
    # Sort Angle_index row-wise by atomic indices
    if len(Angle_index) > 0:
        # Sort the entire array: first by middle atom (column 1), 
        # then by lowest bonded atom (column 0), then by highest bonded atom (column 2)
        sorted_indices = np.lexsort((Angle_index[:, 2], Angle_index[:, 0], Angle_index[:, 1]))
        Angle_index = Angle_index[sorted_indices]
    
    return atoms, Bond_index, Angle_index

def bond_angle_with_cell_list(atoms, Box_dim, rmaxH=1.2, rmaxM=2.45, same_element_bonds=False, same_molecule_only=True):
    """Implementation using cell_list_dist_matrix"""
    # Get the number of atoms
    N = len(atoms)

    # Get all atom pairs within the bond cutoff using cell_list approach
    # This calculates distances applying PBC and different cutoffs for H and other atoms
    dmat, dx, dy, dz, precalc_bond_list, dist_list = cell_list_dist_matrix(
        atoms, cutoff=rmaxM, Box_dim=Box_dim, rmaxH=rmaxH, H_type='H'
    )

    # Initialize lists for all atoms
    for i in range(N):
        atoms[i]['neigh'] = []
        atoms[i]['bonds'] = []
        atoms[i]['angles'] = []
    
    # Filter bonds based on element types and molecule IDs
    bond_pairs = []  # Store bonds as (atom1_idx, atom2_idx, distance)
    
    # Process the precalculated bonds from cell_list_dist_matrix
    if len(precalc_bond_list) > 0:
        for k in range(len(precalc_bond_list)):
            i, j = precalc_bond_list[k]
            distance = dist_list[k]
            
            # Ensure i < j for consistency - smaller index always in first column
            if i > j:
                i, j = j, i
                
            el_i = atoms[i].get('element','X')
            el_j = atoms[j].get('element','X')
            
            # Get molecule IDs if available, otherwise use None
            molid_i = atoms[i].get('molid', None)
            molid_j = atoms[j].get('molid', None)
            
            # Apply element check and molecule check if needed
            molecule_condition = True if not same_molecule_only else (molid_i == molid_j)
            element_condition = same_element_bonds or el_i != el_j
            
            if element_condition and molecule_condition:
                # Add to both atoms' neighbor and bond lists
                atoms[i]['neigh'].append(j)
                atoms[i]['bonds'].append((j, distance))
                
                atoms[j]['neigh'].append(i)
                atoms[j]['bonds'].append((i, distance))
                
                # Store bond information as tuple (low_idx, high_idx, distance)
                bond_pairs.append((i, j, distance))
    
    # Calculate angles for atoms with bonds
    angle_data = []  # Store angle data
    
    for i in range(N):
        # Skip if atom has less than 2 bonds
        if len(atoms[i]['neigh']) < 2:
            continue
            
        # Compute angles for each pair of bonded neighbors
        for m in range(len(atoms[i]['neigh'])):
            for n in range(m+1, len(atoms[i]['neigh'])):
                j = atoms[i]['neigh'][m]
                k = atoms[i]['neigh'][n]
                
                # Get vectors from atom i to atoms j and k with PBC correction
                rij = np.array([dx[i, j], dy[i, j], dz[i, j]])
                rik = np.array([dx[i, k], dy[i, k], dz[i, k]])
                
                # Normalize vectors
                rij_norm = np.linalg.norm(rij)
                rik_norm = np.linalg.norm(rik)
                
                # Calculate angle using dot product
                cos_angle = np.dot(rij, rik) / (rij_norm * rik_norm)
                
                # Clamp to valid range to prevent numerical errors
                cos_angle = max(min(cos_angle, 1.0), -1.0)
                angle = np.degrees(np.arccos(cos_angle))
                
                # Add angle to atom's data
                atoms[i]['angles'].append(((j, k), angle))
                
                # Store angle data with proper ordering for Angle_index
                # Ensure first atom has lower index than third atom
                if j < k:
                    atom1, atom3 = j, k
                    # Vector from middle atom (i) to lowest index atom (j)
                    dx12, dy12, dz12 = dx[i, j], dy[i, j], dz[i, j]
                    # Vector from middle atom (i) to highest index atom (k)
                    dx23, dy23, dz23 = dx[i, k], dy[i, k], dz[i, k]
                else:
                    atom1, atom3 = k, j
                    # Vector from middle atom (i) to lowest index atom (k)
                    dx12, dy12, dz12 = dx[i, k], dy[i, k], dz[i, k]
                    # Vector from middle atom (i) to highest index atom (j)
                    dx23, dy23, dz23 = dx[i, j], dy[i, j], dz[i, j]
                
                # Store the angle data in consistent format
                angle_data.append((atom1, i, atom3, angle, 
                                 dx12, dy12, dz12, 
                                 dx23, dy23, dz23))
    
    # Convert bond_pairs list to Nx3 numpy array
    Bond_index = np.array(bond_pairs)
    
    # For each bond, ensure the smaller atom index is in the first column
    if len(Bond_index) > 0:
        # This should already be taken care of when creating bond_pairs,
        # but let's make sure by doing a final check
        for i in range(len(Bond_index)):
            if Bond_index[i, 0] > Bond_index[i, 1]:
                # Swap indices to put smaller first
                Bond_index[i, 0], Bond_index[i, 1] = Bond_index[i, 1], Bond_index[i, 0]
        
        # Now sort rows based on first column (atom1_idx) and then second column (atom2_idx)
        sorted_indices = np.lexsort((Bond_index[:, 1], Bond_index[:, 0]))
        Bond_index = Bond_index[sorted_indices]
    
    # Convert angle_data list to Mx10 numpy array
    Angle_index = np.array(angle_data)
    
    # Sort Angle_index row-wise by atomic indices
    if len(Angle_index) > 0:
        # Sort the entire array: first by middle atom (column 1), 
        # then by lowest bonded atom (column 0), then by highest bonded atom (column 2)
        sorted_indices = np.lexsort((Angle_index[:, 2], Angle_index[:, 0], Angle_index[:, 1]))
        Angle_index = Angle_index[sorted_indices]
    
    return atoms, Bond_index, Angle_index

def compare_arrays(arr1, arr2, name):
    """Compare two arrays and return statistics about their differences"""
    if len(arr1) != len(arr2):
        print(f"{name} lengths differ: {len(arr1)} vs {len(arr2)}")
        return False, None
    
    if len(arr1) == 0:
        print(f"Both {name} arrays are empty")
        return True, None
    
    # Check if shapes match
    if arr1.shape != arr2.shape:
        print(f"{name} shapes differ: {arr1.shape} vs {arr2.shape}")
        return False, None
    
    # Calculate absolute differences
    abs_diff = np.abs(arr1 - arr2)
    
    # Check if arrays are identical
    if np.all(abs_diff < 1e-10):  # Using small threshold for floating point comparison
        print(f"{name} arrays are identical")
        return True, None
    
    # Calculate statistics
    max_diff = np.max(abs_diff)
    mean_diff = np.mean(abs_diff)
    median_diff = np.median(abs_diff)
    std_diff = np.std(abs_diff)
    
    # Count number of significant differences (greater than small threshold)
    significant_diff_count = np.sum(abs_diff > 1e-6)
    
    return False, {
        'max_diff': max_diff,
        'mean_diff': mean_diff,
        'median_diff': median_diff,
        'std_diff': std_diff,
        'significant_diff_count': significant_diff_count,
        'total_elements': arr1.size
    }

def run_comparison(gro_file, runs=3, subset_size=None):
    """Run comparison between the two implementations and report results"""
    print(f"Loading structure from {gro_file}...")
    atoms, box_dim = ap.import_gro(gro_file)
    
    # Process atoms to ensure they have all required fields
    print("Preprocessing atoms for element assignment...")
    for i in range(len(atoms)):
        # Transfer atname to type field for better element detection
        if 'atname' in atoms[i] and not atoms[i].get('type', ''):
            atoms[i]['type'] = atoms[i]['atname']
    
    # Apply element name extraction using atomipy.element
    print("Assigning elements to atoms...")
    # ap.element expects atom dictionaries, not strings
    atoms = ap.element(atoms)
    
    # Check if molid exists in the atoms
    has_molid = any('molid' in atom for atom in atoms)
    print(f"Molecule IDs found in data: {'Yes' if has_molid else 'No'}")
    
    if subset_size and subset_size < len(atoms):
        print(f"Using subset of {subset_size} atoms for testing")
        atoms = atoms[:subset_size]
    
    print(f"Structure contains {len(atoms)} atoms")
    print("-" * 60)

    # Times for each implementation
    dist_matrix_times = []
    cell_list_times = []
    
    # Run multiple times to get average performance
    for run in range(1, runs+1):
        print(f"Run {run}/{runs}")
        
        # Make deep copies of atoms to avoid any side effects between runs
        atoms_copy1 = deepcopy(atoms)
        atoms_copy2 = deepcopy(atoms)
        
        # Run dist_matrix implementation
        print("\nRunning bond_angle with dist_matrix...")
        start_time = time.time()
        atoms_dm, bond_index_dm, angle_index_dm = bond_angle_with_dist_matrix(
            atoms_copy1, box_dim, rmaxH=1.2, rmaxM=2.45, same_element_bonds=False, same_molecule_only=True
        )
        end_time = time.time()
        dist_matrix_time = end_time - start_time
        dist_matrix_times.append(dist_matrix_time)
        print(f"Execution time: {dist_matrix_time:.4f} seconds")
        print(f"Found {len(bond_index_dm)} bonds and {len(angle_index_dm)} angles")
        
        # Run cell_list implementation
        print("\nRunning bond_angle with cell_list_dist_matrix...")
        start_time = time.time()
        atoms_cl, bond_index_cl, angle_index_cl = bond_angle_with_cell_list(
            atoms_copy2, box_dim, rmaxH=1.2, rmaxM=2.45, same_element_bonds=False, same_molecule_only=True
        )
        end_time = time.time()
        cell_list_time = end_time - start_time
        cell_list_times.append(cell_list_time)
        print(f"Execution time: {cell_list_time:.4f} seconds")
        print(f"Found {len(bond_index_cl)} bonds and {len(angle_index_cl)} angles")
        
        # Compare bond indices
        print("\nComparing Bond_index arrays:")
        identical_bonds, bond_stats = compare_arrays(bond_index_dm, bond_index_cl, "Bond_index")
        if not identical_bonds and bond_stats:
            print(f"  Max difference: {bond_stats['max_diff']}")
            print(f"  Mean difference: {bond_stats['mean_diff']}")
            print(f"  Significant differences: {bond_stats['significant_diff_count']} out of {bond_stats['total_elements']} elements")
        
        # Compare angle indices
        print("\nComparing Angle_index arrays:")
        identical_angles, angle_stats = compare_arrays(angle_index_dm, angle_index_cl, "Angle_index")
        if not identical_angles and angle_stats:
            print(f"  Max difference: {angle_stats['max_diff']}")
            print(f"  Mean difference: {angle_stats['mean_diff']}")
            print(f"  Significant differences: {angle_stats['significant_diff_count']} out of {angle_stats['total_elements']} elements")
        
        print("-" * 60)
    
    # Calculate average times
    avg_dist_matrix_time = sum(dist_matrix_times) / runs
    avg_cell_list_time = sum(cell_list_times) / runs
    
    # Print performance summary
    print("\nPerformance Summary:")
    print(f"Average dist_matrix time: {avg_dist_matrix_time:.4f} seconds")
    print(f"Average cell_list_dist_matrix time: {avg_cell_list_time:.4f} seconds")
    speedup = avg_dist_matrix_time / avg_cell_list_time if avg_cell_list_time > 0 else float('inf')
    if speedup > 1:
        print(f"cell_list_dist_matrix is {speedup:.2f}x faster")
    else:
        print(f"dist_matrix is {1/speedup:.2f}x faster")
    
    # Create performance comparison chart
    plt.figure(figsize=(10, 6))
    plt.bar(['dist_matrix', 'cell_list_dist_matrix'], 
            [avg_dist_matrix_time, avg_cell_list_time],
            color=['blue', 'orange'])
    plt.title('Performance Comparison: Bond Angle Calculation')
    plt.ylabel('Time (seconds)')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    for i, v in enumerate([avg_dist_matrix_time, avg_cell_list_time]):
        plt.text(i, v + 0.1, f"{v:.4f}s", ha='center')
    
    plt.savefig('bond_angle_performance_comparison.png')
    print(f"Performance comparison chart saved to: bond_angle_performance_comparison.png")

if __name__ == "__main__":
    # Check if input file is provided
    if len(sys.argv) > 1:
        gro_file = sys.argv[1]
    else:
        # Default to the file used in run_bond_angle_save_mat.py
        gro_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), '4xpreem.gro')
    
    # Check if number of runs is provided
    runs = 3
    if len(sys.argv) > 2:
        try:
            runs = int(sys.argv[2])
        except ValueError:
            print(f"Invalid number of runs '{sys.argv[2]}', using default: {runs}")
    
    # Check if subset size is provided
    subset_size = None
    if len(sys.argv) > 3:
        try:
            subset_size = int(sys.argv[3])
        except ValueError:
            print(f"Invalid subset size '{sys.argv[3]}', using full structure")
    
    # Run the comparison
    run_comparison(gro_file, runs, subset_size)
