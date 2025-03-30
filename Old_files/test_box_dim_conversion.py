#!/usr/bin/env python
"""
Test script to verify that the cell_utils.Box_dim2Cell implementation
is consistent with the previous box_dim_to_cell implementation.
"""

import numpy as np
import sys
import os

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# Import both implementations for comparison
from atomipy.cell_utils import Box_dim2Cell

# Define some test box dimensions
test_boxes = [
    # Orthogonal box
    [10.0, 10.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    # Triclinic box
    [10.0, 10.0, 10.0, 0.0, 0.0, 2.0, 0.0, 1.0, 1.5],
    # Realistic simulation box
    [60.5, 60.5, 60.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    # Complex triclinic box
    [45.2, 37.8, 52.1, 0.0, 0.0, 5.2, 0.0, 3.1, 2.7]
]

def original_box_dim_to_cell(Box_dim):
    """Original implementation from cell_list_dist_matrix.py"""
    lx, ly, lz = Box_dim[0], Box_dim[1], Box_dim[2]
    xy, xz, yz = Box_dim[5], Box_dim[7], Box_dim[8]
    
    a = lx
    b = np.sqrt(ly**2 + xy**2)
    c = np.sqrt(lz**2 + xz**2 + yz**2)
    
    alpha = np.degrees(np.arccos((xy*xz + ly*yz) / (b*c)))
    beta = np.degrees(np.arccos(xz / c))
    gamma = np.degrees(np.arccos(xy / b))
    
    return [a, b, c, alpha, beta, gamma]

def main():
    """Run tests to compare the original and new implementations"""
    print("Testing box dimension to cell parameter conversion...")
    print("-" * 70)
    print("Box Dimensions | Original Implementation | cell_utils.Box_dim2Cell")
    print("-" * 70)
    
    # Test each box
    for i, box in enumerate(test_boxes):
        # Get cell parameters using original implementation
        original_cell = original_box_dim_to_cell(box)
        
        # Get cell parameters using new implementation
        new_cell = Box_dim2Cell(box)
        
        # Print results
        box_str = f"Box {i+1}"
        orig_str = f"[{', '.join([f'{x:.3f}' for x in original_cell])}]"
        new_str = f"[{', '.join([f'{x:.3f}' for x in new_cell])}]"
        
        print(f"{box_str:15} | {orig_str:40} | {new_str:40}")
        
        # Check for significant differences
        abs_diff = np.abs(np.array(original_cell) - np.array(new_cell))
        if np.any(abs_diff > 1e-10):
            print(f"  WARNING: Differences detected in box {i+1}:")
            for j, (o, n, d) in enumerate(zip(original_cell, new_cell, abs_diff)):
                if d > 1e-10:
                    print(f"    Parameter {j+1}: {o:.6f} vs {n:.6f} (diff: {d:.6e})")
        else:
            print(f"  ✓ Implementations match for box {i+1}")
        
        print("-" * 70)
    
    print("\nTest completed!")

if __name__ == "__main__":
    main()
