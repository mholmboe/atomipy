import os
import time
import numpy as np
import atomipy as ap
from atomipy.dist_matrix import dist_matrix
from atomipy.cell_list_dist_matrix import cell_list_dist_matrix, neighbor_list_fast

def benchmark_distance_matrix():
    # Use the systems from DifferentSystemSizes_v2
    base_dir = "DifferentSystemSizes_v2"
    targets = [2000, 4000, 10000, 25000, 50000, 100000]
    results = []
    
    print("="*60)
    print("DISTANCE MATRIX PERFORMANCE BENCHMARK")
    print("Direct (O(N^2)) vs Cell List (O(N))")
    print("="*60)
    
    for t in targets:
        gro_file = os.path.join(base_dir, f"System{t}.gro")
        if not os.path.exists(gro_file):
            print(f"Skipping System{t}: File not found")
            continue
            
        print(f"\nBenchmarking System {t} atoms...")
        atoms, box_input = ap.import_auto(gro_file)
        # normalize_box returns (Box_dim, Cell)
        box_dim, cell = ap.cell_utils.normalize_box(box_input)
        n_atoms = len(atoms)
        
        # 1. Direct Method (Up to 25k)
        t_direct = None
        if t <= 25000:
            print(f"  Running Direct Method (N={n_atoms})...")
            try:
                start = time.perf_counter()
                _ = dist_matrix(atoms, box_dim)
                t_direct = time.perf_counter() - start
                print(f"    Direct: {t_direct:.3f} s")
            except MemoryError:
                print("    Direct: Memory Error!")
                t_direct = "MemoryError"
        else:
            print("  Direct Method: Skipped (too large)")
            t_direct = "Skipped"
            
        # 2. Cell List Method (Full Matrix where possible)
        t_cell_full = None
        if t <= 25000:
            print(f"  Running Cell List Method (Full Matrix, N={n_atoms})...")
            try:
                start = time.perf_counter()
                _ = cell_list_dist_matrix(atoms, box_dim, cutoff=15.0) # Larger cutoff to make it work harder
                t_cell_full = time.perf_counter() - start
                print(f"    Cell List (Full): {t_cell_full:.3f} s")
            except MemoryError:
                print("    Cell List (Full): Memory Error!")
                t_cell_full = "MemoryError"
        else:
            t_cell_full = "Skipped"

        # 3. Cell List Method (Sparse/Neighbor List - Scales to 100k)
        print(f"  Running Cell List Method (Sparse, N={n_atoms})...")
        start = time.perf_counter()
        _ = neighbor_list_fast(atoms, box_dim, cutoff=15.0)
        t_cell_sparse = time.perf_counter() - start
        print(f"    Cell List (Sparse): {t_cell_sparse:.3f} s")
        
        results.append({
            'target': t,
            'actual': n_atoms,
            'direct': t_direct,
            'cell_full': t_cell_full,
            'cell_sparse': t_cell_sparse
        })
        
    # Output the table
    print("\n\n" + "="*80)
    print("DISTANCE CALCULATION PERFORMANCE SUMMARY (Seconds)")
    print("="*80)
    
    header = "| Target Size | Actual Atoms | Direct O(N²) | Cell List (Full NxN) | Cell List (Sparse) |"
    sep = "|---|---|---|---|---|"
    print(header)
    print(sep)
    
    for r in results:
        direct_str = f"{r['direct']:.3f}" if isinstance(r['direct'], float) else str(r['direct'])
        full_str = f"{r['cell_full']:.3f}" if isinstance(r['cell_full'], float) else str(r['cell_full'])
        sparse_str = f"{r['cell_sparse']:.3f}" if isinstance(r['cell_sparse'], float) else str(r['cell_sparse'])
        
        print(f"| {r['target']} | {r['actual']} | {direct_str} | {full_str} | {sparse_str} |")

if __name__ == "__main__":
    benchmark_distance_matrix()
