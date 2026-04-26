import os
import time
import sys
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import self_distance_array, capped_distance
import atomipy as ap
from atomipy.distances import dist_matrix, neighbor_list_fast

def benchmark_vs_mda():
    base_dir = "DifferentSystemSizes"
    log_file = os.path.join(base_dir, "performance_comparison.log")
    targets = [2000, 4000, 10000, 25000, 50000, 100000]
    cutoff = 10.0
    results = []
    
    # Custom logger to print to both console and file
    class Logger:
        def __init__(self, filename):
            self.terminal = sys.stdout
            self.log = open(filename, "w")
        def write(self, message):
            self.terminal.write(message)
            self.log.write(message)
        def flush(self):
            self.terminal.flush()
            self.log.flush()

    sys.stdout = Logger(log_file)
    
    print("="*80)
    print("PERFORMANCE COMPARISON: atomipy vs MDAnalysis")
    print(f"Cutoff for sparse methods: {cutoff} Å")
    print(f"Log file: {log_file}")
    print("="*80)
    
    for t in targets:
        gro_file = os.path.join(base_dir, f"System{t}.gro")
        if not os.path.exists(gro_file):
            print(f"Skipping System{t}: File not found")
            continue
            
        print(f"\nBenchmarking System {t} atoms...")
        atoms, box_input = ap.import_auto(gro_file)
        box_dim, cell = ap.cell_utils.normalize_box(box_input)
        n_atoms = len(atoms)
        
        # Prepare coordinates for MDAnalysis
        coords = np.array([[a['x'], a['y'], a['z']] for a in atoms], dtype=np.float32)
        mda_box = np.array(cell, dtype=np.float32)
        
        # 1. Full Distance Matrix (Up to 25k)
        t_ap_full = "Skipped"
        t_mda_full = "Skipped"
        
        if n_atoms <= 26000: # Slightly over 25k to catch the ~25k target
            print(f"  [Full Matrix] N={n_atoms}")
            # atomipy
            start = time.perf_counter()
            _ = dist_matrix(atoms, box_dim)
            t_ap_full = time.perf_counter() - start
            print(f"    atomipy: {t_ap_full:.3f} s")
            
            # MDAnalysis
            start = time.perf_counter()
            _ = self_distance_array(coords, box=mda_box)
            t_mda_full = time.perf_counter() - start
            print(f"    MDAnalysis: {t_mda_full:.3f} s")
        
        # 2. Sparse / Capped Distance (Neighbor List)
        print(f"  [Sparse/Capped] N={n_atoms}, Cutoff={cutoff}")
        # atomipy
        start = time.perf_counter()
        _ = neighbor_list_fast(atoms, box_dim, cutoff=cutoff)
        t_ap_sparse = time.perf_counter() - start
        print(f"    atomipy: {t_ap_sparse:.3f} s")
        
        # MDAnalysis
        start = time.perf_counter()
        _ = capped_distance(coords, coords, max_cutoff=cutoff, box=mda_box)
        t_mda_sparse = time.perf_counter() - start
        print(f"    MDAnalysis: {t_mda_sparse:.3f} s")
        
        results.append({
            'target': t,
            'actual': n_atoms,
            'ap_full': t_ap_full,
            'mda_full': t_mda_full,
            'ap_sparse': t_ap_sparse,
            'mda_sparse': t_mda_sparse
        })
        
    # Output the table
    print("\n\n" + "="*100)
    print("PERFORMANCE SUMMARY: atomipy (ap) vs MDAnalysis (mda) [Time in Seconds]")
    print("="*100)
    
    header = "| Target | Actual | ap Full | mda Full | ap Sparse | mda Sparse | Ratio Sparse (ap/mda) |"
    sep = "|---|---|---|---|---|---|---|"
    print(header)
    print(sep)
    
    for r in results:
        full_ap = f"{r['ap_full']:.3f}" if isinstance(r['ap_full'], float) else str(r['ap_full'])
        full_mda = f"{r['mda_full']:.3f}" if isinstance(r['mda_full'], float) else str(r['mda_full'])
        sparse_ap = f"{r['ap_sparse']:.3f}"
        sparse_mda = f"{r['mda_sparse']:.3f}"
        ratio = f"{r['ap_sparse']/r['mda_sparse']:.2f}x"
        
        print(f"| {r['target']} | {r['actual']} | {full_ap} | {full_mda} | {sparse_ap} | {sparse_mda} | {ratio} |")

    # Output MATLAB Commands
    print("\n\n" + "="*80)
    print("MATLAB COMMANDS FOR atom_MATLAB_toolbox")
    print("="*80)
    print("% Run these commands in MATLAB with the atom_MATLAB_toolbox in your path")
    print(f"base_dir = '{base_dir}';")
    print("targets = [2000, 4000, 10000, 25000, 50000, 100000];")
    print("cutoff = 10.0;")
    print("\nfor t = targets")
    print("    filename = fullfile(base_dir, sprintf('System%d.gro', t));")
    print("    if exist(filename, 'file')")
    print("        fprintf('\\nBenchmarking System %d atoms...\\n', t);")
    print("        atom = import_atom(filename);")
    print("        Box_dim = Box_dim; % Box_dim is usually imported with import_atom as a global or local variable")
    print("        ")
    print("        if numel(atom) <= 25000")
    print("            tic;")
    print("            [dist_matrix, dx, dy, dz] = dist_matrix_particle_MATLAB(atom, Box_dim);")
    print("            t_full = toc;")
    print("            fprintf('  MATLAB Full Matrix: %.3f s\\n', t_full);")
    print("        end")
    print("        ")
    print("        tic;")
    print("        [bond_list, dist_list] = neighbor_list_particle_MATLAB(atom, Box_dim, 1.2, cutoff);")
    print("        t_sparse = toc;")
    print("        fprintf('  MATLAB Neighbor List (Sparse): %.3f s\\n', t_sparse);")
    print("    end")
    print("end")

    # Restore stdout
    sys.stdout.log.close()
    sys.stdout = sys.stdout.terminal

if __name__ == "__main__":
    benchmark_vs_mda()
