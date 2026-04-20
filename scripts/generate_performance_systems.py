#!/usr/bin/env python3
"""
Generate mixed clay + ion + water systems of various sizes for performance benchmarking.
Standard recipe: Montmorillonite bilayer with Na counter-ions and SPCE water.
"""

import os
import sys
import numpy as np

# Add repo root to path
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

import time

import atomipy as ap

def build_performance_system(target_size, output_prefix):
    print(f"\n--- Generating system with target size: {target_size} atoms ---")
    
    timings = {}
    
    # Unit cell info
    input_conf = "Pyrophyllite_GII_0.071.pdb"
    atoms, cell = ap.import_auto(input_conf)
    box_dim = ap.Cell2Box_dim(cell)
    
    # Calculate replication factors
    base_uc_per_layer = 24
    base_total = 4112
    atoms_per_uc_layer = base_total / base_uc_per_layer # ~171.3
    
    target_uc_per_layer = target_size / atoms_per_uc_layer
    nx = int(np.round(np.sqrt(target_uc_per_layer)))
    if nx < 1: nx = 1
    ny = int(np.round(target_uc_per_layer / nx))
    if ny < 1: ny = 1
    
    actual_uc_per_layer = nx * ny
    print(f"Using replication: {nx}x{ny}x1 (Total unit cells in bilayer: {actual_uc_per_layer * 2})")
    
    # Proportional scaling
    n_subst_target = int(np.round(actual_uc_per_layer * (16/24)))
    n_waters = int(np.round(actual_uc_per_layer * (360/24)))
    
    print(f"Targeting: {n_subst_target} Na per interlayer, {n_waters} H2O per interlayer")
    
    # 1. Replication Setup
    start = time.perf_counter()
    atoms = ap.center(atoms, box_dim, dim="z")
    atoms = ap.translate(atoms, [0, 0, -box_dim[2]/2])
    atoms_rep, box_rep, _ = ap.replicate_system(atoms, box_dim, [nx, ny, 1])
    timings['replication'] = time.perf_counter() - start
    
    # 2. Substitution (with up to 5 retries)
    timings['substitution'] = 0
    actual_subst1 = 0
    actual_subst2 = 0
    
    for attempt in range(5):
        print(f"Substitution attempt {attempt+1}/5 for Layer 1...")
        sub_start = time.perf_counter()
        # Use lo_limit/hi_limit to avoid auto-split issues with single flat layer
        mmt1_atoms, _, _ = ap.substitute(
            atoms_rep, box_rep, num_oct_subst=n_subst_target,
            o1_type='Alo', o2_type='Mgo', min_o2o2_dist=5.2,
            lo_limit=-100, hi_limit=100 # Large limits just to bypass auto-split logic
        )
        actual_subst1 = sum(1 for a in mmt1_atoms if a['type'] == 'Mgo')
        if actual_subst1 >= n_subst_target or attempt == 4:
            timings['substitution'] += (time.perf_counter() - sub_start)
            break
        print(f"  Warning: Only placed {actual_subst1}/{n_subst_target} sites. Retrying...")

    for attempt in range(5):
        print(f"Substitution attempt {attempt+1}/5 for Layer 2...")
        sub_start = time.perf_counter()
        mmt2_atoms, _, _ = ap.substitute(
            atoms_rep, box_rep, num_oct_subst=n_subst_target,
            o1_type='Alo', o2_type='Mgo', min_o2o2_dist=5.2,
            lo_limit=-100, hi_limit=100
        )
        actual_subst2 = sum(1 for a in mmt2_atoms if a['type'] == 'Mgo')
        if actual_subst2 >= n_subst_target or attempt == 4:
            timings['substitution'] += (time.perf_counter() - sub_start)
            break
        print(f"  Warning: Only placed {actual_subst2}/{n_subst_target} sites. Retrying...")

    for atom in mmt1_atoms:
        atom['resname'] = 'MMT1'
        atom['molid'] = 1
    for atom in mmt2_atoms:
        atom['resname'] = 'MMT2'
        atom['molid'] = 2
        
    # Combine layers for next steps
    basal_spacing = 22.0
    full_box = box_rep.copy()
    full_box[2] = 2 * basal_spacing
    mmt2_atoms = ap.translate(mmt2_atoms, [0, 0, basal_spacing])
    combined_atoms = ap.update(mmt1_atoms, mmt2_atoms)
    combined_atoms = ap.wrap(combined_atoms, full_box, return_type='cartesian')
    
    # 3. Add Ions (matching actual successful substitutions per layer)
    start = time.perf_counter()
    lower_limits = [0, 0, 0, full_box[0], full_box[1], basal_spacing]
    upper_limits = [0, 0, basal_spacing, full_box[0], full_box[1], full_box[2]]
    
    print(f"Ionizing: adding {actual_subst1} Na to lower and {actual_subst2} Na to upper interlayer")
    na_lower = ap.ionize(ion_type='Na', resname='Na', limits=lower_limits, 
                         num_ions=actual_subst1, Box=full_box, min_distance=3.0, solute_atoms=combined_atoms)
    na_upper = ap.ionize(ion_type='Na', resname='Na', limits=upper_limits, 
                         num_ions=actual_subst2, Box=full_box, min_distance=3.0, solute_atoms=combined_atoms)
    
    all_atoms = ap.update(combined_atoms, na_lower, na_upper)
    all_atoms = ap.wrap(all_atoms, full_box, return_type='cartesian')
    timings['addition of ions'] = time.perf_counter() - start
    
    # ... rests of steps ...
    # (keeping them as they were, but ensuring they are wrapped correctly in the function)
    
    # 4. Solvation
    start = time.perf_counter()
    water_lower = ap.solvate(limits=lower_limits, max_solvent=n_waters, 
                             solute_atoms=all_atoms, Box=full_box, min_distance=2.0, solvent_type='spce')
    water_upper = ap.solvate(limits=upper_limits, max_solvent=n_waters, 
                             solute_atoms=all_atoms, Box=full_box, min_distance=2.0, solvent_type='spce')
    
    final_atoms = ap.update(all_atoms, water_lower, water_upper)
    timings['solvation'] = time.perf_counter() - start
    
    # 5. Atomtype assignment
    start = time.perf_counter()
    if len(full_box) > 6:
        cell = ap.Box_dim2Cell(full_box)
        minff_atoms = ap.minff(final_atoms, cell)
    else:
        minff_atoms = ap.minff(final_atoms, full_box)
    timings['atomtype assignment'] = time.perf_counter() - start
        
    print(f"Final atom count: {len(minff_atoms)}")
    
    # 6-9 Writing files
    start = time.perf_counter()
    ap.write_gro(minff_atoms, Box=full_box, file_path=f"{output_prefix}.gro")
    timings['writing gro'] = time.perf_counter() - start
    
    start = time.perf_counter()
    ap.write_pdb(minff_atoms, Box=full_box, file_path=f"{output_prefix}.pdb")
    timings['writing pdb'] = time.perf_counter() - start
    
    start = time.perf_counter()
    if len(full_box) > 6:
        cell = ap.Box_dim2Cell(full_box)
        ap.write_itp(minff_atoms, cell, f"{output_prefix}.itp")
    else:
        ap.write_itp(minff_atoms, full_box, f"{output_prefix}.itp")
    timings['writing gromacs itp'] = time.perf_counter() - start
        
    start = time.perf_counter()
    if len(full_box) > 6:
        cell = ap.Box_dim2Cell(full_box)
        ap.write_lmp(minff_atoms, Box=cell, file_path=f"{output_prefix}.data", detect_bimodal=True)
    else:
        ap.write_lmp(minff_atoms, Box=full_box, file_path=f"{output_prefix}.data", detect_bimodal=True)
    timings['writing lammps data'] = time.perf_counter() - start
    
    # 10. Bond Stats
    start = time.perf_counter()
    ap.get_structure_stats(minff_atoms, Box=full_box, log_file=f"{output_prefix}.log")
    timings['get bond statistics and write a log file'] = time.perf_counter() - start
    
    timings['total'] = sum(v for k, v in timings.items() if k != 'total' and isinstance(v, (int, float)))
    timings['actual_size'] = len(minff_atoms)
    
    print(f"Saved: {output_prefix}.gro, .pdb, .itp, .data, .log")
    return timings

if __name__ == "__main__":
    targets = [2000, 4000, 10000, 25000, 50000, 100000]
    all_results = []
    
    for t in targets:
        filename = os.path.join("DifferentSystemSizes_v2", f"System{t}")
        res = build_performance_system(t, filename)
        res['target'] = t
        all_results.append(res)
        
    print("\n" + "="*80)
    print("PERFORMANCE BENCHMARK SUMMARY (Run Times in Seconds)")
    print("="*80)
    
    categories = [
        'replication', 'substitution', 'addition of ions', 'solvation',
        'atomtype assignment', 'writing gro', 'writing pdb', 
        'writing gromacs itp', 'writing lammps data', 
        'get bond statistics and write a log file'
    ]
    
    # Headers - shorter versions for table
    short_cats = ['Repl', 'Subst', 'Ions', 'Solv', 'Atomtype', 'GRO', 'PDB', 'ITP', 'LAMMPS', 'Stats']
    
    header = "| Target | Actual | " + " | ".join(short_cats) + " | **Total** |"
    sep = "|---" + "|---"*len(short_cats) + "|---:|"
    print(header)
    print(sep)
    
    for res in all_results:
        row = f"| {res['target']} | {res['actual_size']} | "
        for cat in categories:
            row += f"{res[cat]:.3f} | "
        row += f"**{res['total']:.3f}** |"
        print(row)
    
    print("\nAll systems generated and benchmarked successfully!")
    
    print("\nAll systems generated and benchmarked successfully!")
