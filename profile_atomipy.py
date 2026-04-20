#!/usr/bin/env python3
"""Profile test_atomipy.py for speed and memory usage."""
import time
import tracemalloc
import sys
import os

# Setup paths
REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, REPO_ROOT)

import atomipy as ap

input_file = sys.argv[1] if len(sys.argv) > 1 else "Clay_system.gro"

# Start memory tracking
tracemalloc.start()

timings = {}

def timed(label, func, *args, **kwargs):
    """Run func and record wall-clock time + peak memory delta."""
    snapshot_before = tracemalloc.take_snapshot()
    mem_before = sum(stat.size for stat in snapshot_before.statistics('filename'))
    t0 = time.perf_counter()
    result = func(*args, **kwargs)
    t1 = time.perf_counter()
    snapshot_after = tracemalloc.take_snapshot()
    mem_after = sum(stat.size for stat in snapshot_after.statistics('filename'))
    elapsed = t1 - t0
    mem_delta_mb = (mem_after - mem_before) / (1024 * 1024)
    timings[label] = {"time": elapsed, "mem_delta_mb": mem_delta_mb}
    print(f"  [{label}] {elapsed:.3f}s  | mem delta: {mem_delta_mb:+.2f} MB")
    return result

print(f"\n{'='*70}")
print(f"  PROFILING atomipy with {input_file}")
print(f"{'='*70}\n")

# 1. Import
atoms, Box_dim = timed("import_auto", ap.import_auto, input_file)
print(f"  Loaded {len(atoms)} atoms\n")

# 2. Element assignment
atoms = timed("element", ap.element, atoms)

# 3. find_H2O
SOL, noSOL = timed("find_H2O", ap.find_H2O, atoms, Box_dim)
print(f"  SOL: {len(SOL)}, noSOL: {len(noSOL)}")

# 4. assign_resname
noSOL = timed("assign_resname", ap.assign_resname, noSOL)

# 5. Separate
IONS = [a for a in noSOL if a.get('resname') == 'ION']
MIN = [a for a in noSOL if a.get('resname') == 'MIN']
print(f"  IONS: {len(IONS)}, MIN: {len(MIN)}")

# 6. molecule
if MIN:
    MIN = timed("molecule", ap.molecule, MIN, molid=1, resname='MIN')

# 7. update
System = timed("update", ap.update, MIN, IONS, SOL)
print(f"  System: {len(System)} atoms")

# 8. minff (the big one)
System = timed("minff", ap.minff, System, Box_dim)

# 9. Extract MIN
MIN = [a for a in System if a.get('resname') == 'MIN']

# 10. write_itp
timed("write_itp", ap.write_itp, MIN, Box=Box_dim, file_path='minff.itp')

# 11. write_psf
timed("write_psf", ap.write_psf, System, Box=Box_dim, file_path='minff.psf', detect_bimodal=True, max_angle=150)

# 12. load_forcefield
ff_params = timed("load_forcefield", ap.load_forcefield, 'GMINFF/gminff_all.json', blocks=['GMINFF_k500', 'OPC3_HFE_LM', 'OPC3'])

# 13. write_lmp
timed("write_lmp", ap.write_lmp, System, Box=Box_dim, file_path='minff.data', forcefield=ff_params, detect_bimodal=True)

# 14. write_lmp (no150angles) 
timed("write_lmp_no150", ap.write_lmp, System, Box=Box_dim, file_path='minff_no150angles.data', forcefield=ff_params, detect_bimodal=True, max_angle=150)

# 15. write_gro
timed("write_gro", ap.write_gro, System, Box=Box_dim, file_path='preem.gro')

# 16. write_pdb
timed("write_pdb", ap.write_pdb, System, Box=Box_dim, file_path='preem.pdb')

# 17. get_structure_stats
timed("get_structure_stats", ap.get_structure_stats, System, Box=Box_dim)

# Summary
print(f"\n{'='*70}")
print(f"  SUMMARY")
print(f"{'='*70}")
print(f"  {'Step':<25} {'Time (s)':>10} {'Mem Delta (MB)':>15}")
print(f"  {'-'*50}")

total_time = 0
for label, data in timings.items():
    total_time += data['time']
    print(f"  {label:<25} {data['time']:>10.3f} {data['mem_delta_mb']:>15.2f}")

print(f"  {'-'*50}")
print(f"  {'TOTAL':<25} {total_time:>10.3f}")

# Peak memory
current, peak = tracemalloc.get_traced_memory()
tracemalloc.stop()
print(f"\n  Current memory: {current / (1024*1024):.1f} MB")
print(f"  Peak memory:    {peak / (1024*1024):.1f} MB")
print(f"{'='*70}\n")
