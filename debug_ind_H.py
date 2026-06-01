import atomipy as ap
struct_path = '/Users/miho0052/Dropbox/Coding/Windsurf/atominpython/Pyrophyllite_GII_0.071.gro'
atoms, Box = ap.import_gro(struct_path)
atoms, Box, _ = ap.replicate_system(atoms, Box, replicate=[6, 4, 3])
atoms = ap.minff(atoms, Box=Box)

ind_H = {i for i, atom in enumerate(atoms, 1) if atom.get('type', '').startswith('H')}
print(f"ind_H contains {len(ind_H)} atoms. First 10: {list(ind_H)[:10]}")
for i in [1, 26, 32]:
    print(f"Atom {i}: type='{atoms[i-1].get('type', '')}', fftype='{atoms[i-1].get('fftype', '')}'")
