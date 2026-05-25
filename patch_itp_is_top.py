import sys

with open('atomipy/write_top.py', 'r') as f:
    content = f.read()

old_header = """        # Write header
        f.write("; Gromacs topology file\\n")"""

new_header = """        # Write header
        f.write("; Gromacs topology file\\n")
        if is_top:
            f.write('#include "min.ff/forcefield.itp"\\n\\n')"""

content = content.replace(old_header, new_header)

old_footer = """        # End of moleculetype
        f.write("\\n")
        
        # Add position restraints block placeholder"""

new_footer = """        # End of moleculetype
        f.write("\\n")
        
        if is_top:
            f.write("\\n[ system ]\\n")
            f.write(f"; Name\\n{molecule_name}\\n\\n")
            f.write("[ molecules ]\\n")
            f.write(f"; Compound        #mols\\n{molecule_name}               1\\n\\n")
            
        # Add position restraints block placeholder"""

content = content.replace(old_footer, new_footer)

with open('atomipy/write_top.py', 'w') as f:
    f.write(content)

