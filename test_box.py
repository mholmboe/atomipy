def _isint(s):
    try:
        int(s)
        return True
    except:
        return False
        
def _isfloat(s):
    try:
        float(s)
        return True
    except:
        return False

def _is_gro_coord(line):
    fields = []
    fields.append(line[16:20].strip()) # atom number
    fields.append(line[21:28].strip()) # x coord
    fields.append(line[29:36].strip()) # y coord
    fields.append(line[37:44].strip()) # z coord
    if (all([f != '' for f in fields])): # check for empty fields
        return all([_isint(fields[0]), _isfloat(fields[1]), _isfloat(fields[2]), _isfloat(fields[3])])
    else:
        return 0

vals = [3.09600, 3.58632, 2.75697, 0.00000, 0.00000, 0.02252, 0.00000, -0.50907, -0.05454]

# Atomipy original
s1 = '   '.join(f"{val:.5f}" for val in vals) + "\n"
print("Original:", _is_gro_coord(s1))

# Gromacs format
s2 = "".join(f" {val:10.5f}" for val in vals) + "\n"
print("Gromacs:", _is_gro_coord(s2))

