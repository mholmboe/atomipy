import sys

with open('atomipy/write_top.py', 'r') as f:
    content = f.read()

target_start = "    # Filter Angle_index to only include angles defined in the force field parameters (e.g. MINFF/CLAYFF)"
target_end = "            print(\"write_psf: Warning - Could not determine valid angles from parameter file, keeping all angles.\")"

start_idx = content.find(target_start)
end_idx = content.find(target_end) + len(target_end)

if start_idx != -1 and end_idx != -1:
    replacement = """    # Filter Angle_index to only include angles with at least one hydrogen atom (since CLAYFF/MINFF only define M-O-H and H-O-H angles)
    if Angle_index is not None and len(Angle_index) > 0:
        total_angles = len(Angle_index)
        filtered_angles = []
        for angle in Angle_index:
            a1, a2, a3 = int(angle[0]), int(angle[1]), int(angle[2])
            if a1 in ind_H or a2 in ind_H or a3 in ind_H:
                filtered_angles.append(angle)
        Angle_index = filtered_angles
        print(f"write_itp: Filtered to {len(Angle_index)} hydrogen-containing angles (from {total_angles} total angles)")"""
    new_content = content[:start_idx] + replacement + content[end_idx:]
    with open('atomipy/write_top.py', 'w') as f:
        f.write(new_content)
    print("Patched successfully")
else:
    print("Could not find target strings")
