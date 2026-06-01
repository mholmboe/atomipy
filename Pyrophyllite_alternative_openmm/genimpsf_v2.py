#!/usr/bin/python
import sys
import os
wd = os.getcwd()

fname= os.path.join(wd, sys.argv[1].lower())
inname = fname[:fname.rfind('.')]
inname = inname.lower()
gname = os.path.join(wd,"image.psf")
output = inname + '_imc.psf'
#######################################################################
#
# read original psf file
#
#######################################################################
read_atom = False
read_bond = False
read_angl = False
read_dih  = False
read_auto = False
f = open(fname)
atoms_orig = {}
lines = {'HEADER': [], 'TITLE':[], 'ATOM':[], 'BOND':[], 'ANGL':[], 'DIHE':[], 'AUTO':[], 'OTHER':[]}
atidxdic = {}
iarr_orig = []
jarr_orig = []
iarr2_orig = []
jarr2_orig = []
karr2_orig = []
iarr3_orig = []
jarr3_orig = []
karr3_orig = []
larr3_orig = []
auto_arr   = []
nbond_orig = 0
nangl_orig = 0
ndih_orig = 0
nbond_last = 0
is_after_title = False
is_after_atom = False
is_after_bond = False
for line in f:
    if not line.strip() and read_atom:
        read_atom = False
        is_after_atom = True
    if not line.strip() and read_bond:
        read_bond = False
        is_after_bond = True
    if not line.strip() and read_angl:
        read_angl = False
    if not line.strip() and read_dih:
        read_dih = False
    if 'PSF EXT CMAP CHEQ XPLOR' in line:
        lines['HEADER'].append(line)
        continue
    if '!NATOM' in line:
        read_atom = True
        lines['ATOM'].append(line)
        is_after_title = True
        continue
    elif '!NBOND' in line:
        read_bond = True
        nbond_orig = int(line.split()[0])
        continue
    elif '!NTHETA' in line:
        read_angl = True
        nangl_orig = int(line.split()[0])
        continue
    elif '!NPHI' in line:
        read_dih = True
        ndih_orig = int(line.split()[0])
        continue
    elif '!AUTOGEN' in line:
        read_auto = True
        continue
    if read_atom:
        index = int(line.split()[0])
        resname = line.split()[1]
        ires = int(line.split()[2])
        atomname = line.split()[4]
        atoms_orig[(ires, resname, atomname)] = index
        lines['ATOM'].append(line)
        continue
    if read_bond: 
        if len(line.split()) == 8:
            lines['BOND'].append(line)
        else:
            for x,element in enumerate(line.split()):
                nbond_last += 1
                if x % 2 == 0:
                    iarr_orig.append(element)
                else:
                    jarr_orig.append(element)
        continue
    if read_angl:
        if len(line.split()) == 9:
            lines['ANGL'].append(line)
        else:
            for x,element in enumerate(line.split()):
                if x % 3 == 0:
                    iarr2_orig.append(element)
                elif x % 3 == 1:
                    jarr2_orig.append(element)
                else:
                    karr2_orig.append(element)
        continue
    if read_dih:
        if len(line.split()) == 8:
            lines['DIHE'].append(line)
        else:
            for x,element in enumerate(line.split()):
                if x % 4 == 0:
                    iarr3_orig.append(element)
                elif x % 4 == 1:
                    jarr3_orig.append(element)
                elif x % 4 == 2:
                    karr3_orig.append(element)
                else:
                    larr3_orig.append(element)
        continue
    if read_auto:
        if len(line.split()) == 8 and line.split()[0] != '0':
            lines['AUTO'].append(line)
        #else:
        #    for x,element in enumerate(line.split()):
        #        auto_arr.append(element)
        continue
    if is_after_title:
        lines['OTHER'].append(line)
    else:
        lines['TITLE'].append(line)
f.close()

if lines['HEADER'] and 'AUTOG' in lines['HEADER'][0]:
    lines['HEADER'] = lines['HEADER'][0][:-7] + '\n'

# 
# read image file
#
#######################################################################
f = open(gname)
read_num = False
read_bond = False
read_angl = False
read_dih = False
atoms_image = {}
title = []
numarr = []
nimbond = 0 # bond
iarr = [] # bond
jarr = []
nang = 0 # angle
angarr = [] # angle
iarr2 = []
jarr2 = []
karr2 = []
iarr3 = []
jarr3 = []
karr3 = []
larr3 = []
ndih = 0
diharr = []
for line in f:
    if '*' in line:
        title.append(line)
    if 'NTRANS' in line:
        read_num = True
        continue
    if read_num:
        numarr.append(line)
        nimbond = int(line.split()[1])
        nang = line.split()[2]
        ndih = line.split()[3]
        read_num = False
    if 'BOND ARRAY' in line:
        read_bond = True
        continue
    if read_bond and len(line.split()) != False:
        resi_name = line.split()[2]
        resi = line.split()[3]
        at1  = line.split()[4]
        resj_name = line.split()[6]
        resj = line.split()[7]
        at2  = line.split()[8]
        try:
            iarr.append(atoms_orig[int(resi), resi_name, at1])
            jarr.append(atoms_orig[int(resj), resj_name, at2])
        except KeyError:
            continue
    if len(line.split()) == False:
       read_bond = False
    if 'ANGLE ARRAY' in line:
        read_angl = True
        continue
    if read_angl and len(line.split()) != False:
        resi_name = line.split()[2]
        resi = line.split()[3]
        at1  = line.split()[4]
        resj_name = line.split()[6]
        resj = line.split()[7]
        at2  = line.split()[8]
        resk_name = line.split()[10]
        resk = line.split()[11]
        at3  = line.split()[12]
        try:
            iarr2.append(atoms_orig[int(resi), resi_name, at1])
            jarr2.append(atoms_orig[int(resj), resj_name, at2])
            karr2.append(atoms_orig[int(resk), resk_name, at3])
        except KeyError:
            continue
    if len(line.split()) == False:
       read_angl = False
    if 'DIHEDRAL ARRAY' in line:
        read_dih = True
        continue
    if read_dih and len(line.split()) != False:
        resi_name = line.split()[2]
        resi = line.split()[3]
        at1  = line.split()[4]
        resj_name = line.split()[6]
        resj = line.split()[7]
        at2  = line.split()[8]
        resk_name = line.split()[10]
        resk = line.split()[11]
        at3  = line.split()[12]
        resl_name = line.split()[14]
        resl = line.split()[15]
        at4  = line.split()[16]
        try:
            iarr3.append(atoms_orig[int(resi), resi_name, at1])
            jarr3.append(atoms_orig[int(resj), resj_name, at2])
            karr3.append(atoms_orig[int(resk), resk_name, at3])
            larr3.append(atoms_orig[int(resl), resl_name, at4])
        except KeyError:
            continue
    if len(line.split()) == False:
       read_dih = False
#######################################################################
#
# convert image index into original index
#
#######################################################################
new_iarr = []
new_jarr = []
new_iarr2 = []
new_jarr2 = []
new_karr2 = []
new_iarr3 = []
new_jarr3 = []
new_karr3 = []
new_larr3 = []
# bond
for xx in iarr:
    new_iarr.append(int(xx))
for xx in jarr:
    new_jarr.append(int(xx))
# angle
for xx in iarr2:
    new_iarr2.append(int(xx))
for xx in jarr2:
    new_jarr2.append(int(xx))
for xx in karr2:
    new_karr2.append(int(xx))
# dihe
for xx in iarr3:
    new_iarr3.append(int(xx))
for xx in jarr3:
    new_jarr3.append(int(xx))
for xx in karr3:
    new_karr3.append(int(xx))
for xx in larr3:
    new_larr3.append(int(xx))
#######################################################################
#
# write output psf file
#
#######################################################################
f = open(output, 'w')
f.writelines(lines['HEADER'])
f.writelines(lines['TITLE'])
f.writelines(lines['ATOM'])
# combine between last original array and image array
iarr_orig += new_iarr
jarr_orig += new_jarr 
iarr2_orig += new_iarr2
jarr2_orig += new_jarr2
karr2_orig += new_karr2
# remove duplicate bond
iarr3_orig += new_iarr3
jarr3_orig += new_jarr3
karr3_orig += new_karr3
larr3_orig += new_larr3
remove_duplicate = set()
for x in range(len(iarr_orig)):
    i = int(iarr_orig[x])
    j = int(jarr_orig[x])
    remove_duplicate.add((i,j))
nbond_image = len(remove_duplicate)
# remove duplicate angle

# write nbond header
remove_duplicate = set()
xbond = 0
tot_nbond = nbond_orig + nimbond
nbond_orig += nbond_image
nbond_orig -= nbond_last
norig_bond = 0
for line in lines['BOND']:
    for x,xbond in enumerate(line.split()):
        if x % 2 == 0:
            i = int(xbond)
        else:
            j = int(xbond)
            remove_duplicate.add((i,j))
            norig_bond += 1
 
#f.write('\n%9d !NBOND: bonds\n' % nbond_orig)
f.write('\n%9d !NBOND: bonds\n' % (norig_bond+nbond_image))
f.writelines(lines['BOND'])

# write bond information
xbond = 0
remove_duplicate = set()
for x in range(len(iarr_orig)):
    if xbond > 3:
        f.write('\n')
        xbond = 0
    i = int(iarr_orig[x])
    j = int(jarr_orig[x])
    if (i,j) in remove_duplicate: continue
    if (j,i) in remove_duplicate: continue
    f.write('%10d%10d' % (i,j))
    remove_duplicate.add((i,j))
    xbond += 1

if xbond != 0:
    f.write('\n')

remove_duplicate2 = set()
norig_angle = 0
for line in lines['ANGL']:
    for x,xangl in enumerate(line.split()):
        if x % 3 == 0:
            i = int(xangl)
        elif x % 3 == 1:
            j = int(xangl)
        else:
            k = int(xangl)
            remove_duplicate2.add((i,j,k))
            norig_angle += 1

remove_duplicate3 = set()
norig_dih = 0
for line in lines['DIHE']:
    for x,xdih in enumerate(line.split()):
        if x % 4 == 0:
            i = int(xdih)
        elif x % 4 == 1:
            j = int(xdih)
        elif x % 4 == 2:
            k = int(xdih)
        else:
            l = int(xdih)
            remove_duplicate3.add((i,j,k,l))
            norig_dih += 1

# write angle information
xangl = 0
nangl_image = 0
lines_angle = []
for x in range(len(iarr2_orig)):
    if xangl > 2:
        #f.write('\n')
        lines_angle.append('\n')
        xangl = 0
    try:
        i = int(iarr2_orig[x])
        j = int(jarr2_orig[x])
        k = int(karr2_orig[x])
    except IndexError:
        continue
    if (i,j,k) in remove_duplicate2: continue
    if (k,j,i) in remove_duplicate2: continue
    remove_duplicate2.add((i,j,k))
    lines_angle.append('%10d%10d%10d' % (i,j,k))
    nangl_image += 1
    xangl += 1

if xangl != 0:
    lines_angle.append('\n');

nangl_image += norig_angle
f.write('\n%9d !NTHETA: angles\n' % nangl_image)
f.writelines(lines['ANGL'])
f.writelines(lines_angle)

# write dihedral information
xdih = 0
ndih_image = 0
lines_dih = []
for x in range(len(iarr3_orig)):
    if xdih > 1:
        lines_dih.append('\n')
        xdih = 0
    try:
        i = int(iarr3_orig[x])
        j = int(jarr3_orig[x])
        k = int(karr3_orig[x])
        l = int(larr3_orig[x])
    except IndexError:
        continue
    if (i,j,k,l) in remove_duplicate3: continue
    if (l,k,j,i) in remove_duplicate3: continue
    remove_duplicate3.add((i,j,k,l))
    lines_dih.append('%10d%10d%10d%10d' % (i,j,k, l))
    ndih_image += 1
    xdih += 1

ndih_image += norig_dih
f.write('\n%9d !NPHI: dihedrals\n' % ndih_image)
f.writelines(lines['DIHE'])
if len(lines_dih) > 0:
    if lines_dih[-1]=='\n':
        lines_dih = lines_dih[:-1]
f.writelines(lines_dih)

f.writelines(lines['OTHER'][2:])

#f.write('%10d !AUTOGEN\n' % 0)
#f.writelines(lines['AUTO'])
## add extra 0 for autogen
#auto_arr = ['0'] * tot_nbond
#print(tot_nbond)
#xauto = 0
#for x in range(len(auto_arr)):
#    if xauto > 7:
#        f.write('\n')
#        xauto = 0
#    i = int(auto_arr[x])
#    f.write('%10d' % (i))
#    xauto += 1
#if xauto != 0:
#    f.write('\n')
#

f.close()
#######################################################################

