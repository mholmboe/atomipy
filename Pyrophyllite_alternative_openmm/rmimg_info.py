#!/usr/bin/python
import sys
import os

wd = os.getcwd()

#fname= sys.argv[1].lower()
fname= os.path.join(wd, sys.argv[1].lower())
#fname= 'test.psf'
inname = fname.split('.')[0]
#inname = inname.split('_')[0]+'_'+inname.split('_')[1]+'_'+inname.split('_')[2]
#filename = inname + ".psf"
gname = "image.psf"
output = os.path.join(wd, sys.argv[2].lower())

#######################################################################
#
# read original psf file
#
#######################################################################
read_atom = False
read_bond = False
read_angl = False
f = open(fname, 'r')
atoms_orig = {}
lines = {'TITLE':[], 'ATOM':[], 'BOND':[], 'ANGL':[], 'OTHER':[]}
atidxdic = {}
iarr_orig = []
jarr_orig = []
iarr2_orig = []
jarr2_orig = []
karr2_orig = []
nbond_orig = 0
nangl_orig = 0
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
    if read_atom:
        index = int(line.split()[0])
        ires = int(line.split()[2])
        atomname = line.split()[4]
        atoms_orig[(ires, atomname)] = index
        lines['ATOM'].append(line)
        continue
    if read_bond: 
        for x,element in enumerate(line.split()):
            nbond_last += 1
            if x % 2 == 0:
                iarr_orig.append(element)
            else:
                jarr_orig.append(element)
        continue
    if read_angl:
        for x,element in enumerate(line.split()):
            if x % 3 == 0:
                iarr2_orig.append(element)
            elif x % 3 == 1:
                jarr2_orig.append(element)
            else:
                karr2_orig.append(element)
        continue
    if is_after_title:
        lines['OTHER'].append(line)
    else:
        lines['TITLE'].append(line)
f.close()
#print(len(iarr2_orig))
#print(len(jarr2_orig))
#print(len(karr2_orig))
#input()

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
nbond = 0 # bond
iarrrm = [] # bond
jarrrm = []
nang = 0 # angle
angarr = [] # angle
iarrrm2 = []
jarrrm2 = []
karrrm2 = []
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
        nbond = line.split()[1]
        nang = line.split()[2]
        ndih = line.split()[3]
        read_num = False
        #print(nbond)
    if 'BOND ARRAY' in line:
        read_bond = True
        continue
    if read_bond and len(line.split()) != False:
        resi = line.split()[3]
        at1  = line.split()[4]
        resj = line.split()[7]
        at2  = line.split()[8]
        try:
            iarrrm.append(atoms_orig[int(resi), at1])
            jarrrm.append(atoms_orig[int(resj), at2])
        except KeyError:
            continue
    if len(line.split()) == False:
       read_bond = False
    if 'ANGLE ARRAY' in line:
        read_angl = True
        continue
    if read_angl and len(line.split()) != False:
        #print(line)
        resi = line.split()[3]
        at1  = line.split()[4]
        resj = line.split()[7]
        at2  = line.split()[8]
        resk = line.split()[11]
        at3  = line.split()[12]
        try:
            iarrrm2.append(atoms_orig[int(resi), at1])
            jarrrm2.append(atoms_orig[int(resj), at2])
            karrrm2.append(atoms_orig[int(resk), at3])
        except KeyError:
            continue
    if len(line.split()) == False:
       read_angl = False
    if 'DIHEDRAL ARRAY' in line:
        read_dih = True
        continue
    if read_dih and len(line.split()) != False:
        resi = line.split()[3]
        at1  = line.split()[4]
        resj = line.split()[7]
        at2  = line.split()[8]
        resk = line.split()[11]
        at3  = line.split()[12]
        resl = line.split()[15]
        at4  = line.split()[16]
        #print(resi, at1, resl, at4)
    if len(line.split()) == False:
       read_dih = False
f.close()
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
# bond
trmbn = []
trman = []
for i in range(len(iarrrm)):
    trmbn.append([iarrrm[i], jarrrm[i]])
    trmbn.append([jarrrm[i], iarrrm[i]])
for k in range(len(iarr_orig)):
    temp = [int(iarr_orig[k]), int(jarr_orig[k])]
    if temp not in trmbn:
       new_iarr.append(iarr_orig[k])
       new_jarr.append(jarr_orig[k])
for m in range(len(iarrrm2)):
    try:
        trman.append([iarrrm2[m], jarrrm2[m], karrrm2[m]])
        trman.append([karrrm2[m], jarrrm2[m], iarrrm2[m]])
    except IndexError:
        continue
for k in range(len(iarr2_orig)):
    temp = [int(iarr2_orig[k]), int(jarr2_orig[k]), int(karr2_orig[k])]
    if temp not in trman:
       new_iarr2.append(iarr2_orig[k])
       new_jarr2.append(jarr2_orig[k])
       new_karr2.append(karr2_orig[k])

#for xx in iarr:
#    xx = int(xx)
#    if xx in atoms_image:
#        ires_orig, atomname = atoms_image[xx]
#        xx = atoms_orig[(ires_orig, atomname)]
#    new_iarr.append(xx)
#for xx in jarr:
#    xx = int(xx)
#    if xx in atoms_image:
#        ires_orig, atomname = atoms_image[xx]
#        xx = atoms_orig[(ires_orig, atomname)]
#    new_jarr.append(xx)
## angle
#for xx in iarr2:
#    xx = int(xx)
#    if xx in atoms_image:
#        ires_orig, atomname = atoms_image[xx]
#        xx = atoms_orig[(ires_orig, atomname)]
#    new_iarr2.append(xx)
#for xx in jarr2:
#    xx = int(xx)
#    if xx in atoms_image:
#        ires_orig, atomname = atoms_image[xx]
#        xx = atoms_orig[(ires_orig, atomname)]
#    new_jarr2.append(xx)
#for xx in karr2:
#    xx = int(xx)
#    if xx in atoms_image:
#        ires_orig, atomname = atoms_image[xx]
#        xx = atoms_orig[(ires_orig, atomname)]
#    new_karr2.append(xx)
#######################################################################

#
# write output psf file
#
#######################################################################
f = open(output, 'w')
f.writelines(lines['TITLE'])
#f.writelines(lines)
f.writelines(lines['ATOM'])
# combine between last original array and image array
iarr_orig = new_iarr
jarr_orig = new_jarr 
#print(len(new_iarr))
#print(len(new_jarr))
iarr2_orig = new_iarr2
jarr2_orig = new_jarr2
karr2_orig = new_karr2
# remove duplicate bond
remove_duplicate = set()
for x in range(len(iarr_orig)):
    i = int(iarr_orig[x])
    j = int(jarr_orig[x])
    remove_duplicate.add((i,j))
nbond_image = len(remove_duplicate)
# remove duplicate angle
#remove_duplicate2 = set()
#for x in range(len(iarr2_orig)):
#    i = int(iarr2_orig[x])
#    j = int(jarr2_orig[x])
#    k = int(karr2_orig[x])
#    remove_duplicate2.add((i,j,k))
#nangl_image = len(remove_duplicate2)

# write nbond header
remove_duplicate = set()
xbond = 0
nbond_orig = len(iarr_orig)
nbond_orig = len(iarr2_orig)
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
    #else:
    #   i = int(iarr_orig[x])
    #   j = int(jarr_orig[x])
    #   if (i,j) in remove_duplicate: continue
    #   if (j,i) in remove_duplicate: continue
    #   f.write('%10d%10d' % (i,j))
    #   remove_duplicate.add((i,j))
    #   xbond += 1
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
        
# write angle information
xangl = 0
nangl_image = 0
lines_angle = []
for x in range(len(iarr2_orig)):
    if xangl > 2:
        #f.write('\n')
        lines_angle.append('\n')
        xangl = 0
#    else:
#       i = int(iarr2_orig[x])
#       j = int(jarr2_orig[x])
#       k = int(karr2_orig[x])
#       if (i,j,k) in remove_duplicate2: continue
#       if (k,j,i) in remove_duplicate2: continue
#       #f.write('%10d%10d%10d' % (i,j,k))
#       remove_duplicate2.add((i,j,k))
#       angle_lines.append('%10d%10d%10d' % (i,j,k))
#       xangl += 1
#       all_angle += 1
    i = int(iarr2_orig[x])
    j = int(jarr2_orig[x])
    k = int(karr2_orig[x])
    if (i,j,k) in remove_duplicate2: continue
    if (k,j,i) in remove_duplicate2: continue
    remove_duplicate2.add((i,j,k))
    lines_angle.append('%10d%10d%10d' % (i,j,k))
    nangl_image += 1
    xangl += 1

if xangl != 0:
    lines_angle.append('\n');

# write nangle header
#for line in lines_angle:
#    if not line.strip(): continue
#    print(line.rstrip())
#nangl_orig += nangl_image
nangl_image += norig_angle
#f.write('\n%9d !NTHETA: angles\n' % nangl_orig)
f.write('\n%9d !NTHETA: angles\n' % nangl_image)
f.writelines(lines['ANGL'])
f.writelines(lines_angle)
#print(lines_angle)
#print(lines['OTHER'][:10])
#input()

f.writelines(lines['OTHER'][2:])
f.close()
#######################################################################
