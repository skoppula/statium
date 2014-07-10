import os
from util import *

in_dir = 'data/culled_90'
lib_pdbs = [os.path.join(in_dir, pdb) for pdb in os.listdir(in_dir)][0:100]
tetra_errors = list()
scwrl_errors = list()
for (i, lib_pdb_path) in enumerate(lib_pdbs):
	print lib_pdb_path
	print 'starting tetra analysis'
	lib_pdb = get_pdb_info(lib_pdb_path)
	for aa in lib_pdb:
		if aa.string_name == 'ALA':
			try:
				true_x = aa.atom_dict['CB'].x
			except KeyError:
				continue
			true_y = aa.atom_dict['CB'].y
			true_z = aa.atom_dict['CB'].z
			del aa.atom_dict['CB']
			for i, atom in enumerate(aa.atoms):
				if atom.name is 'CB':
					aa.atoms.pop(i)
					aa.atom_names.pop(i)
					break
			aa.correct()
			error = aa.atom_dict['CB'].distanceTo(Atom('CB',true_x,true_y,true_z,0))
			tetra_errors.append(error)
	print 'starting scwrl analysis'
			
	seq_path = 'scwrl-seq.temp'
	with open(seq_path, 'w') as f:
		for aa in lib_pdb:
			out = 'a' if aa.string_name == 'GLY' else aa.char_name
			f.write(out)

	out_pdb_temp = 'out_temp.pdb'
	os.system('/home/potapov/bin/scwrl4/Scwrl4 -i ' + lib_pdb_path + ' -o ' + out_pdb_temp + ' -s ' + seq_path + ' -h -t >/dev/null 2>&1')

	scwrl_pdb = get_pdb_info(out_pdb_temp)
	lib_pdb = get_pdb_info(lib_pdb_path)

	for i, aa in enumerate(lib_pdb):
		if aa.string_name == 'ALA':
			try:
				error = aa.atom_dict['CB'].distanceTo(scwrl_pdb[i].atom_dict['CB'])
			except:
				continue
			scwrl_errors.append(error)

	print 'tetra_errors', tetra_errors[-10:]
	print 'scwrl_errors', scwrl_errors[-10:]


import matplotlib.pyplot as plt
import numpy as np

hist, bins = np.histogram(tetra_errors, bins=50)
width = 0.7 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2
plt.title('tetra error')
plt.xlabel('error (Angstroms)')
plt.ylabel('frequency')
plt.axis([0,5,0,2000])
plt.bar(center, hist, align='center', width=width)
plt.show()


hist, bins = np.histogram(scwrl_errors, bins=50)
width = 0.7 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2
plt.title('scwrl error')
plt.xlabel('error (Angstroms)')
plt.ylabel('frequency')
plt.axis([0,5,0,2000])
plt.bar(center, hist, align='center', width=width)
plt.show()
#print tetra_errors
#print scwrl_errors
		
