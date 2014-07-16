from util import *
import os
lib_pdb_path = '/home/skoppula/statium/testing/1mhp_AHL_renum.pdb'
lib_pdb = get_pdb_info(lib_pdb_path)
seq_path = 'scwrl-seq.temp'
with open(seq_path, 'w') as f:
	for aa in lib_pdb:
		out = 'a' if aa.string_name == 'GLY' else aa.char_name
		f.write(out)

out_pdb_temp = 'out_temp.pdb'
os.system('/home/potapov/bin/scwrl4/Scwrl4 -i ' + lib_pdb_path + ' -o ' + out_pdb_temp + ' -s ' + seq_path + ' -h -t >/dev/null 2>&1')

