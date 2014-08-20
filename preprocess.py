from util import AAint2char
from analysis import get_pdb_info
from analysis import get_dist_matrix_and_IPs_peptide
from analysis import filter_sc_dists

def preprocess(in_dir, out_dir, ip_dist_cutoff, backbone=False, filter_sidechains=False, verbose):
    pdb_paths = os.listdir(in_dir)
    num_paths = str(len(pdb_paths))

    files = list()
    for i in range(20):
        aa = AAint2char(i)
        files.append(open(aa + '.dists', 'w'))


    for i, pdb_path in enumerate(pdb_paths):
	if verbose: print 'Reading ' + pdb_path + '...' + str(i) + '/' + num_paths
        if filter_sidechain: pdb = get_pdb_info(pdb_path, filter_sidechains = True)
        else: pdb = get_pdb_info(in_pdb)
	for res in pdb:
	    if res.string_name == 'GLY':
	        res.correct()
            if not backbone:
		res.strip_backbone()
	pdbSize = len(pdb)	
        info = get_pdb_info(pdb)

	if verbose: print 'Computing inter-atomic distances ...\n'
	(distance_matrix, use_indices) = get_dist_matrix_and_IPs_peptide(pdbI, residues, ip_dist_cutoff)

	for (p1,p2) in use_indices:
            aa1 = pdb[p1].int_name
            aa2 = pdb[p2].int_name
            if lib_pdb[lib_pos2].stubIntact:
		p2_base_chain = ['CA', 'CB'] 
                dists = filter_sc_dists(lib_pdb[lib_pos1].atom_names, p2_base_chain, lib_distance_matrix[lib_pos1][lib_pos2-lib_pos1-1], True)
                files[aa1].write(str(dists.keys()) + '|' + str(dists.values()))

            if lib_pdb[lib_pos1].stubIntact:
		p1_base_chain = ['CA', 'CB'] 
		lib_dist = filter_sc_dists(p1_base_chain, lib_pdb[lib_pos2].atom_names, lib_distance_matrix[lib_pos1][lib_pos2-lib_pos1-1], False) 
                files[aa2].write(str(dists.keys()) + '|' + str(dists.values()))
    for i in range(20):
        aa = AAint2char(i)
        files[i].close()

