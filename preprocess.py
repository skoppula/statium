from util import AAint2char
from analysis import get_pdb_info
from analysis import filter_sc_dists
import os

def get_dist_matrix_and_IPs_lib(pdb, cutoff):
	N = len(pdb)
	distance_matrix = [[None]*(N-i-1) for i in xrange(N)]
	ips = set()
	first = True
	for i in xrange(N):
	    for j in xrange(i+1, N):
	        result = pdb[i].fastFilteredDistancesTo(pdb[j], cutoff)
	        distance_matrix[i][j-i-1] = result
		if result is not None:
		    ips.add((i,j))

	return (distance_matrix, ips)
#AA at end
#remove []
#round
#ensure order is same
#too many ips?
def preprocess(in_dir, out_dir, ip_dist_cutoff, backbone=False, filter_sidechains=False, correct=False, verbose=True):
    pdb_paths = os.listdir(in_dir)
    num_paths = str(len(pdb_paths))

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    files = list()
    for i in range(20):
        aa = AAint2char(i)
        files.append(open(os.path.join(out_dir, aa + '.dists'), 'w'))


    for i, pdb_name in enumerate(pdb_paths):
        pdb_path = os.path.join(in_dir, pdb_name)
	if verbose: print 'Reading ' + pdb_path + '...' + str(i) + '/' + num_paths
        if filter_sidechains: pdb = get_pdb_info(pdb_path, filter_sidechains = True)
        else: pdb = get_pdb_info(pdb_path)
        if correct:
            for res in pdb:
                if res.string_name == 'GLY':
                    res.correct()
                if not backbone:
                    res.strip_backbone()
        elif not backbone:
            for res in pdb:
                res.strip_backbone()
	pdbSize = len(pdb)	

	(distance_matrix, use_indices) = get_dist_matrix_and_IPs_lib(pdb, ip_dist_cutoff)

	for (p1,p2) in use_indices:
            aa1 = pdb[p1].int_name
            aa2 = pdb[p2].int_name
            if pdb[p2].stubIntact:
		p2_base_chain = ['CA', 'CB'] 
                dists = filter_sc_dists(pdb[p1].atom_names, p2_base_chain, distance_matrix[p1][p2-p1-1], True)
                files[aa1].write(str(dists.keys()) + ';' + str(dists.values()) + '\n')

            if pdb[p1].stubIntact:
		p1_base_chain = ['CA', 'CB'] 
		lib_dist = filter_sc_dists(p1_base_chain, pdb[p2].atom_names, distance_matrix[p1][p2-p1-1], False) 
                files[aa2].write(str(dists.keys()) + ';' + str(dists.values()) + '\n')
    for i in range(20):
        aa = AAint2char(i)
        files[i].close()

