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
			if j-i <= 4: continue
			result = pdb[i].fastFilteredDistancesTo(pdb[j], cutoff)
			distance_matrix[i][j-i-1] = result
		if result is not None:
			ips.add((i,j))

	return (distance_matrix, ips)

def get_dist_matrix_and_IPs_lib_intra(pdb, cutoff):
	N = len(pdb)
	distance_matrix = [[None]*(4) for i in xrange(N)]
	ips = set()
	first = True
	for i in xrange(N):
		for j in xrange(i+1, i+5):
			if j-i <= 4: continue
			result = pdb[i].fastDistancesTo(pdb[j], cutoff)
			distance_matrix[i][j-i-1] = result
			ips.add((i,j))

	return (distance_matrix, ips)

#NO-FILTER OPTION IS NOT SUPPORTED
def preprocess(in_dir, out_dir, preprocess_type, ip_dist_cutoff, backbone=False, filter_sidechains=True, correct=False, verbose=True):

	pdb_paths = os.listdir(in_dir)
	num_paths = str(len(pdb_paths))

	if not os.path.exists(out_dir):
		os.mkdir(out_dir)

	files = list()
	if 'full' in preprocess_type:
		for i in range(20):
			aa = AAint2char(i)
			files.append(open(os.path.join(out_dir, aa + '.dists'), 'w'))
	elif 'cacb' in preprocess_type:
		outfile = open(os.path.join(out_dir, 'cacb.dists'), 'w')
	else:
		for i in range(20):
			aa = AAint2char(i)
			files.append(open(os.path.join(out_dir, aa + '.dists'), 'w'))
		outfile = open(os.path.join(out_dir, 'cacb.dists'), 'w')


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

	if 'long' in preprocess_type:
		(distance_matrix, use_indices) = get_dist_matrix_and_IPs_lib(pdb, ip_dist_cutoff)
	elif 'short' in preprocess_type:
		(distance_matrix, use_indices) = get_dist_matrix_and_IPs_lib_intra(pdb, ip_dist_cutoff)
	else:
		(distance_matrix1, use_indices1) = get_dist_matrix_and_IPs_lib(pdb, ip_dist_cutoff)
		(distance_matrix2, use_indices2) = get_dist_matrix_and_IPs_lib(pdb, ip_dist_cutoff)

	for (p1,p2) in use_indices:
		aa1 = pdb[p1].int_name
		aa2 = pdb[p2].int_name
		if pdb[p2].stubIntact:
			p2_base_chain = ['CA', 'CB'] 
			if 'full' in preprocess_type:
				dists = filter_sc_dists(pdb[p1].atom_names, p2_base_chain, distance_matrix[p1][p2-p1-1], True)
				files[aa1].write(str([round(val, 2) for val in dists.values()])[1:-1] + ';' + str(aa2) + '\n')
			elif 'cacb' in preprocess_type:
				dists = filter_sc_dists(p2_base_chain, p2_base_chain, distance_matrix[p1][p2-p1-1], True)
				outfile.write(str([round(val, 2) for val in dists.values()])[1:-1] + ';' + str(aa2) + '\n')
			else:
				dists1 = filter_sc_dists(pdb[p1].atom_names, p2_base_chain, distance_matrix[p1][p2-p1-1], True)
				dists2 = filter_sc_dists(p2_base_chain, p2_base_chain, distance_matrix[p1][p2-p1-1], True)
				files[aa1].write(str([round(val, 2) for val in dists1.values()])[1:-1] + ';' + str(aa2) + '\n')
				outfile.write(str([round(val, 2) for val in dists2.values()])[1:-1] + ';' + str(aa2) + '\n')

		if pdb[p1].stubIntact:
			p1_base_chain = ['CA', 'CB'] 
			if 'full' in preprocess_type:
				dists = filter_sc_dists(p1_base_chain, pdb[p2].atom_names, distance_matrix[p1][p2-p1-1], False)
				files[aa2].write(str([round(val, 2) for val in dists.values()])[1:-1] + ';' + str(aa1) + '\n')
			elif 'cacb' in preprocess_type:
				dists = filter_sc_dists(p1_base_chain, p1_base_chain, distance_matrix[p1][p2-p1-1], False)
				outfile.write(str([round(val, 2) for val in dists.values()])[1:-1] + ';' + str(aa1) + '\n')
			else:
				dists1 = filter_sc_dists(p1_base_chain, pdb[p2].atom_names, distance_matrix[p1][p2-p1-1], False)
				files[aa2].write(str([round(val, 2) for val in dists1.values()])[1:-1] + ';' + str(aa1) + '\n')
				dists2 = filter_sc_dists(p1_base_chain, p1_base_chain, distance_matrix[p1][p2-p1-1], False)
				outfile.write(str([round(val, 2) for val in dists2.values()])[1:-1] + ';' + str(aa1) + '\n')
	if 'full' in preprocess_type:
		for i in range(20):
			aa = AAint2char(i)
			files[i].close()

	elif 'cacb' in preprocess_type:
		outfile.close()

	else:
		for i in range(20):
			aa = AAint2char(i)
			files[i].close()
		outfile.close()
	

