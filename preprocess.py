from util import AAint2char
from analysis import get_pdb_info
from analysis import filter_sc_dists
import os

#/home/bartolo/statium_cpp_params/intrapeptide_longrange/
#for complete libraries

#NO-FILTER OPTION IS NOT SUPPORTED
def preprocess(in_dir, out_dir, preprocess_type, ip_dist_cutoff, backbone=False, filter_sidechains=True, correct=False, verbose=True):
	if 'all' == preprocess_type:
		if verbose: print 'Creating all four library types...'
		process_all(in_dir, out_dir, ip_dist_cutoff, backbone, filter_sidechains, correct, verbose)

	elif 'long_fixed_variable' == preprocess_type:
		if verbose: print 'Creating long-range fixed-variable library...'
		long_fixed_variable(in_dir, out_dir, ip_dist_cutoff, backbone, filter_sidechains, correct, verbose)

	elif 'short_fixed_variable' == preprocess_type:
		if verbose: print 'Creating short-range fixed-variable library...'
		short_fixed_variable(in_dir, out_dir, ip_dist_cutoff, backbone, filter_sidechains, correct, verbose)

	elif 'long_variable_variable' == preprocess_type:
		if verbose: print 'Creating short-range variable-variable library...'
		short_variable_variable(in_dir, out_dir, ip_dist_cutoff, backbone, filter_sidechains, correct, verbose)

	elif 'short_variable_variable' == preprocess_type:
		if verbose: print 'Creating short-range variable-variable library...'
		short_variable_variable(in_dir, out_dir, ip_dist_cutoff, backbone, filter_sidechains, correct, verbose)
		
def setup(in_dir, correct, backbone, filter_sidechains):

	pdb_paths = os.listdir(in_dir)
	num_paths = str(len(pdb_paths))
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

	return pdb


def long_fixed_variable(in_dir, out_dir, ip_dist_cutoff, backbone=False, filter_sidechains=True, correct=False, verbose=True):
	out = os.path.join(out_dir,'long_range_fixed-variable')
	if not os.path.exists(out):
		os.mkdir(out)
	for i in range(20):
		aa = AAint2char(i)
		files.append(open(os.path.join(out, aa + '.dists'), 'w'))
	
	pdb = setup(in_dir, correct, backbone, filter_sidechains)
	pdbSize = len(pdb)

	(distance_matrix, use_indices) = get_dist_and_IPs(pdb, ip_dist_cutoff)

	for (p1,p2) in use_indices:
		aa1 = pdb[p1].int_name
		aa2 = pdb[p2].int_name
		if pdb[p2].stubIntact:
			p2_base_chain = ['CA', 'CB'] 
			dists = filter_sc_dists(pdb[p1].atom_names, p2_base_chain, distance_matrix[p1][p2-p1-1], True)
			files[aa1].write("".join([round(val, 2) for val in dists.values()])[1:-1] + ';' + str(aa2) + '\n')

		if pdb[p1].stubIntact:
			p1_base_chain = ['CA', 'CB'] 
			dists = filter_sc_dists(p1_base_chain, pdb[p2].atom_names, distance_matrix[p1][p2-p1-1], False)
			files[aa2].write("".join([round(val, 2) for val in dists.values()])[1:-1] + ';' + str(aa1) + '\n')

	for f in files: f.close()


def short_fixed_variable(in_dir, out_dir, ip_dist_cutoff, backbone=False, filter_sidechains=True, correct=False, verbose=True):
	out = os.path.join(out_dir,'short_range_fixed-variable')
	if not os.path.exists(out):
		os.mkdir(out)
	
	files = list()
	for i in range(20):
		files.append([])
		for j in range(1,5):
			files[i].append([])
			for k in ['rev', 'for']: #why do we need to seperate out the 'rev' and 'for' dists?
				aa = AAint2char(i)
				files[i][j].append(open(os.path.join(out, aa + '_' + j + '_' + k + '.dists'), 'w'))
	
	pdb = setup(in_dir, correct, backbone, filter_sidechains)
	pdbSize = len(pdb)
	(distance_matrix, use_indices) = get_dist_and_IPs_intra_full(pdb, ip_dist_cutoff)

	for (p1,p2) in use_indices:
		aa1 = pdb[p1].int_name
		aa2 = pdb[p2].int_name
		if pdb[p2].stubIntact: #need to check whether both stubs are intact?
			p2_base_chain = ['CA', 'CB'] 
			dists = filter_sc_dists(pdb[p1].atom_names, p2_base_chain, distance_matrix[p1][p2-p1-1], True)
			files[aa1][p2-p1-1][0].write("".join([round(val, 2) for val in dists.values()])[1:-1] + ';' + str(aa2) + '\n')

		if pdb[p1].stubIntact:
			p1_base_chain = ['CA', 'CB'] 
			dists = filter_sc_dists(p1_base_chain, pdb[p2].atom_names, distance_matrix[p1][p2-p1-1], False)
			files[aa2][p2-p1-1][1].write("".join([round(val, 2) for val in dists.values()])[1:-1] + ';' + str(aa1) + '\n')


	for i in range(20):
		for j in range(1,5):
			for k in ['rev', 'for']:
				files[i][j][k].close()
	

def short_variable_variable(in_dir, out_dir, ip_dist_cutoff, backbone=False, filter_sidechains=True, correct=False, verbose=True):
	out = os.path.join(out_dir,'short_range_variable-variable')
	for i in range(1,5):
		files.append(open(os.path.join(out, i + '.dists'), 'w'))

	pdb = setup(in_dir, correct, backbone, filter_sidechains)
	pdbSize = len(pdb)

	(distance_matrix, use_indices) = get_dist_and_IPs_intra_CACBOnly(pdb, ip_dist_cutoff)

	for (p1,p2) in use_indices:
		aa1 = pdb[p1].int_name
		aa2 = pdb[p2].int_name
		if pdb[p2].stubIntact and pdb[p2].stubIntact:
			dists = filter_sc_dists(['CA', 'CB'], ['CA', 'CB'], distance_matrix[p1][p2-p1-1], True)
			files[p2-p1-1].write("".join([round(val, 2) for val in dists.values()])[1:-1] + ';' + str(aa1) + ';' + str(aa2) + '\n')

	for f in files: f.close()

def long_variable_variable(in_dir, out_dir, ip_dist_cutoff, backbone=False, filter_sidechains=True, correct=False, verbose=True):
	outfile = open(os.path.join(out_dir, 'short_range_fixed-variable.dists'), 'w')

	pdb = setup(in_dir, correct, backbone, filter_sidechains)
	pdbSize = len(pdb)

	(distance_matrix, use_indices) = get_dist_and_IPs(pdb, ip_dist_cutoff)

	for (p1,p2) in use_indices:
		aa1 = pdb[p1].int_name
		aa2 = pdb[p2].int_name
		if pdb[p2].stubIntact and pdb[p2].stubIntact:
			dists = filter_sc_dists(['CA', 'CB'], ['CA', 'CB'], distance_matrix[p1][p2-p1-1], True)
			outfile.write(str([round(val, 2) for val in dists.values()])[1:-1] + ';' + str(aa1) + ';' + str(aa2) + '\n')
	
	outfile.close()

def process_all(in_dir, out_dir, ip_dist_cutoff, backbone=False, filter_sidechains=True, correct=False, verbose=True):
	long_fixed_variable()
	short_fixed_variable()
	long_variable_variable()
	long_variable_variable()


def get_dists_and_IPs(pdb, cutoff):
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


def get_dist_and_IPs_intra_CACBOnly(pdb, cutoff):
	N = len(pdb)
	distance_matrix = [[None]*(4) for i in xrange(N)]
	ips = set()
	first = True
	for i in xrange(N):
		for j in xrange(i+1, i+5):
			if j > libN-1: continue
			result = pdb[i].fastDistancesToCACBOnly(pdb[j], cutoff)
			distance_matrix[i][j-i-1] = result
			ips.add((i,j))

	return (distance_matrix, ips)


def get_dist_and_IPs_intra_full(pdb, cutoff):
	N = len(pdb)
	distance_matrix = [[None]*(4) for i in xrange(N)]
	ips = set()
	first = True
	for i in xrange(N):
		for j in xrange(i+1, i+5):
			if j > libN-1: continue
			result = pdb[i].fastDistancesTo(pdb[j], cutoff)
			distance_matrix[i][j-i-1] = result
			ips.add((i,j))

	return (distance_matrix, ips)
