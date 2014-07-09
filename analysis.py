import os
import sys
import math
import heapq
import itertools
import timeit
import itertools
import pickle
from util import *


def get_dist_matrix_and_IPs(pdb, cutoff):
	N = len(pdb)
	distance_matrix = [[None]*(N-i-1) for i in xrange(N)]
	lib_ips = set()
	first = True
	print '\tOut of %d residues finished:' % N
	for i in xrange(N):
		if i % 5 == 0:
			if first:
				sys.stdout.write('\t')
				first = False
			print i,
			sys.stdout.flush()
		for j in xrange(i+1, N):
			result = pdb[i].filteredDistancesTo(pdb[j], cutoff)
			if result is not None:
				lib_ips.add((i,j))
				distance_matrix[i][j-i-1] = result
	return (distance_matrix, lib_ips)


def preprocess(in_dir, out_dir, ip_dist_cutoff, restart, verbose):
	lib_pdbs = [os.path.join(in_dir, pdb) for pdb in os.listdir(in_dir)]

	for (i, lib_pdb_path) in enumerate(lib_pdbs):

		out_path = os.path.join(out_dir, os.path.split(lib_pdb_path)[1].split('.')[0] + '.pickle') 
                if os.path.exists(out_path) and not restart:
                        if(verbose): print 'Skipping ' + out_path + '. Already in directory.'
                        continue
 
		if verbose: print "\nProcessing library .pdb: " + lib_pdb_path + "\t (" + str(i) + " out of " + str(len(lib_pdbs)) + ")"
		lib_pdb = get_pdb_info(lib_pdb_path)	
		
		if verbose: print '\tComputing inter-atomic distances and finding interacting pairs...'
		(lib_distance_matrix, lib_ips) = get_dist_matrix_and_IPs(lib_pdb, ip_dist_cutoff)

		if verbose: print '\tPreparing directory folders...'
		if not os.path.exists(out_dir): os.makedirs(out_dir)

		if verbose: print '\tPrinting PICKLE file...'
		with open(out_path,'w') as outfile:
			pickle.dump((lib_pdb,lib_ips,lib_distance_matrix), outfile)

def get_dist_matrix_and_IPs_template(pdb, residues, cutoff):
	N = len(pdb)
	distance_matrix = [[None]*(N-i-1) for i in xrange(N)]
	ips = set()
	first = True
	print '\tOut of %d residues finished:' % N
	for i in xrange(N):
		if i % 5 == 0:
			if first:
				sys.stdout.write('\t')
				first = False
			print i,
			sys.stdout.flush()
		for j in xrange(i+1, N):
			if ((i in residues) ^ (j in residues)):
				result = pdb[i].filteredDistancesTo(pdb[j], cutoff)
				if result is not None:
					ips.add((i,j))
					distance_matrix[i][j-i-1] = result
	return (distance_matrix, ips)


def statium(in_res, in_pdb, in_dir, out_dir, ip_cutoff_dist, match_cutoff_dists, counts, verbose):
	
	if verbose: print '\nPreparing directory folders...'
	if not os.path.exists(out_dir): os.makedirs(out_dir)

	if verbose: print 'Starting STATUM analysis...' 
	tic = timeit.default_timer()
	sidechain(in_res, in_pdb, in_dir, out_dir, ip_cutoff_dist, match_cutoff_dists, counts, verbose)
	toc = timeit.default_timer()
	if verbose: print 'Done in ' + str((tic-toc)/60) + 'minutes! Output in: ' + out_dir


def sidechain(in_res, in_pdb, in_preprocess_dir, out_dir, ip_dist_cutoff, match_dist_cutoffs, print_counts, verbose):
	
	if verbose: print 'Extracting residue position from ' + in_res + '...'
	res_lines = filelines2list(in_res)
	residues = [int(line.strip()) for line in res_lines]
	
	if verbose: print 'Extracting information from ' + in_pdb + '...'
	pdbI = get_pdb_info(in_pdb)
	pdbSize = len(pdbI)	

	if verbose: print 'Computing inter-atomic distances and finding interacting pairs...\n'
	(distance_matrix, use_indices) = get_dist_matrix_and_IPs_template(pdbI, residues, ip_dist_cutoff)
	num_ips = len(use_indices)
	
	if verbose: print 'Storing distance information for each interacting pair...'
	distances = dict() 
	for (p1,p2) in use_indices:
		if pdbI[p2].string_name == 'GLY':
			pdbI[p2].correct()
			distance_matrix[p1][p2-p1-1] = pdbI[p1].distancesTo(pdbI[p2])
		if pdbI[p2].stubIntact:
			p2_base_chain = [pdbI[p2].atom_dict['CA'], pdbI[p2].atom_dict['CB']]
			distances[(p1,p2)] = filter_sc_dists(pdbI[p1].atoms, p2_base_chain, distance_matrix[p1][p2-p1-1], 'forward') 
		else:
			print 'Position %d does not have a valid stub' % p2
			sys.exit(1)

	#stores total number of each residue identity across IP
	totals = [0 for i in range(20)] 
	#stores total number of 'matching sidechain'
	counts = [[0 for j in range(20)] for i in range(num_ips)]  
	
	lib_pdb_paths = [os.path.join(in_preprocess_dir, pdb) for pdb in os.listdir(in_preprocess_dir)]
	for (i, lib_pdb_path) in enumerate(lib_pdb_paths):	
		if verbose: print 'Processing ' + lib_pdb_path + ' (' + str(i) + ' out of ' + str(len(lib_pdb_paths)) + ')'
		with open(lib_pdb_path, 'r') as infile:
			(lib_pdb,lib_ips,lib_distance_matrix) = pickle.load(infile)

		for (lib_pos1, lib_pos2) in lib_ips:
			if lib_pos2 - lib_pos1 <= 4: continue
			
			(lib_AA1, lib_AA2) = (lib_pdb[lib_pos1].int_name, lib_pdb[lib_pos2].int_name)
			if lib_AA1 < 0 or lib_AA2 < 0 or lib_AA1 > 19 or lib_AA2 > 19: continue
			
			totals[lib_AA1] += 1
			totals[lib_AA2] += 1
			
			for (j, (pos1, pos2)) in enumerate(use_indices):
				d = distances[(pos1,pos2)]
				if d is None: continue
				AA1 = pdbI[pos1].int_name 
				if lib_pdb[lib_pos2].stubIntact and lib_AA1==AA1:
					p2_base_chain = [lib_pdb[lib_pos2].atom_dict['CA'], lib_pdb[lib_pos2].atom_dict['CB']] 
					lib_dist = filter_sc_dists(lib_pdb[lib_pos1].atoms, p2_base_chain, lib_distance_matrix[lib_pos1][lib_pos2-lib_pos1-1], "forward")
					if matching_sidechain_pair(distances[(pos1,pos2)], lib_dist, match_dist_cutoffs[pdbI[pos1].char_name]):
						counts[j][lib_AA2] += 1
						
				if lib_pdb[lib_pos1].stubIntact and lib_AA2==AA1:
					p1_base_chain = [lib_pdb[lib_pos1].atom_dict['CA'], lib_pdb[lib_pos1].atom_dict['CB']] 
					lib_dist = filter_sc_dists(p1_base_chain, lib_pdb[lib_pos2].atoms, lib_distance_matrix[lib_pos1][lib_pos2-lib_pos1-1], "backward") 
					if matching_sidechain_pair(distances[(pos1,pos2)], lib_dist, match_dist_cutoffs[pdbI[pos1].char_name]):
						counts[j][lib_AA1] += 1
	
	if(verbose): print('Finished processing library .pdb files.')

	if print_counts:
		if verbose: print 'Writing raw library counts to files...'

		counts_dir = out_dir[-1] + '_counts' if out_dir[-1] == '/' else out_dir + '_counts'
		if not os.path.exists(counts_dir):
			os.mkdir(counts_dir)

		total_counts_file = open(os.path.join(counts_dir, 'lib_ip_residue_totals.txt'), 'w')
		for i in range(20): total_counts_file.write(AAint2char(i) + '\t' + str(totals[i]) + '\n')
		total_counts_file.close()

		for (i, pair) in enumerate(use_indices):
			counts_file = open(os.path.join(counts_dir, str(pair[0] + 1) + '_' + str(pair[1] + 1) + '_counts.txt'), 'w')
			for j in range(20): counts_file.write(AAint2char(j) + '\t' + str(counts[i][j]) + '\n')
			counts_file.close()
	
	if(verbose): print("Computing probabilities from counts...")
	determine_probs(use_indices, totals, counts, out_dir, verbose)


def determine_probs(use_indices, totals, counts, out_dir, verbose):
	
	#create the output directory
	if not os.path.exists(out_dir):
		os.mkdir(out_dir)
	
	#total number of residues across all library interacting pairs
	lib_sum = float(sum(totals))
	#frequency of each residue in total library
	lib_total_probs = [x/lib_sum for x in totals]
	
	for (i, pair) in enumerate(use_indices):
		total = float(sum(counts[i]))
		if total > 99:
			path = os.path.join(out_dir, str(pair[0] + 1) + '_' + str(pair[1] + 1) + '_probs.txt')
			with open(path, 'w') as prob_file:
				#probability of each residue occuring at *that* IP (i.e. / by total res's at the IP)
				AA_probs = [(x/total if x != 0 else 1/total) for x in counts[i]]
				for j in range(20):
					e = -1.0 * math.log(AA_probs[j] / lib_total_probs[j])
					prob_file.write(AAint2char(j) + '\t' + str(e) + '\n')
				prob_file.close()
	if(verbose): print("Finished calculating probabilities. Written to: " + out_dir + '_probs')

	
def matching_sidechain_pair(dists1, dists2, cutoff):

	sd = 0.0
	count = 0

	for atom_pair1 in dists1:
		for atom_pair2 in dists2:
			if atom_pair1[0].sameName(atom_pair2[0]) and atom_pair1[1].sameName(atom_pair2[1]):
				sd += (dists1[atom_pair1]-dists2[atom_pair2])**2
				count += 1

	return math.sqrt(sd / count) < cutoff


#returns a more sparse distance matrix, filled on at IP positions
#	uses the library pdb_info
def get_distance_matrix_ip(pdb, ips):
	
	N = len(pdb)
	distance_matrix = [[None]*N for j in range(N)]
	
	for (i,j) in ips:
		try:
			distance_matrix[i][j] = pdb[i].distancesTo(pdb[j])
		except:
			print len(pdb)
			print pdb
			print i,j,N
			
	return distance_matrix

	
#OLD VERSION: select_sidechain_distances
#New version: returns a dictionary, subset of the dictionary at
#		location in distances matrix defined by an interacting pair (R1, R2)
#	 	a subset defined by (R1's atoms,CA or CB):distance
def filter_sc_dists(atomsA, atomsB, interatomic_distances, direction):
	out = dict()
	if direction == 'forward':
		pairs = list(itertools.product(atomsA, atomsB))
	else:
		pairs = list(itertools.product(atomsB, atomsA))
	return {pair: interatomic_distances[pair] for pair in pairs}
	

def check_cutoff(residue_pairs, cutoff):
	return any([dist < cutoff for dist in residue_pairs.values()])


def generate_random_distribution (in_res, in_probs_dir, num_seqs=100000):
	
	sequence_length = len(filelines2list(in_res))

	print 'Calculating distribution of energies for %d random sequences...', num_seqs
	energies = [calc_seq_energy(in_res, in_probs_dir, generate_random_seq(sequence_length)) for _ in num_seqs]
		
	energies.sort()
	avg = mean(energies)
	sd = std(energies)
		
	return (sequence_length, energies, avg, sd)


def calc_seq_energy (in_res_path, probs_dir, seq):

	#Handle two cases: e.g. AAAX,LLL and LXAAM
	seq = ''.join(seq.split(','))
	
	#loading in probability into all_probs
	prob_files = os.listdir(probs_dir)
	all_probs = dict()
	
	for f in prob_files:
		file_path = os.path.join(probs_dir, f)
		lines = filelines2deeplist(file_path)
		probs = [float(prob.split('\t')[0]) for prob in lines]
		ip_res1 = int(f.split('_')[0]) 
		ip_res2 = int(f.split('_')[1])
		all_probs[(ip_res1,ip_res2)] = probs

	energy = 0
	lines = filelines2list(in_res_path)
	residues = [int(line.strip()) for line in lines]

	for i, residue in enumerate(residues):
		filtered = {ip: probs for ip, probs in all_probs.items() if ip[1] == residue}
		AA = seq[i]
		if AA == 'X' or AA == 'G' or filtered == {}: continue
		for ip, probs in filtered.items():
			energy += probs[AAchar2int(AA)]
			
	return energy
	
#Note: need res file, because not all residues 
def calc_top_seqs(in_res_path, probs_dir, num_sequences, outfile):
	 
	#read back from .res file where ligand residues start
	lines = filelines2list(in_res_path)
	residues = [int(line.strip()) for line in lines]
	
   	#loading in probability into all_probs
	prob_files = os.listdir(probs_dir)
	all_probs = OrderedDict()
	
	for f in prob_files:
		file_path = os.path.join(probs_dir, f)
		lines = filelines2deeplist(file_path)
		probs = [float(prob.split('\t')[0]) for prob in lines]
		ip_res1 = int(f.split('_')[0]) 
		ip_res2 = int(f.split('_')[1])
		all_probs[(ip_res1,ip_res2)] = probs

	#the following now fills ordered_probs
	#[[sorted list of AA probs for a residue position: (0.5, 'A'), (0.3, 'C'),...], [like before for residue pos 2], etc...]
	ordered_probs = [] 
	AAs = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
	for residue in residues:
		probs = [prob for ip, prob in all_probs.items() if residue in ip]
		probs_sum = map(sum, zip(*probs))
		sorted_probs = sorted(zip(probs_sum, AAs), key=lambda pair: pair[0])
		ordered_probs.append(sorted_probs)
	
	#enumerate all the balls in urns possibilities
	#note that to maintain the max-heap property in Python's heapq implementation (which only does min-heap), we multiple by -1 before adding to heap
	num_urns = len(residues)
	heap = []
	
	seq = ''
	energy = 0
	for i, c in enumerate(residues):
		if ordered_probs[i]:
			aa = ordered_probs[i][0][1]
			seq += aa
			energy += (0.0 if aa == 'G' else ordered_probs[i][0][0])
		else:
			seq += 'X'
	heapq.heappush(heap, (energy, seq))
	
	max_num_balls = 0
	total = 0
	while total < num_sequences:
		max_num_balls += 1
		total += nCr((max_num_balls+num_urns-1), (max_num_balls))
	
	for num_balls in range(max_num_balls+1)[1:]:
		combo_elements = range(num_balls+num_urns-1)
		combos = list(itertools.combinations(combo_elements, (num_urns-1)))
				
		for combo in combos:
			#print combo
			urn_counts = [0]*num_urns
			for i, position in enumerate(combo):
				if(i == 0):
					urn_counts[i] = position
				elif(i == len(combo)-1):
					urn_counts[i] = position - combo[i-1] - 1
					urn_counts[i+1] = len(combo_elements)-1-position
				else:
					urn_counts[i] = position - combo[i-1] - 1
			
				seq = ''
				energy = 0
				for i, c in enumerate(urn_counts):
					if(ordered_probs[i]):
						aa = ordered_probs[i][c][1]
						seq += aa
						energy += (0.0 if aa == 'G' else ordered_probs[i][c][0])
					else:
						seq += 'X'
			
			if(seq in [i[1] for i in heap]):
				continue
				
			if(len(heap) < num_sequences):
				heapq.heappush(heap, (-1*energy, seq))
			elif(energy < heap[0][0]*-1):
				heapq.heappushpop(heap, (-1*energy, seq))
	
	heap = sorted(heap, reverse=True)
	out = [(seq, energy*-1) for energy, seq in heap]
	return out


