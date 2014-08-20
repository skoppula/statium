import os
import sys
import math
import heapq
import itertools
import timeit
import itertools
import pickle
from util import *
from collections import OrderedDict


def statium(in_res, in_pdb, in_dir, in_ip, out, ip_cutoff_dist, match_cutoff_dists, backbone, filter_sidechain, counts, verbose):
	
	if verbose: print 'Starting STATUM analysis...' 
	tic = timeit.default_timer()
	sidechain(in_res, in_pdb, in_dir, in_ip, out, ip_cutoff_dist, match_cutoff_dists, backbone, filter_sidechain, counts, verbose)
	toc = timeit.default_timer()
	if verbose: print 'Done in ' + str((toc-tic)/60) + ' minutes! Output in: ' + out


def sidechain(in_res, in_pdb, in_pdb_dir, in_ip_dir, out, ip_dist_cutoff, match_dist_cutoffs, backbone, filter_sidechain, print_counts, verbose):
	
	if verbose: print 'Extracting residue position from ' + in_res + '...'
	res_lines = filelines2list(in_res)
	residues = [int(line.strip())-1 for line in res_lines]
	
	if verbose: print 'Extracting information from ' + in_pdb + '...'
        if filter_sidechain: pdbI = get_pdb_info(in_pdb, filter_sidechains = True)
        else: pdbI = get_pdb_info(in_pdb)
	for res in pdbI:
		if res.string_name == 'GLY':
			res.correct()
                if not backbone:
		    res.strip_backbone()
	pdbSize = len(pdbI)	

	if verbose: print 'Computing inter-atomic distances and finding interacting pairs...\n'
	(distance_matrix, use_indices) = get_dist_matrix_and_IPs_peptide(pdbI, residues, ip_dist_cutoff)
        if verbose: print use_indices
	num_ips = len(use_indices)
	
	if verbose: print 'Storing distance information for each interacting pair...'
	distances = dict() 
	for (p1,p2) in use_indices:
			if pdbI[p2].stubIntact:
					p2_base_chain = ['CA','CB']
					distances[(p1,p2)] = filter_sc_dists(pdbI[p1].atom_names, p2_base_chain, distance_matrix[p1][p2-p1-1], 'forward') 
			else:
					print 'Position %d does not have a valid stub' % p2
					sys.exit(1)

	#stores total number of each residue identity across IP
	totals = [0 for i in xrange(20)] 
	#stores total number of 'matching sidechain'
	counts = [[0 for j in xrange(20)] for i in xrange(num_ips)]  
	
	lib_pdb_paths = [os.path.join(in_pdb_dir, pdb) for pdb in os.listdir(in_pdb_dir)]
	for (i, lib_pdb_path) in enumerate(lib_pdb_paths):	
		if verbose: print 'Processing ' + lib_pdb_path + ' (' + str(i) + ' of ' + str(len(lib_pdb_paths)) + ')'
		if in_ip_dir:
			lib_pdb = get_pdb_info(lib_pdb_path)
			ip_path = os.path.join(in_ip_dir, os.path.split(lib_pdb_path)[1].split('.')[0] + '.ip')
			lib_ips = [(int(pair[0]),int(pair[1])) for pair in filelines2deeplist(ip_path) if pair != []]
			lib_distance_matrix = get_lib_dist_matrix(lib_pdb, lib_ips)
			if not lib_distance_matrix:
				print lib_pdb_path
				continue
		else:
			with open(lib_pdb_path, 'r') as infile:
				(lib_pdb,lib_ips,lib_distance_matrix) = pickle.load(infile)

		for (lib_pos1, lib_pos2) in lib_ips:
			if abs(lib_pos2 - lib_pos1) <= 4: continue
			
			(lib_AA1, lib_AA2) = (lib_pdb[lib_pos1].int_name, lib_pdb[lib_pos2].int_name)
			if lib_AA1 < 0 or lib_AA2 < 0 or lib_AA1 > 19 or lib_AA2 > 19: continue
			
			totals[lib_AA1] += 1
			totals[lib_AA2] += 1
			
			for (j, (pos1, pos2)) in enumerate(use_indices):
				AA1 = pdbI[pos1].int_name 
                                if lib_pdb[lib_pos2].stubIntact:
                                    if lib_AA1==AA1:
					p2_base_chain = ['CA', 'CB'] 
                                        lib_dist = filter_sc_dists(lib_pdb[lib_pos1].atom_names, p2_base_chain, lib_distance_matrix[lib_pos1][lib_pos2-lib_pos1-1], True)
					if matching_sidechain_pair(distances[(pos1,pos2)], lib_dist, match_dist_cutoffs[pdbI[pos1].char_name]):
						counts[j][lib_AA2] += 1
						
                                if lib_pdb[lib_pos1].stubIntact:
                                    if lib_AA2==AA1:
					p1_base_chain = ['CA', 'CB'] 
					lib_dist = filter_sc_dists(p1_base_chain, lib_pdb[lib_pos2].atom_names, lib_distance_matrix[lib_pos1][lib_pos2-lib_pos1-1], False) 
					if matching_sidechain_pair(distances[(pos1,pos2)], lib_dist, match_dist_cutoffs[pdbI[pos1].char_name]):
						counts[j][lib_AA1] += 1
	
	if(verbose): print('Finished processing library .pdb files.')

	if print_counts:
		if verbose: print 'Writing raw library counts to files...'

		counts_dir = out + '_counts'
		if not os.path.exists(counts_dir):
			os.mkdir(counts_dir)

		total_counts_file = open(os.path.join(counts_dir, 'lib_ip_residue_totals.txt'), 'w')
		for i in xrange(20): total_counts_file.write(AAint2char(i) + '\t' + str(totals[i]) + '\n')
		total_counts_file.close()

		for (i, pair) in enumerate(use_indices):
			counts_file = open(os.path.join(counts_dir, str(pair[0] + 1) + '_' + str(pair[1] + 1) + '_counts.txt'), 'w')
			for j in xrange(20): counts_file.write(AAint2char(j) + '\t' + str(counts[i][j]) + '\n')
			counts_file.close()
	
	if verbose: print("Computing probabilities from counts...")
	probs = determine_probs(use_indices, totals, counts)

	if len(use_indices) != len(probs):
		print 'Number of IPs is not consistent'
		print use_indices
		print probs
		sys.exit(1)

	write_output(use_indices, probs, out)	
	if(verbose): print("Finished calculating probabilities. Written to: " + out + '_probs')


def get_dist_matrix_and_IPs_peptide(pdb, residues, cutoff):
	N = len(pdb)
	distance_matrix = [[None]*(N-i-1) for i in xrange(N)]
	ips = set()
	first = True
	print '\tOut of %d residues finished:' % N
	for i in xrange(N):
		for j in xrange(i+1, N):
			if (i in residues) ^ (j in residues):
				result = pdb[i].fastFilteredDistancesTo(pdb[j], cutoff)
				distance_matrix[i][j-i-1] = result
				if result is not None:
					if i in residues:
						ips.add((j,i))
					else:
						ips.add((i,j))

	return (distance_matrix, ips)


def get_lib_dist_matrix(pdb, ips):
	N = len(pdb)
	matrix = [[None]*(N-i-1) for i in xrange(N)]

	for (i,j) in ips:
		try:
			matrix[i][j-i-1] = pdb[i].fastDistancesTo(pdb[j])
		except IndexError:
			print 'Invalid indices'
			print i,j
			print len(matrix), len(matrix[0])
			print len(pdb)
			return None
	return matrix

def write_output(use_indices, probs, out):
	AAs = [AAint2char(i) for i in xrange(20)]
	with open(out, 'w') as f:
		s = 'IPs\t' + '\t'.join(AAs)
		f.write(s + '\n')
		for l, (i,j) in enumerate(use_indices):
			s = str(i) + '-' + str(j)
			s2 = '\t'.join([str(probs[l][aa]) for aa in AAs])
			f.write(s + '\t' + s2 + '\n')

def read_output(in_path):
	lines = filelines2deeplist(in_path, skipComments=True, skipEmptyLines=True)
	AAs = [AAint2char(i) for i in xrange(20)]
	probs = OrderedDict()
	for ip_line in lines[1:]:
		ip = tuple(ip_line[0].split('-'))
		use_indices.append(ip)
		aa_probs = dict()
		for i, prob in enumerate(ip_line[1:]):
			aa_probs[AAs[i]] = prob
		probs[ip] = aa_probs

	return probs

def determine_probs(use_indices, totals, counts):
	
	#total number of residues across all library interacting pairs
	lib_sum = float(sum(totals))

	#frequency of each residue in total library
	lib_total_probs = [x/lib_sum for x in totals]
	out = list()
	
	for (i, pair) in enumerate(use_indices):
		total = float(sum(counts[i]))
		if total > 99:
			probs = dict()
			#probability of each residue occuring at *that* IP (i.e. / by total res's at the IP)
			AA_probs = [(x/total if x != 0 else 1/total) for x in counts[i]]
			for j in xrange(20):
				e = -1.0 * math.log(AA_probs[j] / lib_total_probs[j])
				probs[AAint2char(j)] = e
			out.append(probs)
	return out

#Assumes fast & distances stored as dict (atom1_name, atom2_name): distance val
def matching_sidechain_pair(dists1, dists2, cutoff):
  
    sd = 0.0
    count = 0.0

    for pair1, dist1 in dists1.iteritems():
        for pair2, dist2 in dists2.iteritems():
            if pair1 == pair2:
	        sd += ((dist1 - dist2) ** 2)
	        count += 1.0

    if math.sqrt(sd / count) < cutoff: return True
    else: return False

	

#returns a more sparse distance matrix, filled on at IP positions
#	uses the library pdb_info
def get_distance_matrix_ip(pdb, ips):
	
	N = len(pdb)
	distance_matrix = [[None]*N for j in xrange(N)]
	
	for (i,j) in ips:
		try:
			distance_matrix[i][j] = pdb[i].distancesTo(pdb[j])
		except:
			print len(pdb)
			print pdb
			print i,j,N
			
	return distance_matrix
	

#If fast, assumes pos1_list is atom_names (not atom objects), and pair_dists is indexed the same way
def filter_sc_dists(pos1_list, pos2_list, pair_dists, forward, output_dict=True, fast=True):
  
	if output_dict:
			out = dict()	
			if forward:
				for atomi in pos1_list:
					for atomj in pos2_list:
							out[(atomi, atomj)] = pair_dists[(atomi, atomj)]
			else:
				for atomj in pos2_list:
					for atomi in pos1_list:
							out[(atomj, atomi)] = pair_dists[(atomi, atomj)]

			return out

	else:
			pairs = list()
			dists = list()

			if forward:
				for atomi in pos1_list:
					for atomj in pos2_list:
							idx = distance_matrix[0].index((atomi, atomj))
							pairs.append((atomi, atomj))
							dists.append(distance_matrix[1][idx])
			else:
				for atomj in pos2_list:
					for atomi in pos1_list:
							idx = distance_matrix[0].index((atomi, atomj))
							pairs.append((atomj, atomi))
							dists.append(distance_matrix[1][idx])

			return [pairs, dists]



def generate_random_distribution (in_res, in_probs_dir, num_seqs=1000):
	
	sequence_length = len(filelines2list(in_res))

	print 'Calculating distribution of energies for %d random sequences...', num_seqs
	energies = list()
	for i in xrange(num_seqs):
		if i % 100 is 0:
			print i,
			sys.stdout.flush()
		energies.append(calc_seq_energy(in_res, in_probs_dir, generate_random_seq(sequence_length)))
		
	energies.sort()
	avg = mean(energies)
	sd = std(energies)
		
	return (sequence_length, energies, avg, sd)


def calc_seq_energy (in_res_path, in_probs, seq):

	#Handle two cases: e.g. AAAX,LLL and LXAAM
	seq = ''.join(seq.split(','))
	
	#loading in probability into all_probs
	all_probs = read_output(in_probs)
	
	energy = 0
	lines = filelines2list(in_res_path)
	residues = [int(line.strip()) for line in lines]

	for i, residue in enumerate(residues):
		filtered = {ip: probs for ip, probs in all_probs.items() if ip[1] == residue}
		AA = seq[i-1]
		if AA == 'X' or AA == 'G' or filtered == {}: continue
		for ip, probs in filtered.items():
			energy += probs[AAchar2int(AA)]
			
	return energy
	
#Note: need res file, because not all residues 
def calc_top_seqs(in_res_path, in_probs, num_sequences):
	 
	#read back from .res file where ligand residues start
	lines = filelines2list(in_res_path)
	residues = [int(line.strip()) for line in lines]
	
   	#loading in probability into all_probs
	prob_files = os.listdir(probs_dir)
	all_probs = read_output(in_probs)
	
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
	
	for num_balls in xrange(max_num_balls+1)[1:]:
		combo_elements = xrange(num_balls+num_urns-1)
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


