import os
import sys
import math
import heapq
import itertools
import timeit
import itertools
import pickle
from operator import itemgetter
from util import filelines2list
from util import list2file
from util import isAA
from util import AA2char
from util import AAchar2int
from util import AAint2char
from util import get_sidechain_atoms
from util import filelines2deeplist
from util import binary_placement
from util import mean
from util import std
from util import nCr
from util import read_results
from util import Residue
from util import Atom
from util import get_pdb_info
from reformat import get_orig_seq
from reformat import generate_random_seqs


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
 
		if(verbose): print "\nProcessing library .pdb: " + lib_pdb_path + "\t (" + str(i) + " out of " + str(len(lib_pdbs)) + ")"
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
	residues = [int(line.strip()) - 1 for line in res_lines]
	
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
			distances[(p1,p2)] = filter_sc_dists(pdbI[p1].atoms, p2_base_chain, distance_matrix[p1][p2-p1-1]) 
		else:
			print 'Position %d does not have a valid stub' % p2
			sys.exit(1)

	#stores total number of each residue identity across IP
	totals = [0 for i in range(20)] 
	#stores total number of 'matching sidechain'
	counts = [[0 for j in range(20)] for i in range(num_ips)]  
	
	lib_pdb_paths = [os.path.join(in_preprocess_dir, pdb) for pdb in os.listdir(in_preprocess_dir)]
	for (i, lib_pdb_path) in enumerate(lib_pdb_paths):	
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
					lib_dist = filter_sc_dists(lib_pdb[lib_pos1].atoms, p2_base_chain, lib_distance_matrix[lib_pos1][lib_pos2], "forward")
					if matching_sidechain_pair(distances[(pos1,pos2)], lib_dist, match_dist_cutoffs[pdbI[pos1].char_name]):
						counts[j][lib_AA2] += 1
						
				if lib_pdb[lib_pos1].stubIntact and lib_AA2==AA1:
					p1_base_chain = [lib_pdb[lib_pos1].atom_dict['CA'], lib_pdb[lib_pos1].atom_dict['CB']] 
					lib_dist = filter_sc_dists(p1_base_chain, lib_pdb[lib_pos2].atoms, lib_distance_matrix[lib_pos1][lib_pos2], "backward") 
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
	determine_probs(totals, counts, out_dir, verbose)


def determine_probs(totals, counts, out_dir, verbose):
	
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
			prob_file = open(path, 'w')
			#probability of each residue occuring at *that* IP (i.e. / by total res's at the IP)
			AA_probs = [(x/total if x != 0 else 1/total) for x in counts[i]]
			for j in range(20):
				e = -1.0 * math.log(AA_probs[j] / lib_total_probs[j])
				counts_file.write(AAint2char(j) + '\t' + e + '\n')
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
def filter_sc_dists(atomsA, atomsB, interatomic_distances, direction = 'forward'):
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
	sequences = generate_random_seqs(sequence_length, num_seqs)
	
	if(num_seqs > 5000):
		energies = list(zip(*calc_seq_energy(in_res, in_probs_dir, sequences, verbose=True))[0])
	else:
		energies = list(zip(*calc_seq_energy(in_res, in_probs_dir, sequences, verbose=False))[0])
		
	energies.sort()
	avg = mean(energies)
	sd = std(energies)
		
	return (sequence_length, sequences, energies, avg, sd)


def calc_seq_zscore(mean, std, energy):
	zscore = (energy - mean)/std
	return zscore

#binary search on sorted energies list
def calc_seq_percentile(energies, energy):
	i = binary_placement(energies, energy)
	percentile = i*100/len(energies)
	return percentile

def fix_sequence_line(seq, desired_seq_length, in_pdb_orig=None, warning=False):
	parts = seq.split()
	
	if(len(seq) != desired_seq_length and len(parts) == 1 and warning):
		print('NOTE: IRREGULAR SEQUENCE LENGTH WITHOUT A START POSITION FOR ' + seq)
	
	elif(len(parts) > 1):
		if(in_pdb_orig == None):
			print('Using syntax \'' + seq + '\' needs valid file --IN_PDB_ORIG.')
			return parts[1]
		
		start_seq = int(parts[0])
		start_reference = get_orig_seq(in_pdb_orig)[2]
		
		if(start_seq < start_reference): #the input sequence starts earlier than chain B sequence
				seq = parts[1][(start_reference - start_seq):]
							  
				if(seq == ''):
					print('UNRELATED SEQUENCE INPUT: ' + seq)
						
		elif(start_seq > start_reference): #input sequence starts later than chain B sequence
				seq = 'X'*(start_seq - start_reference) + parts[1]
					
	return seq

def calc_seq_energy (in_res_path, probs_dir, seq, in_pdb_orig=None, verbose=False):
	
	#loading in probability into all_probs
	probs_files = os.listdir(probs_dir)
	all_probs = [[], []] #[[[PROBS FOR IP1], [PROBS FOR IP2], ...], [[IP1], [IP2],...]]
	
	lines = filelines2list(in_res_path)
	residue_positions = [int(line.strip()) - 1 for line in lines]
	
	for file in probs_files:
		file_path = os.path.join(probs_dir, file)
		lines = filelines2deeplist(file_path)

		probs = [float(x[1]) for x in lines]
		all_probs[0].append(probs)
		all_probs[1].append([int(file.split('_')[0]) - 1, int(file.split('_')[1]) - 1])
	
	if(isinstance(seq, str)):
		return sum_energy(residue_positions, all_probs, seq, in_pdb_orig)
	
	elif(isinstance(seq, list)):
		out = []
		for (i,x) in enumerate(seq):
			if(verbose and i%1000 == 1):
				print('Calculating energy of ' + str(i) + 'th random sequence for the distribution...')
			out.append(sum_energy(residue_positions, all_probs, x, in_pdb_orig))
		
		return out
	
def sum_energy(residue_positions, all_probs, seq, in_pdb_orig=None):
	
	#deal with irregularly sized sequences
	seq = fix_sequence_line(seq, len(residue_positions), in_pdb_orig)
	
	energy = 0.0
	for i in range(len(all_probs[0])):
		
		(ip_probs, ip_pos) = (all_probs[0][i], all_probs[1][i])
		pos1 = ip_pos[1] #peptide position on ligand
		
		try:
			if seq[pos1 - residue_positions[0]] == 'X': continue
			if pos1 in residue_positions:
				AA = AAchar2int(seq[pos1 - residue_positions[0]])
		except: continue
		
		energy += (ip_probs[AA] if  AAint2char(AA) != 'G' else 0.0)
			
	return (energy, seq)

def calc_top_seqs(in_res_path, probs_dir, num_sequences, outfile):
	probs_files = os.listdir(probs_dir)
	all_probs = [[], []] #[[[PROBS FOR IP1], [PROBS FOR IP2], ...], [[IP1], [IP2],...]]
	 
	#read back from .res file where ligand residues start
	lines = filelines2list(in_res_path)
	residue_positions = [int(line.strip()) - 1 for line in lines]
	
	for file in probs_files:
		file_path = os.path.join(probs_dir, file)
		lines = filelines2deeplist(file_path)

		probs = [float(x[1]) for x in lines]
		all_probs[0].append(probs)
		all_probs[1].append([int(file.split('_')[0]) - 1, int(file.split('_')[1]) - 1])
   
	#the following now fills ordered_probs
	ordered_probs = [] #[[sorted list of AA probabilities for a specific residue position: (0.3, 'A'), (0.5, 'C'),...], [like before for residue pos 2], etc...]
	AAs = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
	for residue in residue_positions:
		indices = [idx for idx in range(len(all_probs[1])) if residue in all_probs[1][idx]]
		probs = [all_probs[0][i] for i in indices]
		probs_sum = map(sum, zip(*probs))
		sorted_probs = sorted(zip(probs_sum, AAs), key=lambda pair: pair[0])
		ordered_probs.append(sorted_probs)
	
	#enumerate all the balls in urns possibilities
	#note that to maintain the max-heap property in Python's heapq implementation (which only does min-heap), we multiple by -1 before adding to heap
	num_urns = len(residue_positions)
	heap = []
	
	seq = ''
	energy = 0
	for i, c in enumerate(residue_positions):
		if(ordered_probs[i]):
			aa = ordered_probs[i][0][1]
			seq += aa
			energy += (0.0 if aa == 'G' else ordered_probs[i][0][0])
		else:
			seq += 'X'
	heapq.heappush(heap, (energy, seq))
	
	max_num_balls = 0
	total = 0
	while(total < num_sequences):
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
	results = [list(t) for t in zip(*heap)]
	out = [results[1][i] + '\t' + str(-1*energy) for i, energy in enumerate(results[0])]
	list2file(out, outfile)

#based on a X/100 threshold, classifies top X% as strong binders, bottom X% as weak binders. Uses results file and last column of values in results file 
#SHOULD IDEALLY USE INVERSE NORM, NOW JUST HARDCODING IN THRESHOLD Z-SCORES FOR ALPHA=0.05
def classify(results_file, outfile, threshold=0.05):

	sorted_results = get_sorted_results(results_file)
	
	z_score_threshold_low = -1.645
	z_score_threshold_high = 1.645
	
	out = []
	for result in sorted_results:
		if(result[0] < z_score_threshold_low):
			out.append(result[1] + '\t' + str(result[0]) + '\t' + 'strong')
		elif(result[0] > z_score_threshold_high):
			out.append(result[1] + '\t' + str(result[0]) + '\t' + 'weak')
		else:
			out.append(result[1] + '\t' + str(result[0]) + '\t' + 'inconclusive')
			
	list2file(out, outfile)

#TODO: IMPLEMENT THIS WITH DICTIONARIES LATER
#TODO: FIX SEMICOLON BUSINESS
def get_confusion_matrix(in_res_path, classification_results_file, true_class_file, in_pdb_orig=None):
	class_results = filelines2deeplist(classification_results_file, skipComments = True, useTuples = False, skipEmptyLines = True)
	true_class = get_true_class(in_res_path, true_class_file, in_pdb_orig=in_pdb_orig)
	
	(TP, FP, TN, FN) = (0, 0, 0, 0)
	for pair in class_results:
		if(pair[2] == 'inconclusive'):
			continue
		
		else:
			if pair[0] in true_class:
				if(('strong' in pair[2]) and ('strong' in true_class[pair[0]])):
					TP += 1
				elif(('weak' in pair[2]) and ('weak' in true_class[pair[0]])):
					TN += 1
				elif(('strong' in pair[2]) and ('weak' in true_class[pair[0]])):
					FP += 1
				elif(('weak' in pair[2]) and ('strong' in true_class[pair[0]])):
					FN += 1
				
	return (TP, FP, TN, FN)

def get_roc_curve_data(in_res_path, results_file, true_class_file, class_results, in_pdb_orig=None):
	true_class = get_true_class(in_res_path, true_class_file, class_results, in_pdb_orig)
	results = read_results(results_file)
	
	roc_data = list()
	
	if(class_results != None):
		class_results_dict = read_results(class_results, valueIsNum=False)
	
	for seq in true_class:
		
		if(class_results != None and class_results_dict[seq] == 'inconclusive'):
			print('Discarded ' + seq)
			continue
		
		if('weak' in true_class[seq].lower()):
			class_type = 1
		elif('strong' in true_class[seq].lower()):
			class_type = 0
		else:
			continue
		
#		print (seq, true_class[seq])
			
		if seq in results:
			roc_data.append((class_type, results[seq]))
		else:
			print 'Error: ' + seq + ' not found in results.'
	
	print roc_data		
	return roc_data

def calc_auroc(in_res_path, results_file, true_class_file, class_results, in_pdb_orig=None):
	import pyroc
	roc_data = get_roc_curve_data(in_res_path, results_file, true_class_file, class_results, in_pdb_orig)
	roc = pyroc.ROCData(roc_data)
	return roc.auc()
	

def plot_roc_curve(in_res_path, results_file, true_class_file, class_results, in_pdb_orig=None, title='ROC CURVE'):
	import pyroc
	roc_data = get_roc_curve_data(in_res_path, results_file, true_class_file, class_results, in_pdb_orig)
	roc = pyroc.ROCData(roc_data)
	roc.plot(title)

#from true classification key (text file of format, SEQ\tCLASSIFICATION of 'weak' or 'strong'),
#returns dictionary
#helper function
def get_true_class(in_res_path, true_class_file, class_results=None, in_pdb_orig=None):
	
	true_class = filelines2deeplist(true_class_file, skipComments = True, useTuples = False, skipEmptyLines = True)
	true_class_dict = dict()
		
	#read back from .res file where Chain B starts
	lines = filelines2list(in_res_path)
	residue_positions = [int(line.strip()) - 1 for line in lines]
			
	for pair in true_class:
		seq = pair[0] if (len(pair) == 2) else pair[0]+' '+pair[1]
		seq = fix_sequence_line(seq, len(residue_positions), in_pdb_orig)
		true_class_dict[seq] = pair[1] if (len(pair) == 2) else pair[2]
		
	return true_class_dict

#helper function
def get_sorted_results(results_file):
	results = filelines2deeplist(results_file, skipComments=True, useTuples=True, skipEmptyLines=True)
	filtered_results = []
	for result in results:
		#remember seq stored in results file is already 'fixed' so we don't need to fix sequence it again
		filtered_results.append((float(result[-1]), result[0]))
	
	sorted_results = sorted(filtered_results)
	
	return sorted_results

#function to combine sequence energy and true classification files so that we can easily read file in format: SEQ	   ENERGY  TRUE_CLASS
def summarize(in_energy_file_path, in_true_class_file_path, in_res_path, in_pdb_orig, out_file = '1ycr_mdm2_summarize.txt'):
	lines = open(in_energy_file_path).readlines()
	lines2 = open(in_true_class_file_path).readlines()
	lines = [line for line in lines if line[0] != '#' and line[0] != '\n']
	lines2 = [line for line in lines2 if line[0] != '#' and line[0] != '\n' and line[0] != '\r']
	lines = [line[:-1] for line in lines]
	table1 = [line.split('\t') for line in lines]
	lines2 = [line[:-2] for line in lines2]
	table2 = [line.split('\t') for line in lines2 if not line.isspace()]

	res_lines = filelines2list(in_res_path)
	residue_positions = [int(line.strip()) - 1 for line in res_lines]
	table2 = [[fix_sequence_line(pair[0], len(residue_positions), in_pdb_orig),pair[1]] for pair in table2]

	out = []
	for pair in table1:
		for i, pair2 in enumerate(table2):
			if pair2[0] == pair[0]:
				s = pair[0] + ' ' + pair[1] + ' ' +  pair2[1]
				out.append(s)
				break

	out_parts = [(triplet.split(' ')[0], float(triplet.split(' ')[1]), triplet.split(' ')[2]) for triplet in out]
	out_parts = sorted(out_parts, key=itemgetter(1), reverse=True) 
	if out_file is not None:
		list2file([part[0] + ' ' + str(part[1]) + ' ' + part[2] for part in out_parts], out_file)
	return out_parts

