import os
import math
import heapq
import itertools
import pyroc
from util import filelines2list
from util import list2file
from util import isAA
from util import AA2char
from util import AAchar2int
from util import AAint2char
from util import get_sidechain_atoms
from util import AA_cutoff_dist
from util import filelines2deeplist
from util import binary_placement
from util import mean
from util import std
from util import nCr
from util import read_results
from statium_reformat import get_orig_seq
from statium_reformat import generate_random_seqs

def statium_pipeline(in_res_path, in_pdb_path, in_pdb_lib_dir, in_ip_lib_dir, out_dir, verbose):
    
    if(verbose): print("Preparing directory folders...")
    lib_pdbs_path = prepare_directory(in_res_path, in_pdb_path, in_pdb_lib_dir, out_dir)
    if(verbose): print("Written list of library PDB paths to: " + lib_pdbs_path)
    
    run_analysis(in_res_path, in_pdb_path, lib_pdbs_path, in_ip_lib_dir, out_dir, 6.0, verbose)
    if(verbose): print("Done.")

def prepare_directory(in_res_path, in_pdb_path, in_pdb_lib_dir, out_dir):

    lib_pdbs = []
    pdbs = os.listdir(in_pdb_lib_dir);
    
    for pdb in pdbs:
        lib_pdbs.append(os.path.join(in_pdb_lib_dir, pdb))
    
    if(not os.path.exists(out_dir)):
        os.mkdir(out_dir)
      
    lib_pdbs_path = os.path.join(out_dir, 'pdbs.txt')
    list2file(lib_pdbs, lib_pdbs_path)
    
    return lib_pdbs_path


def run_analysis(in_res_path, in_pdb_path, lib_pdbs_path, in_ip_lib_dir, out_dir, pair_dist_cutoff, verbose):
    
    res_lines = filelines2list(in_res_path)
    residues = [(int(line.strip()) - 1) for line in res_lines]
    
    pdb_info = get_pdb_info(in_pdb_path)
    if(verbose): print("Finished extracting information from the input PDB file: " + in_pdb_path)

    distance_matrix = compute_distance_matrix(pdb_info)
    if(verbose): print("Finished computing inter-atomic distances for residues in input PDB file")
    
    use_indices = [] #check which residue pairs to use (within interacting distance)
    for i in range(len(pdb_info)):
        for j in range(i+1, len(pdb_info)):
            if ((i in residues and j not in residues) or (j in residues and i not in residues)) and (check_cutoff(distance_matrix[i][j], pair_dist_cutoff)):
                use_indices.append([i, j])
    
    distances = [] #stores distance info for each use_indices
    for pair in use_indices:
        (pos1, pos2) = (pair[0], pair[1])
        (pos1_list, pos2_list) = (pdb_info[pos1][2][0], ['CA', 'CB']) 
        
        if not stub_intact(pdb_info[pos2][2][0]):
            distances.append([])
        else:
            distances.append(select_sidechain_distances(pos1_list, pos2_list, distance_matrix[pos1][pos2], "forward"))

    ip_res = [0 for i in range(20)]
    counts = [[0 for j in range(20)] for i in range(len(use_indices))]  #final total counts of each similar IP
    lib_pdb_paths = filelines2list(lib_pdbs_path)
    
    for (i, pdb_path) in enumerate(lib_pdb_paths):    
        if(verbose): print("Processing library .pdb: " + pdb_path + "\t (" + str(i) + " out of " + str(len(lib_pdb_paths)) + ")")
        lib_ip_path = os.path.join(in_ip_lib_dir, os.path.split(pdb_path)[1].split('.')[0] + '.ip')
        lib_pdb_info = get_pdb_info(pdb_path)        
        lib_use_indices = get_IPs(lib_ip_path)
        lib_distance_matrix = distance_matrix_sidechain_use(lib_pdb_info, lib_use_indices)

        for lib_pair in lib_use_indices:

            (lib_pos1, lib_pos2) = (lib_pair[0], lib_pair[1])
            
            if lib_pos2 - lib_pos1 <= 4: continue
            
            (lib_AA1, lib_AA2) = (int(lib_pdb_info[lib_pos1][1]), int(lib_pdb_info[lib_pos2][1]))

            if lib_AA1 < 0 or lib_AA2 < 0 or lib_AA1 > 19 or lib_AA2 > 19: continue
            
            ip_res[lib_AA1] += 1
            ip_res[lib_AA2] += 1
            
            for (j, pair) in enumerate(use_indices):
                
                (pos1, pos2) = (pair[0], pair[1])
                AA1 = int(pdb_info[pos1][1])
                
                if stub_intact(lib_pdb_info[lib_pos2][2][0]):
                    if lib_AA1 == AA1:
                        lib_dist_forward = select_sidechain_distances(lib_pdb_info[lib_pos1][2][0], ['CA', 'CB'], lib_distance_matrix[lib_pos1][lib_pos2], "forward")
                        if matching_sidechain_pair(distances[j], lib_dist_forward, AA_cutoff_dist(AAint2char(AA1))):
                            counts[j][lib_AA2] += 1
                        
                if stub_intact(lib_pdb_info[lib_pos1][2][0]):
                    if lib_AA2 == AA1:
                        lib_dist_rev = select_sidechain_distances(['CA', 'CB'], lib_pdb_info[lib_pos2][2][0], lib_distance_matrix[lib_pos1][lib_pos2], "backward") 
                        if matching_sidechain_pair(distances[j], lib_dist_rev, AA_cutoff_dist(AAint2char(AA1))):
                            counts[j][lib_AA1] += 1
    
    if(verbose): print("Finished processing library .pdb files. Writing counted results to files in directory: " + out_dir)
              
    counts_file = open(os.path.join(out_dir, 'lib_ip_residue_counts.txt'), 'w')
    for i in range(20): counts_file.write(AAint2char(i) + '\t' + str(ip_res[i]) + '\n')
    counts_file.close()

    for (i, pair) in enumerate(use_indices):
        counts_file = open(os.path.join(out_dir, str(pair[0] + 1) + '_' + str(pair[1] + 1) + '_counts.txt'), 'w')
        for j in range(20): counts_file.write(AAint2char(j) + '\t' + str(counts[i][j]) + '\n')
    counts_file.close()
    
    if(verbose): print("Determining probabilities from counts...")
    determine_probs(out_dir, verbose)

def determine_probs(out_dir, verbose):
    
    #create the _probs output directory
    if(out_dir[-1] == '/'):
        probs_dir = out_dir[-1] + '_probs'
    else:
        probs_dir = out_dir + '_probs'
        
    if not os.path.exists(probs_dir):
        os.mkdir(probs_dir)
    
    if(verbose): print("Reading in the total residue counts of the PDB library: lib_ip_residue_counts.txt")
    output_files = os.listdir(out_dir)
    lib_summary_path = filter(lambda x: 'lib_ip_residue' in x, output_files) #search files for lib_ip_residue_counts.txt
    lib_sum_data = filelines2deeplist(os.path.join(out_dir, lib_summary_path[0]))
    
    lib_pdb_total = float(sum([int(x[1]) for x in lib_sum_data]))
    lib_AA_probs = [float(x[1]) / lib_pdb_total for x in lib_sum_data]
    
    #Transform every count file into a _probs file    
    for file in output_files:
        if('count' in file and len(file.split('_')) == 3):
            path = os.path.join(out_dir, file)
            counts = filelines2deeplist(path)
            out_path = os.path.join(probs_dir, file.split('_')[0] + '_' + file.split('_')[1] + '_probs.txt')

            total = float(sum([int(x[1]) for x in counts]))

            if(total > 99):
                out = open(out_path, 'w')
                AA_probs = [(float(x[1])/total if int(x[1]) != 0 else 1/total) for x in counts]
                
                for i in range(20):
                    e = -1.0 * math.log(AA_probs[i] / lib_AA_probs[i])
                    out.write(AAint2char(i) + '\t' + str(e) + '\n')
                
                out.close()   
    
    if(verbose): print("Finished calculating probabilities. Written to: " + out_dir + '_probs')

    
def matching_sidechain_pair(distances1, distances2, cutoff):
  
    sd = 0.0
    count = 0.0

    for i in range(len(distances1[0])):
        pair_i = distances1[0][i]
        di = distances1[1][i]
        
        for j in range(len(distances2[0])): 
            pair_j = distances2[0][j]
            dj = distances2[1][j]

            if pair_j == pair_i:
                sd += ((di - dj) ** 2)
                count += 1.0

    if math.sqrt(sd / count) < cutoff:
        return True
    else:
        return False
    

#checks that 'CA' and 'CB' are in list of atoms
def stub_intact(atoms):
    if 'CA' in atoms and 'CB' in atoms:
        return True
    else:
        return False

#returns a small distance matrix for the residues with given indices using the input pdb_info structure
def distance_matrix_sidechain_use(pdb_info, use_index):
    
    N = len(pdb_info)
    distance_matrix = [[[[], []] for i in range(N)] for j in range(N)]
    
    for pair in use_index:
        (i, j) = (pair[0], pair[1])
        
        for k in range(len(pdb_info[i][2][0])):
            for l in range(len(pdb_info[j][2][0])):
                (atom1, atom2) = (pdb_info[i][2][0][k], pdb_info[j][2][0][l])
                dist = distance(pdb_info[i][2][1][k], pdb_info[j][2][1][l])
                distance_matrix[i][j][0].append([atom1, atom2])
                distance_matrix[i][j][1].append(dist)
        
    return distance_matrix
    
#gets interacting pairs from .ip file
def get_IPs(ip_file):
    lines = filelines2list(ip_file)
    return map(extractIP, lines)

def extractIP(line):
    items = line.strip().split()
    (pos1, pos2) = (int(items[0]), int(items[1]))
    return [pos1, pos2]

#just choose any distance with similar interacting atom pairs as given
#I don't know if this actually returns meaningful results because it's using
#atom indices to query residues, but shouldn't matter since it's 
#only used in a filler/else/except clause
def select_sidechain_distances(pos1_atoms, pos2_atoms, distance_matrix, direction):
  
    pair_info = [[], []]
    for atomi in pos1_atoms:
        for atomj in pos2_atoms:
            idx = distance_matrix[0].index([atomi, atomj])
            if(direction == 'forward'):
                pair_info[0].append([atomi, atomj])
            else:
                pair_info[0].append([atomj, atomi])
            pair_info[1].append(distance_matrix[1][idx])
    
    return pair_info
    
def check_cutoff(residue_pair, cutoff):
    for k in range(len(residue_pair[1])):
        if residue_pair[1][k] < cutoff: return True
        
    return False


def compute_distance_matrix(pdb_info):
    
    N = len(pdb_info)
    distance_matrix = [ [ [[], []] for i in range(N) ] for j in range(N)]
    
    for i in range(N):
        for j in range(i + 1, N):
            for k in range(len(pdb_info[i][2][0])):
                for l in range(len(pdb_info[j][2][0])):
                    distance_matrix[i][j][0].append([pdb_info[i][2][0][k], pdb_info[j][2][0][l]])
                    distance_matrix[i][j][1].append(distance(pdb_info[i][2][1][k], pdb_info[j][2][1][l]))
        
    return distance_matrix

#just apply distance formula
def distance(C1, C2):
    return math.sqrt(((C1[0] - C2[0])**2) + ((C1[1] - C2[1])**2) + ((C1[2] - C2[2])**2))


#Creates a list with information for each residue:
#e.g for each residue: [[1, ''], '16', [['CA', 'CB', 'OG1', 'CG2'], [[21.142, -19.229, 4.185], [21.957, -18.596, 5.322], [23.023, 17.818, 4.773], [22.547, -19.67, 6.206]]]]
def get_pdb_info(in_pdb_path):

    pdb_info = [] #list of lists: contains info on each AA
    res_list = []
 
    lines = filelines2list(in_pdb_path)
    
    for i, line in enumerate(lines):
        if line[0:4] == 'ATOM' and line[13:15] == 'CA':
            try:
                position = int(line[22:28].strip())
                chainID = line[21:22].strip()
                
                if position in res_list:
                    continue
                res_list.append(position)
                
            except:
                continue
            
            AA = line[17:20]
            if isAA(AA):
                
                pdb_info.append([])
                pdb_info[-1].append([position, chainID])        #Put in position
                pdb_info[-1].append(AAchar2int(AA2char(AA)))    #Put in AA identity
                
                pdb_info[-1].append([[], []])                   #Put in alpha carbon, and coordinates of alpha carbon
                pdb_info[-1][-1][0].append('CA')
                pdb_info[-1][-1][1].append([float(line[30:38]), float(line[39:46]), float(line[47:54])])
                
                atoms_list = get_sidechain_atoms(AA2char(AA))
                found_list = []
                
                for line2 in lines[i: len(lines)]:
                    if line2[0:4] == 'ATOM':
                        try:
                            position2 = int(line2[22:28].strip())
                            chainID2 = line2[21:22].strip()
                        except: continue
                        
                        if position2 > position:
                            break
                        elif position2 == position and chainID2 == chainID:
                            atom_type = line2[13:16].strip()
                        
                            #Put in other atoms, and coordinates of atoms
                            if atom_type in atoms_list and not atom_type in found_list:
                                found_list.append(atom_type)                            
                                pdb_info[-1][-1][0].append(atom_type)
                                pdb_info[-1][-1][1].append([float(line2[30:38]), float(line2[39:46]), float(line2[47:54])])

    for t in pdb_info:
        if len(t) != 3:
            print 'INCORRECT FORMAT: ' + in_pdb_path + ': ' + t
        
        for point in t[2][1]:
            for k in range(3):
                try: float(point[k])
                except: print 'Bad coordinate at ' + str(point[k]) + ' in ' + in_pdb_path
    
    return pdb_info

def generate_random_distribution (in_res, in_probs_dir, num_seqs=4000):
    
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

def fix_sequence_line(seq, desired_seq_length, in_pdb_orig=None):
    parts = seq.split()
    
    if(len(seq) != desired_seq_length and len(parts) == 1):
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
    
    #read back from .res file where Chain B starts
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
        pos1 = ip_pos[1] #peptide position on chain B
        
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
     
    #read back from .res file where Chain B starts
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
    true_class = get_true_class(in_res_path, true_class_file, in_pdb_orig)
    
    (TP, FP, TN, FN) = (0, 0, 0, 0)
    for pair in class_results:
        if(pair[2] == 'inconclusive'):
            continue
        
        else:
            if pair[0] in true_class:
                if(pair[2] == 'strong' and true_class[pair[0]] == 'strong'):
                    TP += 1
                elif(pair[2] == 'weak' and true_class[pair[0]] == 'weak'):
                    TN += 1
                elif(pair[2] == 'strong' and true_class[pair[0]] == 'weak'):
                    FP += 1
                elif(pair[2] == 'weak' and true_class[pair[0]] == 'strong'):
                    FN += 1
            
#            print('FOUND MATCH', true_class[idx], pair)
    
    return (TP, FP, TN, FN)

def calc_auroc(in_res_path, results_file, true_class_file, in_pdb_orig=None):
    true_class = get_true_class(in_res_path, true_class_file, in_pdb_orig)
    results = read_results(results_file)
    
    roc_data = list()
    
    for seq in true_class:
        class_type = 0 if true_class[seq].lower() == 'weak' else 1
        if seq in results:
            roc_data.append((class_type, results[seq]))
        else:
            print 'Error: ' + seq + ' not found in results.'
    
    roc = pyroc.ROCData(roc_data)
    return roc.auc()
    

def plot_roc_curve(in_res_path, results_file, true_class_file, in_pdb_orig=None, title='ROC CURVE'):
    true_class = get_true_class(in_res_path, true_class_file, in_pdb_orig)
    results = read_results(results_file)
    
    roc_data = list()
    
    for seq in true_class:
        class_type = 0 if true_class[seq].lower() == 'weak' else 1
        if seq in results:
            roc_data.append((class_type, results[seq]))
        else:
            print 'Error: ' + seq + ' not found in results.'
    
    roc = pyroc.ROCData(roc_data)
    roc.plot(title)

#from true classification key (text file of format, SEQ\tCLASSIFICATION of 'weak' or 'strong'),
#process input sequences and output lists of [[seq, class], [..],...]
#returns dictionary
#helper function
def get_true_class(in_res_path, true_class_file, in_pdb_orig=None):
    true_class = filelines2deeplist(true_class_file, skipComments = True, useTuples = False, skipEmptyLines = True)
    true_class_dict = dict()
    
    #read back from .res file where Chain B starts
    lines = filelines2list(in_res_path)
    residue_positions = [int(line.strip()) - 1 for line in lines]
    
    for i, pair in enumerate(true_class):
        seq = pair[0] if (len(pair) == 2) else pair[0]+' '+pair[1]
        seq = fix_sequence_line(seq, len(residue_positions), in_pdb_orig)
        true_class_dict[seq] = true_class[i][1]
        
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
