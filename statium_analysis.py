import os
import math
from util import filelines2list
from util import list2file
from util import isAA
from util import AA2char
from util import char2AA
from util import AAchar2int
from util import AAint2char
from util import get_sidechain_atoms

def analysis_pipeline(in_res_path, in_pdb_path, in_pdb_lib_dir, in_ip_lib_dir, out_dir, verbose):
    
    if(verbose): print("Preparing directory folders...")
    lib_pdbs_path = prepare_directory(in_res_path, in_pdb_path, in_pdb_lib_dir, out_dir)
    if(verbose): print("Written list of library PDB paths to: " + lib_pdbs_path)
    
#    Joe uses a chmod'd shell script to do run statium_sidechain (AKA run_analysis):
#        sp = os.path.join(out_dir, 'seq_' + str(1) + '.sh')
#        of.write('#PBS -S /bin/sh\n' + runp + ' -statium_sidechain ' + preset + ' ' + seqp + ' ' + out_dir + ' ' + str(1) + '\n')
#        os.system('chmod u+x ' + sp)
#        os.system(sp)
#
#    v1.0.0's input: bfl1_2vm6    seq_1.txt    bgl1_2vm6_coyote    1
    run_analysis(in_res_path, in_pdb_path, lib_pdbs_path, in_ip_lib_dir, out_dir, 1, 6.0, 4.0, verbose)


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


def run_analysis(in_res_path, in_pdb_path, lib_pdbs_path, in_ip_lib_dir, out_dir, num, pair_dist_cutoff, sidechain_match_cutoff, verbose):
    
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
    
    d = [[0 for j in range(20)] for i in range(len(use_indices))]
    
    template_distances = []
    for pair in use_indices:
        (pos1, pos2) = (pair[0], pair[1])
        (pos1_list, pos2_list) = (pdb_info[pos1][2][0], ['CA', 'CB']) 
        
        if ('CA' in pdb_info[pos2][2][0] and 'CB' in pdb_info[pos2][2][0]):
            template_distances.append([])
        else:
            template_distances.append(select_sidechain_distances(pos1_list, pos2_list, distance_matrix[pos1][pos2]))
    
        
    ip_res = [0 for i in range(20)]
    lib_pdb_paths = filelines2list(lib_pdbs_path)
    
    
    for (i, pdb_path) in enumerate(lib_pdb_paths):    
        if(verbose): print("Processing library .pdb: " + pdb_path + "\t (" + str(i) + " out of " + str(len(lib_pdb_paths)) + ")")
        lib_ip_path = os.path.join(in_ip_lib_dir, os.path.split(pdb_path)[1].split('.')[0] + '.ip')
        lib_pdb_info = get_pdb_info(pdb_path)        
        lib_use_indices = get_IPs(lib_ip_path)
        lib_distance_matrix = distance_matrix_sidechain_use(lib_pdb_info, lib_use_indices)

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
    

def get_IPs(ip_file):
    lines = filelines2list(ip_file)
    return map(extractIP, lines)

def extractIP(line):
    items = line.strip().split()
    (pos1, pos2) = (int(items[0]), int(items[1]))
    return [pos1, pos2]

#just choose any distance with similar interacting atom pairs as given
#I don't think this actually returns meaningful results because it's using
#atom indices to query residues, but shouldn't matter since it's 
#only used in a filler/else/except clause
def select_sidechain_distances(pos1_atoms, pos2_atoms, distance_matrix):
  
    pair_info = [[], []]
    for atomi in pos1_atoms:
        for atomj in pos2_atoms:
            idx = distance_matrix[0].index([atomi, atomj])
            pair_info[0].append([atomi, atomj])
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
                
                if position not in res_list:
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
        
'''
def statium_sidechain(preset_dir, pdb_paths, out_dir, num):
  
    ip = False
    pair_def_cutoff = 6.0
    sidechain_match_cutoff = 0.4
  
    paths = lines2list(pdb_paths)
  
    file_dir = os.path.split(preset_dir)[0]
    file_base = os.path.split(preset_dir)[1]
    
    template_pdb_path = os.path.join(file_dir, file_base + '.pdb')
    residue_path = os.path.join(file_dir, file_base + '.res')
    mask_path = os.path.join(file_dir, file_base + '.mask')
    ip_path = os.path.join(file_dir, file_base + '.ip')
    
    res_vec = []
    res_lines = readlines(residue_path)
    for line in res_lines:
        res_vec.append(int(line.strip()) - 1)
       
    mask_vec = []
    if os.path.exists(mask_path):
        mask_lines = readlines(mask_path)
        for line in mask_lines:
            mask_vec.append([int(line.strip().split()[0]) - 1, int(line.strip().split()[1]) - 1])
    
    pdbinfo_vec = store_pdb_info(template_pdb_path)
    N = len(pdbinfo_vec)
    template_distance_matrix = distance_matrix_sidechain(pdbinfo_vec)

    if os.path.exists(ip_path) and ip:
        use_index = []
        ip_lines = readlines(ip_path)
        for line in ip_lines:
            pos0 = line.split()[0]
            pos1 = line.split()[1]
            use_index.append([int(pos0) - 1, int(pos1) - 1])
    else:
        use_index = []
        for i in range(N):
        for j in range(i + 1, N):
            if not i in res_vec and not j in res_vec: continue
            if [i, j] in mask_vec: continue
            if i in res_vec and j in res_vec: continue
            #if AAChar_fasta(pdbinfo_vec[i][1]) == 'G' or AAChar_fasta(pdbinfo_vec[j][1]) == 'G': continue
            if atoms_within_cutoff(template_distance_matrix[i][j], pair_def_cutoff):
            use_index.append([i, j])
    
    Nuse = len(use_index)

    dcounts = []
    for i in range(Nuse):
    dcounts.append([])
    for j in range(20):
        dcounts[i].append(0)
        
    template_distances = []
    for i in range(Nuse):
    pos1 = use_index[i][0]
    pos2 = use_index[i][1]
    pos1_list = pdbinfo_vec[pos1][2][0]
    pos2_list = ['CA', 'CB']
    if not stub_intact(pdbinfo_vec[pos2][2][0]): template_distances.append([])
    else: template_distances.append(select_sidechain_distances(pos1_list, pos2_list, template_distance_matrix[pos1][pos2], True))

    ip_res = []
    for i in range(20): ip_res.append(0)

    for path in range(len(paths)):
    print path
        pdb_path = paths[path][0]
        lib_ip_path = os.path.join('/home/bartolo/web/PDB/ip_90_wGLY', os.path.split(paths[path][0])[1].split('.')[0] + '.ip')
    lib_pdbinfo_vec = store_pdb_info(pdb_path)
        
        libN = len(lib_pdbinfo_vec)
    
    lib_use_index = upload_interacting_pairs(lib_ip_path)
    lib_distance_matrix = distance_matrix_sidechain_use(lib_pdbinfo_vec, lib_use_index)

    #lib_distance_matrix = distance_matrix_sidechain(lib_pdbinfo_vec)
        #lib_use_index = []
        #ip_file = open(lib_ip_path, 'w')
        #for i in range(libN):
     #   for j in range(i + 1, libN):
            #if AAChar_fasta(lib_pdbinfo_vec[i][1]) == 'G' or AAChar_fasta(lib_pdbinfo_vec[j][1]) == 'G': continue
      #      if atoms_within_cutoff(lib_distance_matrix[i][j], pair_def_cutoff):
    #        lib_use_index.append([i, j])
    #        ip_file.write(str(i) + '\t' + str(j) + '\n')
    #ip_file.close()
    #continue

    libNuse = len(lib_use_index)
    for i in range(libNuse):
        lib_pos1 = lib_use_index[i][0]
        lib_pos2 = lib_use_index[i][1]
        if lib_pos2 - lib_pos1 <= 4: continue
        lib_AA1 = lib_pdbinfo_vec[lib_pos1][1]
        lib_AA2 = lib_pdbinfo_vec[lib_pos2][1]
        if lib_AA1 < 0 or lib_AA2 < 0 or lib_AA1 > 19 or lib_AA2 > 19: continue
            ip_res[lib_AA1] += 1
            ip_res[lib_AA2] += 1
        for j in range(Nuse):
            pos1 = use_index[j][0]
            pos2 = use_index[j][1]
            AA1 = pdbinfo_vec[pos1][1]
            if stub_intact(lib_pdbinfo_vec[lib_pos2][2][0]):
                if lib_AA1 == AA1:
                       lib_dist_for = select_sidechain_distances(lib_pdbinfo_vec[lib_pos1][2][0], ['CA', 'CB'], lib_distance_matrix[lib_pos1][lib_pos2], True)
                if matching_sidechain_pair(template_distances[j], lib_dist_for, sidechain_cutoff_dist(AAChar_fasta(AA1))):
                dcounts[j][lib_AA2] += 1
            if stub_intact(lib_pdbinfo_vec[lib_pos1][2][0]):
                if lib_AA2 == AA1:
                       lib_dist_rev = select_sidechain_distances(['CA', 'CB'], lib_pdbinfo_vec[lib_pos2][2][0], lib_distance_matrix[lib_pos1][lib_pos2], False)
                if matching_sidechain_pair(template_distances[j], lib_dist_rev, sidechain_cutoff_dist(AAChar_fasta(AA1))):
                dcounts[j][lib_AA1] += 1

    counts_file = open(os.path.join(out_dir, 'ip_res.txt_' + str(num)), 'w')
    for i in range(20): counts_file.write(AAChar_fasta(i) + '\t' + str(ip_res[i]) + '\n')
    counts_file.close()

    for pos in range(Nuse):
          
    counts_file = open(os.path.join(out_dir, str(use_index[pos][0] + 1) + '_' + str(use_index[pos][1] + 1) + '_counts.txt_' + str(num)), 'w')
    for i in range(20): counts_file.write(AAChar_fasta(i) + '\t' + str(dcounts[pos][i]) + '\n')
    counts_file.close()
    
'''
#                 if sys.argv[i] == '-build_statium':
#         preset = sys.argv[i + 1]
#                 statium_sidechain_coyote(preset, 'REG', 5000, sys.argv[0])
#                 statium_sidechain_coyote_compile(preset)
#         convert_counts_sidechain(preset)