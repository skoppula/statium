import os
from util import write_list_to_file

def analysis_wrapper(in_cfg_path, in_res_path, in_pdb_path, in_pdb_lib_dir, out_dir, verbose):
    
    if(verbose): print("Preparing directory folders...")
    list_lib_pdbs_path = prepare_directory(in_cfg_path, in_res_path, in_pdb_path, in_pdb_lib_dir, out_dir)
    
#    Joe uses a chmod'd shell script to do run statium_sidechain (AKA run_analysis):
#        sp = os.path.join(out_dir, 'seq_' + str(1) + '.sh')
#        of.write('#PBS -S /bin/sh\n' + runp + ' -statium_sidechain ' + preset + ' ' + seqp + ' ' + out_dir + ' ' + str(1) + '\n')
#        os.system('chmod u+x ' + sp)
#        os.system(sp)
#
#    v1.0.0's input: bfl1_2vm6    seq_1.txt    bgl1_2vm6_coyote    1
    run_analysis(in_cfg_path, in_res_path, in_pdb_path, list_lib_pdbs_path, 1)

def prepare_directory(in_cfg_path, in_res_path, in_pdb_path, in_pdb_lib_dir, out_dir, verbose):

    lib_pdbs = []
    pdbs = os.listdir(in_pdb_lib_dir);
    
    for pdb in pdbs:
        lib_pdbs.append(os.path.join(in_pdb_lib_dir, pdb))
    
    if(not os.path.exists(out_dir)):
        os.mkdir(out_dir)
    
    list_lib_pdbs_path = os.path.join(out_dir, 'pdbs.txt')
    write_list_to_file(lib_pdbs, list_lib_pdbs_path)
    
    return list_lib_pdbs_path

'''
def run_analysis(in_cfg_path, in_res_path, in_pdb_path, list_lib_pdbs_path, 1):
    

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