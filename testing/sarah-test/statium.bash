#!/bin/bash

input_PDB=1mhp_AB.pdb
output_PDB=1mhp_AB_out.pdb
res_file=1mhp_AB_out.res
out_folder="1mhp/"

# python statium_wrapper.py -h 
# renumber the input pdb 
#renumber (IN_PDB) [OUT_PDB] [-v | --verbose]
python statium_wrapper.py renumber  $input_PDB $output_PDB -v       

# create a residue file
#create_res (IN_PDB_ORIG IN_PDB_RENUMBERED) [OUT_RES] [-v | --verbose]
python statium_wrapper.py create_res $input_PDB $output_PDB  $res_file   


#run_statium (IN_RES IN_PDB IN_PDB_LIB_DIR IN_IP_LIB_DIR) [OUT_DIR] [-v | --verbose]
python statium_wrapper.py run_statium $res_file   $output_PDB  \
   /home/bartolo/web/PDB/culled_90s /home/bartolo/web/PDB/ip_90_wGLY  $out_folder  -v
#statium_wrapper.py [-f] calc_energy (IN_RES PROBS_DIR SEQ_OR_FILE) [OUT_FILE] [--IN_PDB_ORIG=None] [-z | --zscores] [-p | --percentiles] [-v | --verbose]
#statium_wrapper.py get_orig_seq (IN_PDB_ORIG) [-v | --verbose]
#statium_wrapper.py [-f] generate_random_seqs (SEQ_LENGTH NUM_SEQS) [--OUT_FILE=None] [--TOTAL_PROTEIN_LIBRARY=None] [-v | --verbose]
#statium_wrapper.py calc_top_seqs (IN_RES PROBS_DIR N) [OUT_FILE] [-v | --verbose]
#statium_wrapper.py classify (RESULTS_FILE) [OUT_FILE] [ALPHA_THRESHOLD] [-v | --verbose]
#statium_wrapper.py get_confusion_matrix (IN_RES CLASS_RESULTS TRUE_CLASS) [OUT_FILE] [--IN_PDB_ORIG=None] [-v | --verbose]

#statium_wrapper.py [-i] calc_auroc (IN_RES RESULTS_FILE TRUE_CLASS) [OUT_FILE] [--IN_PDB_ORIG=None] [--CLASS_RESULTS=None] [-v | --verbose]

#statium_wrapper.py [-i] plot_roc_curve (IN_RES RESULTS_FILE TRUE_CLASS) [--IN_PDB_ORIG=None] [--CLASS_RESULTS=None] [-v | --verbose]

#statium_wrapper.py [-h | --help]



# job ends here; spit time:
echo "Job end time:  `date`" 


