<b>STATIUM: smart scoring. promising proteins.</b><br>
STATIUM is an ongoing project at the Keating Lab to quantitatively understand how amino acid sequences interact. This repository implements a still-under-development algorithm we call STATIUM that effectively scores how well two proteins bind.

<b>Details and Documentation</b>
***
`precompute`:<br>
<i>Template</i>: `python wrapper.py precompute (--in_pdb --in_pdb_lib_dir --in_ip_lib_dir) [--out_dir]`<br>

<i>Example</i>:  `python wrapper.py precompute --in_pdb=1ycr_mdm2.pdb --in_pdb_lib_dir=data/culled_90/ --in_ip_lib_dir=data/ip_90_wGLY --out_dir=testing/output`<BR>

<i>Information</i>: Runs STATIUM analysis to create weights for each possible amino acid at the protein's interacting pair positions. Note that this command is a shortcut combination of the `renumber`,`create_res`, `run_statium`, and `get_orig_seq` commands.
***
'renumber':<br>
<i>Template</i> 'python wrapper.py renumber (--in_pdb) [--out_pdb --SRN --SAN]'<br>

<i>Example</i> 'python wrapper.py renumber testing/1ycr_mdm2_orig.pdb testing/1ycr_mdm2_new.pdb

<i>Information</i>: Internal function. Takes a PDB file, strips away the meta-data, and renumbers the residues and atoms. Renumbering starts on the first valid line of the PDB file, at starting atom number = SAN and starting residue number = SRN. 'Valid line' is any PDB line with 'ATOM' or 'HETATM' with 'MSE' (selenomethionine).
***
<b>Helpful Hints</b>:
+ Verbose output is turned on by default. To turn verbose output off, include the '-nv' or '--noverbose' flag.
+ Arguments wrapped in parenthesis () are required; arguments wrapped in square brackets [] are optional.

<b>Thanks for using STATIUM!</b>:
			   statium_wrapper.py create_res (IN_PDB_ORIG IN_PDB_RENUMBERED) [OUT_RES] [-v | --verbose]
			   statium_wrapper.py run_statium (IN_RES IN_PDB IN_PDB_LIB_DIR IN_IP_LIB_DIR) [OUT_DIR] [-v | --verbose]
			   statium_wrapper.py [-f] calc_energy (IN_RES PROBS_DIR SEQ_OR_FILE) [OUT_FILE] [--IN_PDB_ORIG=None] [-z | --zscores] [-p | --percentiles] [-v | --verbose] [-d | --draw_histogram]
			   statium_wrapper.py get_orig_seq (IN_PDB_ORIG) [-v | --verbose]
			   statium_wrapper.py [-f] generate_random_seqs (SEQ_LENGTH NUM_SEQS) [--OUT_FILE=None] [--TOTAL_PROTEIN_LIBRARY=None] [-v | --verbose]
			   statium_wrapper.py calc_top_seqs (IN_RES PROBS_DIR N) [OUT_FILE] [-v | --verbose]
			   statium_wrapper.py classify (RESULTS_FILE) [OUT_FILE] [ALPHA_THRESHOLD] [-v | --verbose]
			   statium_wrapper.py get_confusion_matrix (IN_RES CLASS_RESULTS TRUE_CLASS) [OUT_FILE] [--IN_PDB_ORIG=None] [-v | --verbose]
			   statium_wrapper.py [-i] calc_auroc (IN_RES RESULTS_FILE TRUE_CLASS) [OUT_FILE] [--IN_PDB_ORIG=None] [--CLASS_RESULTS=None] [-v | --verbose]
			   statium_wrapper.py [-i] plot_roc_curve (IN_RES RESULTS_FILE TRUE_CLASS) [--IN_PDB_ORIG=None] [--CLASS_RESULTS=None] [-v | --verbose]
			   statium_wrapper.py [-h | --help]

An example of an analysis sequence: <br>

	python statium_wrapper.py renumber 
	python statium_wrapper.py create_res testing/1ycr_mdm2_orig.pdb testing/1ycr_mdm2.pdb testing/1ycr_mdm2.res -v
	python statium_wrapper.py run_statium testing/1ycr_mdm2.res testing/1ycr_mdm2.pdb data/culled_90/ data/ip_90_wGLY/ testing/output -v
	python statium_wrapper.py get_orig_seq testing/1ycr_mdm2_orig.pdb
	
	python statium_wrapper.py -v generate_random_seqs 11 10
	python statium_wrapper.py -fv generate_random_seqs 11 10 random_seqs.txt
	python statium_wrapper.py -fv generate_random_seqs 11 10 random_seqs.txt --TOTAL_PROTEIN_LIBRARY=data/all_protein_sequences.txt
	
	python statium_wrapper.py -v calc_top_seqs testing/1ycr_mdm2.res testing/1ycr_mdm2_output_probs/ 10
	
	python statium_wrapper.py calc_energy testing/1ycr_mdm2.res testing/1ycr_mdm2_output_probs/ ETFSDLWKLLP
	python statium_wrapper.py calc_energy -fv testing/1ycr_mdm2.res testing/1ycr_mdm2_output_probs/ testing/1ycr_mdm2_sequences.txt testing/1ycr_mdm2_seq_energies.txt
	
	#-d flag draws histogram (requires a z (z-score) or p (percentile) flag simultaneously)
	python statium_wrapper.py calc_energy -fv testing/1ycr_mdm2.res testing/1ycr_mdm2_output_probs/ testing/1ycr_mdm2_sequences.txt testing/1ycr_mdm2_seq_energies.txt --IN_PDB_ORIG=testing/1ycr_mdm2_orig.pdb
	python statium_wrapper.py calc_energy -fzpv testing/1ycr_mdm2.res testing/1ycr_mdm2_output_probs/ testing/1ycr_mdm2_sequences.txt testing/1ycr_mdm2_seq_zscores_percentiles.txt --IN_PDB_ORIG=testing/1ycr_mdm2_orig.pdb
	python statium_wrapper.py calc_energy -fzv testing/1ycr_mdm2.res testing/1ycr_mdm2_output_probs/ testing/1ycr_mdm2_sequences.txt testing/1ycr_mdm2_seq_zscores.txt --IN_PDB_ORIG=testing/1ycr_mdm2_orig.pdb
	python statium_wrapper.py calc_energy -fzvd testing/1ycr_mdm2.res testing/1ycr_mdm2_output_probs/ testing/1ycr_mdm2_sequences.txt testing/1ycr_mdm2_seq_zscores.txt --IN_PDB_ORIG=testing/1ycr_mdm2_orig.pdb
	
	
	python statium_wrapper.py -v classify testing/1ycr_mdm2_seq_zscores.txt testing/1ycr_mdm2_classify_results_0.05.txt ALPHA_THRESHOLD=0.05
	python statium_wrapper.py -v get_confusion_matrix testing/1ycr_mdm2.res testing/1ycr_mdm2_classify_results_0.05.txt testing/1ycr_mdm2_seqs_true_classification.txt --IN_PDB_ORIG=testing/1ycr_mdm2_orig.pdb
	
	python statium_wrapper.py -v calc_auroc testing/1ycr_mdm2.res testing/1ycr_mdm2_seq_zscores.txt testing/1ycr_mdm2_seqs_true_classification.txt testing/1ycr_mdm2_auroc_zscores.txt --IN_PDB_ORIG=testing/1ycr_mdm2_orig.pdb
	python statium_wrapper.py -v calc_auroc testing/1ycr_mdm2.res testing/1ycr_mdm2_seq_energies.txt testing/1ycr_mdm2_seqs_true_classification.txt testing/1ycr_mdm2_auroc_energies.txt --IN_PDB_ORIG=testing/1ycr_mdm2_orig.pdb
	
	python statium_wrapper.py -iv calc_auroc testing/1ycr_mdm2.res testing/1ycr_mdm2_seq_zscores.txt testing/1ycr_mdm2_seqs_true_classification.txt testing/1ycr_mdm2_auroc_zscores_no_inconclusives.txt --IN_PDB_ORIG=testing/1ycr_mdm2_orig.pdb --CLASS_RESULTS=testing/1ycr_mdm2_classify_results_0.05.txt
	python statium_wrapper.py -iv calc_auroc testing/1ycr_mdm2.res testing/1ycr_mdm2_seq_energies.txt testing/1ycr_mdm2_seqs_true_classification.txt testing/1ycr_mdm2_auroc_energies_no_inconclusives.txt --IN_PDB_ORIG=testing/1ycr_mdm2_orig.pdb --CLASS_RESULTS=testing/1ycr_mdm2_classify_results_0.05.txt
	
	python statium_wrapper.py -v plot_roc_curve testing/1ycr_mdm2.res testing/1ycr_mdm2_seq_zscores.txt testing/1ycr_mdm2_seqs_true_classification.txt --IN_PDB_ORIG=testing/1ycr_mdm2_orig.pdb
	python statium_wrapper.py -v plot_roc_curve testing/1ycr_mdm2.res testing/1ycr_mdm2_seq_energies.txt testing/1ycr_mdm2_seqs_true_classification.txt --IN_PDB_ORIG=testing/1ycr_mdm2_orig.pdb
	
	python statium_wrapper.py -iv plot_roc_curve testing/1ycr_mdm2.res testing/1ycr_mdm2_seq_zscores.txt testing/1ycr_mdm2_seqs_true_classification.txt --IN_PDB_ORIG=testing/1ycr_mdm2_orig.pdb --CLASS_RESULTS=testing/1ycr_mdm2_classify_results_0.05.txt
	python statium_wrapper.py -iv plot_roc_curve testing/1ycr_mdm2.res testing/1ycr_mdm2_seq_energies.txt testing/1ycr_mdm2_seqs_true_classification.txt --IN_PDB_ORIG=testing/1ycr_mdm2_orig.pdb --CLASS_RESULTS=testing/1ycr_mdm2_classify_results_0.05.txt

#import this function and call with appropriate arguments to combine sequence energy and true classification into readable file
summarize(in_energy_file_path, in_true_class_file_path, in_res_path, in_pdb_orig, out_file = '1ycr_mdm2_summarize.txt')
	
