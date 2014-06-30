<b>STATIUM: smart scoring. promising proteins.</b><br>
STATIUM is an ongoing project at the Keating Lab to quantitatively understand how amino acid sequences interact. This repository implements a still-under-development algorithm we call STATIUM that scores how well two or more proteins bind at their interacting positions.

<b>Details and Documentation</b>
***
`renumber`:<br>
<i>Template</i> `python wrapper.py renumber (--in_pdb) [--out_pdb --SRN --SAN --chains]`<br>

<i>Example</i> `python wrapper.py renumber --in_pdb=testing/1mph_HLA.pdb --out_pdb=testing/1mph_AHL.pdb --chains=H,L` <br>

<i>Specifics</i>: Takes a PDB file, strips away the meta-data, and renumbers the residues and atoms, retaining atom coordinate positions and removing the occupancy and temperature factors. Renumbering starts on the first valid line of the PDB file, at starting atom number = SAN and starting residue number = SRN. 'Valid line' is any PDB line with 'ATOM' or 'HETATM' with 'MSE' (selenomethionine).

One requirement of STATIUM as implemented here is that the input PDB has the receptor sequence listed before the peptide/mutant binder sequences in the file. The `renumber` function outputs a PDB that follows this requirement by reading the chains you wish to designate as the ligand sequence from the --chains option. In the example above, the heavy and light chain that we want to mutate were pushed after chain A (*HLA to *AHL).
***
`create_res`:<br>
<i>Template</i> `python wrapper.py create_res (--in_pdb_orig --in_pdb_renum) [--out_res --position_pairs]`<br>

<i>Example</i> `python wrapper.py create_res --in_pdb_orig=testing/1mph_AHL_orig.pdb --in_pdb_renum=testing/1mph_AHL_renum.pdb --out_res=testing/1mph_AHL.res --position_pairs=L1-20,H33`<br>
<i>Example 2</i> `python wrapper.py create_res --in_pdb_orig=testing/1mph_AHL_orig.pdb --in_pdb_renum=testing/1mph_AHL_renum.pdb --out_res=testing/1mph_AHL.res --position_pairs=H`<br>

<i>Specifics</i>: Takes in both the original and renumbered PDB files (see 'renumber'). It translates pairs of (chain identifier, number) uniquely demarcating a residues on the original PDB file to a number uniquely demarcating a residue in the renumbered file. These set of numbers are written to a file and used as the positions to be analyzed by the STATIUM algorithm.

The --position_pairs argument specifies which the set of positions to be included as binder/ligand sequences in the STATIUM analysis. The argument is a set of comma seperated terms which represent continuous sequence of residues to be included in the ligand sequence (inclusive). If you want the entirety of a chain, simply put the name of chain in the list (e.g. --position_pairs=H). In the first example above, residues on the L chain, position 1-20, and a residue on the H chain, position 33 will be included in the output residues file.

If you fail to include a --position_pairs argument, the function will assume you mean to create a *.res file with the entirety of chain <i>B</i>.
***
`preprocess`:<br>
<i>Template</i> `python wrapper.py preprocess (--in_dir) [--out_dir --ip_dist_cutoff]`<br>
<i>Example</i> `python wrapper.py preprocess -r --in_dir=data/culled_90 --out_dir=data/preprocessed_culled_90 --ip_dist_cutoff=5`<br>

<i>Specifics</i> Takes in a directory of library PDBs and outputs a directory of Python *.pickle files, one for each PDB. Each pickle file contains (1) a list of Residue objects parsed from the PDB file, (2) a matrix of inter-residue distances, and (3) a list of all interacting pairs of residues. Two residues are considered interacting if any of their atoms are within the IP cutoff distance (--ip_dist_cutoff). The default IP cutoff distance is 5 Angstroms.

The `-r` (restart) flag rewrite the pickle files previously created in the output directory. Not including it skips preprocessing PDB's whose pickle is already in the output directory. 
***
`run_statium`:<br>
<i>Template</i> `python wrapper.py run_statium (--in_pdb --in_res --pdb_lib) [--out_dir --ip_dist_cutoff --matching_res_dist_cutoffs --counts]`<br>

<i>Example</i> `python wrapper.py run_statium --in_pdb=testing/sarah-test/1mhp_AHL_new.pdb --in_res=testing/sarah-test/1mhp_AHL.res --pdb_lib=data/culled_90/ --ip_lib=data/ip_90_wGLY/` <br>

<i>Specifics</i>: Takes in a renumbered PDB file (see `renumber`), the directory of the preprocessed protein library containing each PDB's *.pickle file (--pdb_lib).

Optional parameters include: --out_dir (the directory where STATIUM outputs its results; default value is value of --in_pdb without the .pdb extension), --counts (whether to print out STATIUM's intermediate analysis outputs; note that this takes no argument; simply including the flag issues printing!), --ip_dist_cutoff (the threshold distance in Angstroms between two atoms, below which the atom's residues are deemed 'interacting'; default is 6.0), and --matching_res_dist_cutoff (a dictionary with all twenty amino acids [in single character representation] each mapped to a cutoff below which a interacting residue pair cannot be deemed 'matching' to a library protein interacting pair. Example of using this parameter [containing the default, recommended dictionary values if you leave this parameter out]: --matching_res_dist_cutoff={'A':2, 'C':6, 'D':6, 'E':6, 'F':6, 'G':2, 'H':6, 'I':6, 'K':6, 'L':6, 'M':6, 'N':6, 'P':6, 'Q':6, 'R':6, 'S':6, 'T':6, 'V':6, 'W':6, 'Y':6, 'X':0}.

The function creates a directory containing a set of files, one file per interacting pair:

Each file contains a set of twenty probabilities (one for each amino acid) describing how likely it is for that identity would exist at that position on the sidechain, given the main chain's amino acid identity at the position.
***
`energy`:<br>
<i>Template</i> `python wrapper.py energy (--in_res --in_probs) [-f] (--in_seqs) [--out] [-z | --zscores] [-p | --percentiles] [-v | --verbose] [-d | --draw_histogram]`<br>
<i>Example One</i> `python wrapper.py --in_res=testing/sarah-test/1mhp_AHL.res --in_probs=testing/sarah-testing/1mhp_AHL_probs --in_seqs=AAAGGGM,LLAA`<br>
<i>Example Two</i> `python wrapper.py --in_res=testing/sarah-test/1mhp_AHL.res --in_probs=testing/sarah-testing/1mhp_AHL_probs -f --in_seqs=testing/sarah-testing/seqs.txt`<br>

<i>Specifics</i>: Calculates STATIUM's binding score for a given sequence of amino acids in the positions listed in the input *.res file (see `create_res`). The `--in_probs` input is the STATIUM probabilities directory computed in `run_statium`. The presence of `-f` indicates that `--in_seqs` is a file (else just [possibly a set of] sequences, corresponding to the chains/position-pairs used to create the *.res file). For example, you might have a --in_seqs=AAA,L if your `--position_pairs` argument in `create_res` was 10-12,13. A file would contain similarly formatted argument, one sequence (set) on each line.

`--out` specifies an output file. If this is option is left out, results will be printed to the console. The presence of the z-score flags finds the z-scores of the input sequences' energy on a distribution of random sequences. 
***
<b>Helpful Hints</b>:
+ Verbose output is turned on by default. To turn verbose output off, include the '-nv' or '--noverbose' flag.
+ Arguments wrapped in parenthesis () are required; arguments wrapped in square brackets [] are optional.
+ `python wrapper.py -h` or `python wrapper.py --help` brings up an in-console summary of program arguments.
***
<b>Thanks for using STATIUM! Feel free to contact skoppula@mit.edu with issues.</b>:
<br>

			   statium_wrapper.py get_orig_seq (IN_PDB_ORIG) [-v | --verbose]
			   statium_wrapper.py [-f] generate_random_seqs (SEQ_LENGTH NUM_SEQS) [--OUT_FILE=None] [--TOTAL_PROTEIN_LIBRARY=None] [-v | --verbose]
			   statium_wrapper.py calc_top_seqs (IN_RES PROBS_DIR N) [OUT_FILE] [-v | --verbose]
			   statium_wrapper.py classify (RESULTS_FILE) [OUT_FILE] [ALPHA_THRESHOLD] [-v | --verbose]
			   statium_wrapper.py get_confusion_matrix (IN_RES CLASS_RESULTS TRUE_CLASS) [OUT_FILE] [--IN_PDB_ORIG=None] [-v | --verbose]
			   statium_wrapper.py [-i] calc_auroc (IN_RES RESULTS_FILE TRUE_CLASS) [OUT_FILE] [--IN_PDB_ORIG=None] [--CLASS_RESULTS=None] [-v | --verbose]
			   statium_wrapper.py [-i] plot_roc_curve (IN_RES RESULTS_FILE TRUE_CLASS) [--IN_PDB_ORIG=None] [--CLASS_RESULTS=None] [-v | --verbose]

	python statium_wrapper.py get_orig_seq testing/1ycr_mdm2_orig.pdb
	
	python statium_wrapper.py -v generate_random_seqs 11 10
	python statium_wrapper.py -fv generate_random_seqs 11 10 random_seqs.txt
	python statium_wrapper.py -fv generate_random_seqs 11 10 random_seqs.txt --TOTAL_PROTEIN_LIBRARY=data/all_protein_sequences.txt
	
	python statium_wrapper.py -v calc_top_seqs testing/1ycr_mdm2.res testing/1ycr_mdm2_output_probs/ 10
	
	
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

import this function and call with appropriate arguments to combine sequence energy and true classification into readable file
summarize(in_energy_file_path, in_true_class_file_path, in_res_path, in_pdb_orig, out_file = '1ycr_mdm2_summarize.txt')
	
