<b>STATIUM: smart scoring. promising proteins.</b><br>
STATIUM is an ongoing project at the Keating Lab to quantitatively understand how amino acid sequences interact. This repository implements a still-under-development algorithm we call STATIUM that scores how well two proteins bind at their interacting positions.

<b>Details and Documentation</b>
***
`precompute`:<br>
<i>Template</i>: `python wrapper.py precompute (--in_pdb --in_pdb_lib --in_ip_lib) [--out_dir]`<br>

<i>Example</i>:  `python wrapper.py precompute --in_pdb=1ycr_mdm2.pdb --in_pdb_lib=data/culled_90/ --in_ip_lib=data/ip_90_wGLY --out_dir=testing/output`<br>

<i>Specifics</i>: Runs STATIUM analysis to create weights for each possible amino acid at the protein's interacting pair positions. Note that this command is a shortcut combination of the `renumber`,`create_res`, and `run_statium` commands.

The easy-to-use `precompute` function assumes that the interacting sequence that you wish to analyze is the entirety of chain <b>B</b> (in the input PDB file), with SRN and SAN as 1 (see `renumber` for more information). If you wish to change these (or other parameters), run the above sequence of four commands with appropriate parameters.

The function takes as input the PDB of the to-be-analyzed complex (--in_pdb), the directory of the total protein library with all protein PDBs (--in_pdb_lib), and the directory containing a list of precomputed interacting pairs for each PDB in the protein library (--in_ip_lib).

***
`renumber`:<br>
<i>Template</i> `python wrapper.py renumber (--in_pdb) [--out_pdb --SRN --SAN --chains]`<br>

<i>Example</i> `python wrapper.py renumber --in_pdb=testing/1mph_HLA.pdb --out_pdb=testing/1mph_AHL.pdb --chains=H,L' <br>

<i>Specifics</i>: Takes a PDB file, strips away the meta-data, and renumbers the residues and atoms, retaining atom coordinate positions and removing the occupancy and temperature factors. Renumbering starts on the first valid line of the PDB file, at starting atom number = SAN and starting residue number = SRN. 'Valid line' is any PDB line with 'ATOM' or 'HETATM' with 'MSE' (selenomethionine).

One requirement of STATIUM as implemented here is that the input PDB has the receptor sequence listed before the peptide/mutant binder sequences in the file. The `renumber` function outputs a PDB that follows this requirement by reading the chains you wish to designate as the ligand sequence from the --chains option. In the example above, the heavy and light chain that we want to mutate were pushed after chain A (*HLA to *AHL).
***
`create_res`:<br>
<i>Template</i> `python wrapper.py create_res (--in_pdb_orig --in_pdb_renum) [--out_res --position_pairs]`<br>

<i>Example</i> `python wrapper.py create_res --in_pdb_orig=testing/1mph_AHL_orig.pdb --in_pdb_renum=testing/1mph_AHL_renum.pdb --out_res=testing/1mph_AHL.res` --position_pairs=L1,L20,H2,H33<br>

<i>Specifics</i>: Takes in both the original and renumbered PDB files (see 'renumber'). It translates pairs of (chain identifier, number) uniquely demarcating a residues on the original PDB file to a number uniquely demarcating a residue in the renumbered file. These set of numbers are written to a file and used as the positions to be analyzed by the STATIUM algorithm.

The --position_pairs argument specifies which the set of positions to be included as binder/ligand sequences in the STATIUM analysis. There must be an even number of comma seperated terms, which in groups of two, represent continuous sequence of residues to be included in the ligand sequence. In the example above, residues on the L chain, position 1-20, and residues on the H chain, positions 2-33 will be included in the output residues file.

If you fail to include a --position_pairs argument, the function will assume you mean to create a *.res file with the entirety of chain <i>B</i>.
***
`run_statium`:<br>
<i>Template</i> `python wrapper.py run_statium (--in_pdb --in_res --in_pdb_lib --in_ip_lib) [--out_dir --dist_cutoff]`<br>

<i>Example</i> `python wrapper.py run_statium --in_pdb=testing/1ycr_mdm2_renum.pdb --in_res=testing/1ycr_mdm2.res --in_pdb_lib=testing/culled90/ --in_ip_lib=testing/ip_90_wGLY/` <br>

<i>Specifics</i>: Takes in a renumbered PDB file (see 'renumber'), the directory of the total protein library with all protein PDBs (--in_pdb_lib), and the directory containing a list of precomputed interacting pairs for each PDB in the protein library (--in_ip_lib).

The function creates a directory containing a set of files, one file per interacting pair. A valid 'interacting pair' is a set of pair residues (one on the main chain and one on the sidechain [in STATIUM sidechain]) that are predicted to interact as per the criteria: 

The optional parameter --dist_cutoff defaults to 6 (Angstroms) and determines the pair distance cutoff at which 

Each file contains a set of probabilities, one for each amino acid identity describing how likely it is that that identity would exist at that position on the sidechain, given the fact that it is in an interacting pair with the main chain.
***
<b>Helpful Hints</b>:
+ Verbose output is turned on by default. To turn verbose output off, include the '-nv' or '--noverbose' flag.
+ Arguments wrapped in parenthesis () are required; arguments wrapped in square brackets [] are optional.
***
<b>Thanks for using STATIUM!</b>:
<br>

			   statium_wrapper.py [-f] calc_energy (IN_RES PROBS_DIR SEQ_OR_FILE) [OUT_FILE] [--IN_PDB_ORIG=None] [-z | --zscores] [-p | --percentiles] [-v | --verbose] [-d | --draw_histogram]
			   statium_wrapper.py get_orig_seq (IN_PDB_ORIG) [-v | --verbose]
			   statium_wrapper.py [-f] generate_random_seqs (SEQ_LENGTH NUM_SEQS) [--OUT_FILE=None] [--TOTAL_PROTEIN_LIBRARY=None] [-v | --verbose]
			   statium_wrapper.py calc_top_seqs (IN_RES PROBS_DIR N) [OUT_FILE] [-v | --verbose]
			   statium_wrapper.py classify (RESULTS_FILE) [OUT_FILE] [ALPHA_THRESHOLD] [-v | --verbose]
			   statium_wrapper.py get_confusion_matrix (IN_RES CLASS_RESULTS TRUE_CLASS) [OUT_FILE] [--IN_PDB_ORIG=None] [-v | --verbose]
			   statium_wrapper.py [-i] calc_auroc (IN_RES RESULTS_FILE TRUE_CLASS) [OUT_FILE] [--IN_PDB_ORIG=None] [--CLASS_RESULTS=None] [-v | --verbose]
			   statium_wrapper.py [-i] plot_roc_curve (IN_RES RESULTS_FILE TRUE_CLASS) [--IN_PDB_ORIG=None] [--CLASS_RESULTS=None] [-v | --verbose]
			   statium_wrapper.py [-h | --help]

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

import this function and call with appropriate arguments to combine sequence energy and true classification into readable file
summarize(in_energy_file_path, in_true_class_file_path, in_res_path, in_pdb_orig, out_file = '1ycr_mdm2_summarize.txt')
	
