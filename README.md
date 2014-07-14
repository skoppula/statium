<b>STATIUM: smart scoring. promising proteins.</b><br>
STATIUM is an ongoing project at the Keating Lab to quantitatively understand how amino acid sequences interact. This repository contains a friendly implemtation of the structure-based statistical-potential STATIUM algorithm that scores how well two or more proteins bind at their interacting positions.

<b>Installation</b><br>
If you are on a *nix machine with `git` installed obtaining STATIUM and all its data should be as simple as: `git clone https://github.com/skoppula/statium.git`. Then, to extract the library data: interacting-pair data can be extracted by `tar -zxvf data/ip_90_wGLY.tar.gz`. To extract and recombine the library PDB files, you can run `mkdir culled_90` followed by `i=0; for i in {0..9}; do tar -zxf culled_90_$i.tar.gz; mv culled_90_$i/* culled_90/; done`.
***
<b>Quickstart! (Example Workflow)</b><br>
If you, for example, wanted to score sequences for chain B of some protein described in 1YCR.pdb, you could run:
python renumber --in_pdb=1YCR.pdb --out_pdb=1YCR_renumbered.pdb --chains=B
python create_res --in_pdb_orig=1YCR.pdb --in_pdb_renum=1YCR_renumbered.pdb --out_res=1YCR.res --position_pairs=B
python run_statium --in_pdb=1YCR_renumbered.pdb --in_res=1YCR.res --pdb_lib=culled_90/ --ip_lib=ip_90_wGLY/ --out_dir=1YCR/
python energy --in_res=1YCR.res --probs_dir=1YCR --in_seqs=AMLTGTMMXX<br>

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

<i>Example</i> `python wrapper.py run_statium --in_pdb=testing/1mhp_AHL_new.pdb --in_res=testing/1mhp_AHL.res --pdb_lib=data/culled_90/ --ip_lib=data/ip_90_wGLY/` <br>

<i>Specifics</i>: Takes in a renumbered PDB file (see `renumber`), the directory of the preprocessed protein library containing each PDB's *.pickle file (--pdb_lib).

Optional parameters include: --out_dir (the directory where STATIUM outputs its results; default value is value of --in_pdb without the .pdb extension), --counts (whether to print out STATIUM's intermediate analysis outputs; note that this takes no argument; simply including the flag issues printing!), --ip_dist_cutoff (the threshold distance in Angstroms between two atoms, below which the atom's residues are deemed 'interacting'; default is 6.0), and --matching_res_dist_cutoff (a dictionary with all twenty amino acids [in single character representation] each mapped to a cutoff below which a interacting residue pair cannot be deemed 'matching' to a library protein interacting pair. Example of using this parameter [containing the default, recommended dictionary values if you leave this parameter out]: --matching_res_dist_cutoff={'A':2, 'C':6, 'D':6, 'E':6, 'F':6, 'G':2, 'H':6, 'I':6, 'K':6, 'L':6, 'M':6, 'N':6, 'P':6, 'Q':6, 'R':6, 'S':6, 'T':6, 'V':6, 'W':6, 'Y':6, 'X':0}.

The function creates a directory containing a set of files, one file per interacting pair:

Each file contains a set of twenty probabilities (one for each amino acid) describing how likely it is for that identity would exist at that position on the sidechain, given the main chain's amino acid identity at the position.<br>

<i>Dependencies</i>: `fsolve` from `scipy.optimize` if there are glycine residues in any of the interacting pair positions
***
`random`:<br>
<i>Template</i> `python wrapper.py random (--seq_length --num_seqs) [--out]`<br>
<i>Example</i>`python wrapper.py random --seq_length=8 --num_seqs=10 --out=testing/random-sequences.txt`<br>

<i>Specifics</i>: Generates `--num_seqs` random sequences of '--seq_length' length. If you include a `--out=X` option, the random sequences will be printed to the specified file. Sequences are randomly drawn from the collection of all known protein sequences contained in `data/all_protein_sequences.txt'. If you choose to modify this (e.g. adjust it so that certain amino acids occur with certain frequencies, ensure that only amino acid in their character representation are present).
***
`get_orig_seq`:<br>
<i>Template</i> `python wrapper.py get_orig_seq (--in_res --in_pdb_orig --in_pdb_renum)`<br>
<i>Example</i> `python wrapper.py get_orig_seq ---in_res=testing/1mph_AHL.res -in_pdb_orig=testing/1mph_AHL_orig.pdb --in_pdb_renum=testing/1mph_AHL_renum.pdb 

Reverse of the `renumber` function. From *.res file and the renumbered and original PDBs (see `renumber`) outputs the list of residues with original chain and position information.
***
`energy`:<br>
<i>Template</i> `python wrapper.py energy (--in_res --in_probs) [-f] (--in_seqs) [--out] [-z | --zscores] [--histogram]`<br>
<i>Example One</i> `python wrapper.py --in_res=testing/1mhp_AHL.res --in_probs=testing/1mhp_AHL_probs --in_seqs=AAAGGGM,LLAA -z --histogram='hist.jpg'`<br>
<i>Example Two</i> `python wrapper.py --in_res=testing/1mhp_AHL.res --in_probs=testing/1mhp_AHL_probs -f --in_seqs=testing/seqs.txt`<br>

<i>Specifics</i>: Calculates STATIUM's binding score for a given sequence of amino acids in the positions listed in the input *.res file (see `create_res`). The `--in_probs` input is the STATIUM probabilities directory computed in `run_statium`. The presence of `-f` indicates that `--in_seqs` is a file (else just [possibly a set of] sequences, corresponding to the chains/position-pairs used to create the *.res file). For example, you might have a --in_seqs=AAA,L if your `--position_pairs` argument in `create_res` was 10-12,13 (note that an in_seqs without a comma is also acceptable: e.g. AAAL). A file would contain similarly formatted argument, one sequence (set) on each line.

`--out` specifies an output file. If this is option is left out, results will be printed to the console. The presence of the z-score flags finds the z-scores of the input sequences' energy on a distribution of random sequences. The presence of the `--histogram=X` saves a histogram of the random distribution of scores generated for use in z-score and percentile calculations.

If you wish to adjust the number of random sequences in the distribution, you can modify the `generate_random_distribution` function. Currently, 100,000 random-sequence scores are used to create the distribution.

<i>Dependencies</i>: `matplotlib.pyplot` and `numpy` in order to use `--histogram`
***
`calc_top_seqs`<br>
<i>Template</i> `python wrapper.py calc_top_seqs (--in_res --probs_dir --N) [--out]`<br>
<i>Example</i> `python wrapper.py calc_top_seqs --probs_dir=testing/1mph_AHL_probs --out=testing/top_100_seqs.txt --N=100`

<i>Specifics</i>: Calculates the top `N` sequences with the lowest STATIUM energies (best predicted binders). `--probs_dir` is the output of the `run-statium` function and `--in_res` the *.res file produced by `create_res`.
***
`roc`<br>
<i>Template</i> `python wrapper.py roc (--scores --true) [--curve --auroc]`<br>
<i>Example</i> `python wrapper.py roc --scores=testing/energies.txt --true=testing/true-classification.txt --auroc=testing/auroc.txt --curve=testing/roc-curve.png`<br>

<i>Specifics</i>: Given a file (output of `energy`) with containing ligand sequences and their energies (`--scores`), a file with those sequences' true binding classification as strong, weak, or inconclusive (`--true`), outputs the corresponding ROC curve (if `--curve` is listed with an output file path) and/or the AUROC (if `--auroc` is listed with an output file path).<br>

<i>Dependencies</i>: `matplotlib` for plotting ROC curves
***
`print_merged`<br>
<i>Template</i> `python wrapper.py print_merged (--scores --true) [--out]'<br>
<i>Example</i> 'python wrapper.py print_merged --scores=testing/energies.txt --true=testing/true-classifcation.txt --out=testing/scores-and-true.txt`<br>

<i>Specifics</i>: Combines STATIUM scores and true binding classification files into one output file, with each line containing a sequence, its score, and true classification. Sequences in the scores file but not in the classification file will not appear in the file. Conversely, scores in the classification but not in the scores file will be listed as with score infinity.
***
<b>Helpful Hints</b>:
+ Verbose output is turned on by default. To turn verbose output off, include the '-nv' or '--noverbose' flag.
+ Arguments wrapped in parenthesis () are required; arguments wrapped in square brackets [] are optional.
+ `python wrapper.py -h` or `python wrapper.py --help` brings up an in-console summary of program arguments. <br>
+ Make sure that all sequences you query are in the right frame as the residue positions specified in the *.res file! You can insert X's in the beginning of a sequence to shift it to the correct frame. <br>
+ The PDB libraries are named `culled_90` because they are a curated set of PDBs with at maximum ninety-percent sequence similarity.

<br>
<b>Thanks for using STATIUM! Feel free to contact skoppula@mit.edu with issues.</b>
<br>
