statium!
=======

Implementation of the STATIUM algorithm [protein-protein binding affinity scoring]

		usage: statium_wrapper.py renumber (IN_PDB) [OUT_PDB] [-v | --verbose]
			   statium_wrapper.py create_res (IN_PDB_ORIG IN_PDB_RENUMBERED) [OUT_RES] [-v | --verbose]
			   statium_wrapper.py run_statium (IN_RES IN_PDB IN_PDB_LIB_DIR IN_IP_LIB_DIR) [OUT_DIR] [-v | --verbose]
			   statium_wrapper.py [-f] calc_energy (IN_RES PROBS_DIR SEQ_OR_FILE) [OUT_FILE] [--IN_PDB_ORIG=None] [-z | --zscores] [-p | --percentiles] [-v | --verbose]
			   statium_wrapper.py get_orig_seq (IN_PDB_ORIG) [-v | --verbose]
			   statium_wrapper.py [-f] generate_random_seqs (SEQ_LENGTH NUM_SEQS) [--OUT_FILE=None] [--TOTAL_PROTEIN_LIBRARY=None] [-v | --verbose]                       
			   statium_wrapper.py calc_top_seqs (IN_RES PROBS_DIR N) [OUT_FILE MAX_TIME]
			   statium_wrapper.py [-h | --help]

An example of an analysis sequence: <br>

	python statium_wrapper.py renumber testing/4hfz_orig.pdb testing/4hfz.pdb -v
	python statium_wrapper.py create_res testing/4hfz_orig.pdb testing/4hfz.pdb testing/4hfz.res -v
	python statium_wrapper.py run_statium testing/4hfz.res testing/4hfz.pdb data/culled_90/ data/ip_90_wGLY/ testing/output -v
	python statium_wrapper.py get_orig_seq testing/4hfz_orig.pdb
	
	python statium_wrapper.py -v generate_random_seqs 11 10
	python statium_wrapper.py -fv generate_random_seqs 11 10 random_seqs.txt
	python statium_wrapper.py -fv generate_random_seqs 11 10 random_seqs.txt --TOTAL_PROTEIN_LIBRARY=data/all_protein_sequences.txt
	
	python statium_wrapper.py calc_energy testing/4hfz.res testing/4hfz_output_probs/ ETFSDLWKLLP
	python statium_wrapper.py calc_energy -f testing/4hfz.res testing/4hfz_output_probs/ testing/4hfz_sequences.txt testing/4hfz_seq_energies.txt
	
	#if sequences in file are of irregular length (different from original chain length), include start position in front of irregular sequences (e.g. '17 EALL') and use parameters like follows
	python statium_wrapper.py calc_energy -fv testing/4hfz.res testing/4hfz_output_probs/ testing/4hfz_sequences.txt testing/4hfz_seq_energies.txt --IN_PDB_ORIG=testing/4hfz_orig.pdb
	
	#outputs energy and corresponding zscores as compared to a distribution of random sequences
	python statium_wrapper.py calc_energy -fzv testing/4hfz.res testing/4hfz_output_probs/ testing/4hfz_sequences.txt testing/4hfz_seq_zscores.txt --IN_PDB_ORIG=testing/4hfz_orig.pdb
	#outputs energy and corresponding zscores and percentiles as compared to a energy distribution of random sequences
	python statium_wrapper.py calc_energy -fzpv testing/4hfz.res testing/4hfz_output_probs/ testing/4hfz_sequences.txt testing/4hfz_seq_zscores_percentiles.txt --IN_PDB_ORIG=testing/4hfz_orig.pdb
	
	python statium_wrapper.py calc_top_seqs testing/4hfz.res testing/4hfz_output_probs/ 10