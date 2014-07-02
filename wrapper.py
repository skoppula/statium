import sys import ast
from docopt import docopt
from reformat import renumber
from reformat import create_res
from reformat import get_orig_seq
from analysis import preprocess
from analysis import statium
from analysis import calc_seq_energy
from analysis import generate_random_distribution
from analysis import calc_top_seqs
from analysis import classify
from analysis import get_confusion_matrix
from analysis import calc_auroc
from analysis import plot_roc_curve
from util import calc_seq_zscore
from util import generate_random_seqs
from util import list2file
from util import filelines2list

def main(argv):
	
	helpdoc =   	"""usage: wrapper.py renumber (--in_pdb=A) [--out_pdb=B --chains=C --SRN=1 --SAN=1] [--noverbose]
				wrapper.py create_res (--in_pdb_orig=A --in_pdb_renum=B) [--out_res=C --position_pairs=D] [--noverbose]
				wrapper.py preprocess (--in_dir=A) [--out_dir=B --ip_dist_cutoff=C] [--noverbose] [-r]
				wrapper.py run_statium (--in_res=A --in_pdb=B --pdb_lib=C) [--out=D --ip_dist_cutoff=E --matching_res_dist_cutoffs=F --counts] [--noverbose]
				wrapper.py [-f] energy (--in_res=A --in_probs=B --in_seqs=C) [--out=D] [-z | --zscore] [--histogram=E] [--noverbose]
				wrapper.py random (--seq_length=A --num_seqs=B) [--out=C] [--noverbose]
				wrapper.py get_orig_seq (--in_res=A --in_pdb_orig=B --in_pdb_renum=C) [--noverbose]
				wrapper.py calc_top_seqs (--in_res=A --probs_dir=B --N=C) [--out=D]

				wrapper.py calc_top_seqs (IN_RES PROBS_DIR N) [OUT_FILE] [--noverbose]
				wrapper.py classify (RESULTS_FILE) [OUT_FILE] [ALPHA_THRESHOLD] [--noverbose]
				wrapper.py get_confusion_matrix (IN_RES CLASS_RESULTS TRUE_CLASS) [OUT_FILE] [--IN_PDB_ORIG] [--noverbose]
				wrapper.py [-i] calc_auroc (IN_RES RESULTS_FILE TRUE_CLASS) [OUT_FILE] [--IN_PDB_ORIG] [--CLASS_RESULTS] [--noverbose]
				wrapper.py [-i] plot_roc_curve (IN_RES RESULTS_FILE TRUE_CLASS) [--IN_PDB_ORIG] [--CLASS_RESULTS] [--noverbose]
				wrapper.py [-h | --help]
			Options:

				--in_pdb=A	Input PDB file path
				--out_pdb=B	Output PDB file path
				--chains=C	Chosen ligand chains
				--SRN=1		Starting residue number
				--SAN=1		Starting atom number

				--in_pdb_orig=A	Input PDB file path (original)
				--in_pdb_renum=B	Input PDB file path (renumbered)
				--out_res=C	Output RES file path
				--position_pairs=D	Positions to include in the ligand

				--in_dir=A	Directory containing library PDBs
				--out_dir=B	Output directory for JSON objects
				--ip_dist_cutoff=C	Threshold for interacting pair designation

				--in_res=A	Input .res file path
				--in_pdb=B	Input renumbered PDB path
				--pdb_lib=C	Input preprocessed PDB library directory
				--out=D		Output directory
				--ip_dist_cutoff=E	Threshold for interacting pair designation
				--matching_res_dist_cutoffs=F	Thresholds for matching IP designation

				--in_res=A	Input .res file path
				--in_probs=B	Input STATIUM probabilities directory
				--in_seqs=C	Sequence patter or path to a file of sequence patterns to be scored
				--out=D		File path to output score (if -f flag is present)
				--histogram=E	File path to output histogram (absence outputs nothing)

				--seq_length=A	Length of the random sequences
				--num_seqs=B	Number of random sequences
				--out=C		Output file path

				--in_res=A	Input .res file path
				--in_pdb_orig=B	Input PDB file path (original)
				--in_pdb_renum-C	Input PDB file path (renumbered)


				--probs_dir=A	Input STATIUM probabilities directory
				--N		Number of sequences to be found
				--out		Output file path
			"""
	
	options = docopt(helpdoc, argv, help = True, version = "3.0.0", options_first=False)
	verbose = not options['--noverbose']
 
	if options['renumber']:
		in_pdb = options['--in_pdb']
		out_pdb = options['--out_pdb'] if options['--out_pdb'] is not None else in_pdb[:-4]+'_renumbered.pdb'
		SRN = 1 if options['--SRN'] == None else int(options['--SRN']) 
		SAN = 1 if options['--SAN'] == None else int(options['--SAN'])
		chains =  {'B'} if options['--chains'] == None else set(options['--chains'].split(','))

		if verbose: print "Renumbering PDB file: " + in_pdb
		renumber(SRN, SAN, chains, in_pdb, out_pdb)
		if verbose: print "Done. Renumbered file: " + out_pdb
	
		
	elif options['create_res']:
		pdb_orig = options['--in_pdb_orig']
		pdb_renum = options['--in_pdb_renum']
		res = pdb_orig[:-4]+'.res' if options['--out_res'] is None else options['--out_res']

		if not options['--position_pairs']:
			positions = {'B'}
		else:
			raw = options['--position_pairs'].split(',')
			positions = set()

			list_iter = iter(raw)
			for i, term in enumerate(list_iter):
				parts = term.split('-')
				if len(parts) == 1:
					positions.add(term)
					chain = term[0]
				elif len(parts) == 2:
					chain = parts[0][0]
					num = parts[0][1:]
					num2 = parts[1]
					positions.add((chain, num, num2))	
				else:
					sys.exit('Invalid position pairs')	

		if verbose: print "Creating .res file using: " + pdb_orig + " and " + pdb_renum
		create_res(pdb_orig, pdb_renum, res, positions)
		if verbose: print "Done. .res file: " + res

	elif options['preprocess']:
		in_dir = options['--in_dir']
		out_dir = options['--out_dir'] if options['--out_dir'] else in_dir + '_JSON_preprocessed'
		ip_dist = float(options['--ip_dist_cutoff']) if options['--ip_dist_cutoff'] is not None else 6.0
		restart = options['-r']

		if verbose: print 'Preprocessing library: %s' % in_dir
		preprocess(in_dir, out_dir, ip_dist, restart, verbose)
		if verbose: print 'Done: %s' % out_dir

	
	elif options['run_statium']:
		res = options['--in_res']
		pdb = options['--in_pdb']
		pdb_lib = options['--pdb_lib']
		out_dir = options['--out'] if options['--out'] is not None else res[:-4]
		ip_dist = float(options['--ip_dist_cutoff']) if options['--ip_dist_cutoff'] is not None else 6.0
		
		default = {'A':2, 'C':6, 'D':6, 'E':6, 'F':6, 'G':2, 'H':6, 'I':6, 'K':6, 'L':6, 'M':6, 'N':6, 'P':6, 'Q':6, 'R':6, 'S':6, 'T':6, 'V':6, 'W':6, 'Y':6, 'X':0}
		match_dist = ast.literal_eval(options['--matching_res_dist_cutoffs']) if options['--matching_res_dist_cutoffs'] else default
		count = True if options['--counts'] is not None else False 
		
		if verbose: print("\nRunning STATIUM with: " + pdb + " " + res + " " + pdb_lib+ " " + ip_lib)
		statium(res, pdb, pdb_lib, out_dir, ip_dist, match_dist, count, verbose)
		if verbose: print("Done. STATIUM probabilities in output directory: " + out_dir)

	elif options['energy']:
		zscores = options['-z'] or options['--zscore']
		histogram = options['--histogram']
 
		in_res = options['--in_res']
		probs_dir = options['--in_probs']
		isfile = options['-f']
		in_seqs = options['--in_seqs']
		out = options['--out']
		
		if verbose: print "Writing to file: " + isfile + ". Calculating z-score: " + zscores
		
		if zscores or percentiles:
			if verbose: print 'Generating random distribution of energies...'
			distribution = generate_random_distribution(in_res, probs_dir)
			if histogram:
				if verbose: print 'Drawing histogram...'
				import matplotlib.pyplot as plt
				import numpy as np

				hist, bins = np.histogram(distribution[2], bins=50)
				width = 0.7 * (bins[1] - bins[0])
				center = (bins[:-1] + bins[1:]) / 2
				plt.bar(center, hist, align='center', width=width)
				plt.savefig(histogram)

			if verbose: print 'Done generating random distribution.'
					
		if isfile:
			if verbose: print 'Reading sequences from file...'
			lines = filelines2list(in_seqs)
			out_lines = []
			
			for line in lines:
				if line != '' and line[0] != '#':
					seq = line.strip()
					energy = calc_seq_energy(in_res, probs_dir, seq)
					out = seq + "\t" + str(energy)
					
					if(zscores):
						out += "\t" + str(calc_seq_zscore(distribution[3], distribution[4], energy))
					
				out_lines.append(line)
					
			list2file(out_lines, outfile)
			print('Done.')
		
		else:
			energy = calc_seq_energy(in_res, probs_dir, in_seqs)
			print("Sequence energy for " + in_seqs + " is: " + str(energy))
			
			if(zscores):
				zscore = calc_seq_zscore(distribution[3], distribution[4], energy)
				print('Z-score is ' + str(zscore))

	#Get the original AA sequence of chain B, along with stats like the length and position of that chain
	elif options['get_orig_seq']:
		res = options['--in_res']
		orig_pdb = options['--in_pdb_orig']
		renum_pdb = options['--in_pdb_renum']
		get_orig_seq(res, orig_pdb, renum_pdb)
	
	elif options['generate_random_seqs']:

		num_seqs = int(options['--num_seqs'])
		seq_length = int(options['--seq_length'])
		out = options['--out']
		
		if verbose: print("Generating " + num_seqs + " random sequences of length " + seq_length)
		sequences = [generate_random_seq(seq_length) for _ in range(num_seqs)]
		
		if out is None:
			for sequence in sequences:
				print(sequence)
		else:
			list2file(sequences, out)
			if(verbose): print("Random sequences written to " + out)
   
	elif options['calc_top_seqs']:
		in_res = options['--in_res']
		probs_dir = options['--probs_dir']
		N = int(options['--N'])
i
		if verbose: print 'Calculating ' + N + ' sequences with lowest energy.'
		results = calc_top_seqs(in_res, probs_dir, N)
		if probs_dir:
			out = [seq + '\t' + str(energy) + '\n' for seq, energy in results]
			out_path = 'top_' + str(N) + '_sequences.txt' if probs_dir == None else probs_dir
			list2file(out, out_path)
		else:
			print sequences
	
	elif(options['classify']):
		threshold = 0.05 if options['ALPHA_THRESHOLD'] == None else float(options['ALPHA_THRESHOLD'].split('=')[1])
		if(verbose): print('Classifying ' + options['RESULTS_FILE'] + ' with dummy (i.e. not considered right now in analysis) threshold ' + str(threshold))
		outfile = 'classify_results_' + options['RESULTS_FILE'] + '_' + str(options['N']) + '.txt' if (options['OUT_FILE'] == None) else options['OUT_FILE']
		classify(options['RESULTS_FILE'], outfile, threshold)
		if(verbose): print('Done. Results in ' + outfile)
	
	elif(options['get_confusion_matrix']):
		if(verbose): print('Calculating confusion matrix for ' + options['CLASS_RESULTS'] + ' with true classifications in ' + options['TRUE_CLASS'])
		(TP, FP, TN, FN) = get_confusion_matrix(options['IN_RES'], options['CLASS_RESULTS'], options['TRUE_CLASS'], options['--IN_PDB_ORIG'])
		outfile = options['CLASS_RESULTS'] + '_confusion_matrix.txt' if (options['OUT_FILE'] == None) else options['OUT_FILE']
		out_str = 'TP: ' + str(TP) + '\t FN: ' + str(FN) + '\nFP: ' + str(FP) + '\tTN: ' + str(TN)
		out = [out_str]
		list2file(out, outfile)
		print(out_str)
		if(verbose): print('Confusion matrix written out to ' + outfile)
	
	#by default, analysis includes tentative inconclusives. the -i flag removes them from analysis 
	elif(options['calc_auroc']):
		if(verbose): print('Calculating AUROC for ' + options['RESULTS_FILE'] + ' with true classifications in ' + options['TRUE_CLASS'])
		
		class_results = options['--CLASS_RESULTS']
		if(verbose):
			if(options['-i']):
				print('Discarding tentative \'inconclusive sequences\' in ROC analysis')
			else:
				print('Including tentative \'inconclusive sequences\' in ROC analysis')
		
			  
		auroc = calc_auroc(options['IN_RES'], options['RESULTS_FILE'], options['TRUE_CLASS'], class_results, options['--IN_PDB_ORIG'])
		outfile = options['CLASS_RESULTS'] + '_auroc.txt' if (options['OUT_FILE'] == None) else options['OUT_FILE']
		list2file([str(auroc)], outfile)
		print(auroc)
		if(verbose): print('AUROC written out to ' + outfile)
		
	elif(options['plot_roc_curve']):
		if(verbose): print('Plotting ROC curve for ' + options['RESULTS_FILE'] + ' with true classifications in ' + options['TRUE_CLASS'])
		if(verbose):
			if(options['-i']):
				print('Discarding tentative \'inconclusive sequences\' in ROC analysis')
				class_results = options['--CLASS_RESULTS']
			else:
				print('Including tentative \'inconclusive sequences\' in ROC analysis')
				class_results = None
				
		plot_roc_curve(options['IN_RES'], options['RESULTS_FILE'], options['TRUE_CLASS'], class_results, options['--IN_PDB_ORIG'])
		if(verbose): print('Done.')
	
if __name__ == "__main__":
	main(sys.argv[1:])
