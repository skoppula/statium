import pstats
import cProfile
import sys
import ast
from docopt import docopt
from reformat import renumber
from reformat import create_res
from reformat import get_orig_seq
from analysis import statium
from analysis import calc_seq_energy
from analysis import generate_random_distribution
from analysis import calc_top_seqs
from util import calc_seq_zscore
from util import generate_random_seq
from util import list2file
from util import filelines2list
from verify import roc
from verify import print_merged

def parse_position_pairs(in_str):
	if not in_str:
		positions = {'B'}
	else:
		raw = in_str.split(',')
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
	return positions

def main(argv):
	
	helpdoc =   	"""usage: wrapper.py quickrun (--in_pdb=A --position_pairs=B --pdb_lib=C --ip_lib=D) [--out=E] [--noverbose]
				wrapper.py renumber (--in_pdb=A) [--out_pdb=B --chains=C --SRN=1 --SAN=1] [--noverbose]
				wrapper.py create_res (--in_pdb_orig=A --in_pdb_renum=B) [--out_res=C --position_pairs=D] [--noverbose]
				wrapper.py preprocess (--in_dir=A) [--out_dir=B --ip_dist_cutoff=C] [--noverbose] [-r]
				wrapper.py run_statium (--in_res=A --in_pdb=B --pdb_lib=C) [--ip_lib=D --out=E --ip_dist_cutoff=F --matching_res_dist_cutoffs=G --backbone --filter --counts] [--noverbose]
				wrapper.py [-f] energy (--in_res=A | --in_pdb=B) (--in_probs=C --in_seqs=D) [--out=E] [-z | --zscore] [--histogram=E] [--noverbose]
				wrapper.py random (--seq_length=A --num_seqs=B) [--out=C] [--noverbose]
				wrapper.py get_orig_seq (--in_res=A --in_pdb_orig=B --in_pdb_renum=C) [--noverbose]
				wrapper.py calc_top_seqs (--in_res=A --in_probs=B --N=C) [--out=D] [--noverbose]
				wrapper.py roc (--scores=A --true=B) [--curve=C --auroc=D] [-noverbose]
				wrapper.py print_merged (--scores=A --true=B) [--out=C] [--noverbose]
				wrapper.py [-h | --help]
			Options:

				--in_pdb=A	Input PDB file path
				--position_pairs=B	Positions to include in the binding sequence
				--pdb_lib=C	Input directory of library PDB files
				--ip_lib=D	Input directory of library IP files
				--out=E	Output directory

				--in_pdb=A	Input PDB file path
				--out_pdb=B	Output PDB file path
				--chains=C	Chosen ligand chains
				--SRN=1		Starting residue number
				--SAN=1		Starting atom number

				--in_pdb_orig=A	Input PDB file path (original)
				--in_pdb_renum=B	Input PDB file path (renumbered)
				--out_res=C	Output RES file path
				--position_pairs=D	Positions to include in the binding sequence

				--in_dir=A	Directory containing library PDBs
				--out_dir=B	Output directory for JSON objects
				--ip_dist_cutoff=C	Threshold for interacting pair designation

				--in_res=A	Input .res file path
				--in_pdb=B	Input renumbered PDB path
				--pdb_lib=C	Input directory of library PDB files
				--ip_lib=D	Input directory of library IP files
				--out=E		Output directory
				--ip_dist_cutoff=F	Threshold for interacting pair designation
				--matching_res_dist_cutoffs=G	Thresholds for matching IP designation

				--in_res=A	Input .res file path
				--in_pdb=B	Input .pdb file path
				--in_probs=C	STATIUM probabilities file
				--in_seqs=D	Sequence patter or path to a file of sequence patterns to be scored
				--out=E		File path to output score (if -f flag is present)
				--histogram=F	File path to output histogram (absence outputs nothing)

				--seq_length=A	Length of the random sequences
				--num_seqs=B	Number of random sequences
				--out=C		Output file path

				--in_res=A	Input .res file path
				--in_pdb_orig=B	Input PDB file path (original)
				--in_pdb_renum-C	Input PDB file path (renumbered)

				--in_res=A
				--in_probs=B	STATIUM output file
				--N=C		Number of sequences to be found
				--out=D		Output file path

				--scores=A	Sequences w/ energy file path
				--true=B	Sequences' true binding classification file path
				--auroc=C	File path to output auroc
				--curve=D	File path to output ROC curve

				--scores=A	Sequences w/ energy file path
				--true=B	Sequences' true binding classification file path
				--out=C		Output file path
			"""
	
	options = docopt(helpdoc, argv, help = True, version = "3.0.0", options_first=False)
	verbose = not options['--noverbose']

	if options['quickrun']:
		in_pdb = options['--in_pdb']
		stem = in_pdb[:-4]
		renum_pdb = stem + '_renumbered.pdb'
		res = stem + '.res'

		pdb_lib = options['--pdb_lib']
		ip_lib = options['--ip_lib']
		out_dir = options['--out'] if options['--out'] is not None else stem

		positions = parse_position_pairs(options['--position_pairs'])
		chains = [term[0] for term in positions]
		default_match_dist = {'A':2, 'C':6, 'D':6, 'E':6, 'F':6, 'G':2, 'H':6, 'I':6, 'K':6, 'L':6, 'M':6, 'N':6, 'P':6, 'Q':6, 'R':6, 'S':6, 'T':6, 'V':6, 'W':6, 'Y':6, 'X':0}
		ip_dist = 6.0

		if verbose: print "Renumbering PDB file: " + in_pdb
		renumber(1, 1, chains, in_pdb, renum_pdb)
		if verbose: print "Creating .res file using: " + in_pdb + " and " + renum_pdb
		create_res(in_pdb, renum_pdb, res, positions)
		if verbose: print "Running STATIUM with: " + renum_pdb + " " + res + " " + pdb_lib + ' and IP lib: ' + ip_lib
		statium(res, renum_pdb, pdb_lib, ip_lib, out_dir, ip_dist, default_match_dist, False, False, False, verbose)
		if verbose: print 'Done'


	elif options['renumber']:
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
		position_pairs = options['--position_pairs'] if options['--position_pairs'] else 'B'


		if verbose: print "Creating .res file using: " + pdb_orig + " and " + pdb_renum
		create_res(pdb_orig, pdb_renum, res, position_pairs)
		if verbose: print "Done. .res file: " + res

	elif options['preprocess']:
		in_dir = options['--in_dir']
		out_dir = options['--out_dir'] if options['--out_dir'] else in_dir + '_JSON_preprocessed'
		ip_dist = float(options['--ip_dist_cutoff']) if options['--ip_dist_cutoff'] is not None else 5.0
		restart = options['-r']

		if verbose: print 'Preprocessing library: %s' % in_dir
		preprocess(in_dir, out_dir, ip_dist, restart, verbose)
		if verbose: print 'Done: %s' % out_dir

	
	elif options['run_statium']:
		res = options['--in_res']
		pdb = options['--in_pdb']
		pdb_lib = options['--pdb_lib']
		ip_lib = options['--ip_lib']
		out = options['--out'] if options['--out'] is not None else res[:-4] + '.out'
		ip_dist = float(options['--ip_dist_cutoff']) if options['--ip_dist_cutoff'] is not None else 6.0
		
 		default = {'A':0.2, 'C':0.4, 'D':0.4, 'E':0.4, 'F':0.4, 'G':0.2, 'H':0.4, 'I':0.4, 'K':0.4, 'L':0.4, 'M':0.4, 'N':0.4, 'P':0.4, 'Q':0.4, 'R':0.4, 'S':0.4, 'T':0.4, 'V':0.4, 'W':0.4, 'Y':0.4}
		match_dist = ast.literal_eval(options['--matching_res_dist_cutoffs']) if options['--matching_res_dist_cutoffs'] else default
		backbone = options['--backbone']
		filter_sidechain = options['--filter']
		count = options['--counts']
		
		if verbose: print "\nRunning STATIUM with: " + pdb + " " + res + " " + pdb_lib + ' and IP lib: ' + str(ip_lib)
		statium(res, pdb, pdb_lib, ip_lib, out, ip_dist, match_dist, backbone, filter_sidechain, count, verbose)
		if verbose: print "Done. STATIUM probabilities in output directory: " + out_dir

	elif options['energy']:
		zscores = options['-z'] or options['--zscore']
		histogram = options['--histogram']
 
		in_res = options['--in_res']
		in_probs = options['--in_probs']
		isfile = options['-f']
		in_seqs = options['--in_seqs']
		outfile = options['--out']
		
		if verbose: print "Writing to file: " + str(isfile) + ". Calculating z-score: " + str(zscores)
		
		if zscores:
			if verbose: print 'Generating random distribution of energies...'
			distribution = generate_random_distribution(in_res, in_probs)
			if histogram:
				if verbose: print 'Drawing histogram...'
				import matplotlib.pyplot as plt
				import numpy as np

				hist, bins = np.histogram(distribution[1], bins=50)
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
					energy = calc_seq_energy(in_res, in_probs, seq)
					out = seq + "\t" + str(energy)
					
					if(zscores):
						out += "\t" + str(calc_seq_zscore(distribution[2], distribution[3], energy))
					
				out_lines.append(out)
					
			list2file(out_lines, outfile)
			print('Done.')
		
		else:
			energy = calc_seq_energy(in_res, in_probs, in_seqs)
			print("Sequence energy for " + in_seqs + " is: " + str(energy))
			
			if(zscores):
				zscore = calc_seq_zscore(distribution[2], distribution[3], energy)
				print('Z-score is ' + str(zscore))

	#Get the original AA sequence of chain B, along with stats like the length and position of that chain
	elif options['get_orig_seq']:
		res = options['--in_res']
		orig_pdb = options['--in_pdb_orig']
		renum_pdb = options['--in_pdb_renum']
		get_orig_seq(res, orig_pdb, renum_pdb)
	
	elif options['random']:
		num_seqs = int(options['--num_seqs'])
		seq_length = int(options['--seq_length'])
		out = options['--out']
		
		if verbose: print("Generating " + str(num_seqs) + " random sequences of length " + str(seq_length))
		sequences = [generate_random_seq(seq_length) for _ in range(num_seqs)]
		
		if out is None:
			for sequence in sequences:
				print(sequence)
		else:
			list2file(sequences, out)
			if(verbose): print("Random sequences written to " + out)
   
	elif options['calc_top_seqs']:
		in_res = options['--in_res']
		probs_dir = options['--in_probs']
		N = int(options['--N'])
		out = options['--out']

		if verbose: print 'Calculating ' + str(N) + ' sequences with lowest energy.'
		results = calc_top_seqs(in_res, in_probs, N)
		if out:
			to_print = [seq + '\t' + str(energy) + '\n' for seq, energy in results]
			list2file(to_print, out)
		else:
			print results

	elif options['roc']:
		scores = options['--scores']
		true = options['--true']
		auroc = options['--auroc']
		curve = options['--curve']
		if verbose: print 'Calculating AUROC for ' + scores + ' with true classifications in ' + true
		roc(scores, true, auroc, curve)
		if verbose: print 'Done. Printed AUROC to: ' + str(auroc) + ' and ROC curve to ' + str(curve)

	elif options['print_merged']:
		scores = options['--scores']
		true = options['--true']
		out = options['--out'] if options['--out'] is not None else 'merged.txt'
		print_merged(scores, true, out)
		print 'Done. Printed to ' + out
	
if __name__ == "__main__":
        timing_path = 'timing_analysis.txt'
	cProfile.run('main(sys.argv[1:])', timing_path)
        

