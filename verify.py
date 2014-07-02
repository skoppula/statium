from util import filelines2deeplist

#based on a X/100 threshold, classifies top X% as strong binders, bottom X% as weak binders. Uses results file and last column of values in results file 
#SHOULD IDEALLY USE INVERSE NORM, NOW JUST HARDCODING IN THRESHOLD Z-SCORES FOR ALPHA=0.05
def classify(results_file, outfile, threshold=0.05):

	sorted_results = get_sorted_results(results_file)
	
	z_score_threshold_low = -1.645
	z_score_threshold_high = 1.645
	
	out = []
	for result in sorted_results:
		if(result[0] < z_score_threshold_low):
			out.append(result[1] + '\t' + str(result[0]) + '\t' + 'strong')
		elif(result[0] > z_score_threshold_high):
			out.append(result[1] + '\t' + str(result[0]) + '\t' + 'weak')
		else:
			out.append(result[1] + '\t' + str(result[0]) + '\t' + 'inconclusive')
			
	list2file(out, outfile)

def roc(scores_path, true_path, auroc_path, curve_path):
	#Read in scores data
	lines = filelines2deeplist(scores_path, skipComments=True, skipEmptyLines=True)
	scores = {pair[0]:pair[1] for pair in lines}

	#Read in true classification data
	lines = filelines2deeplist(true_path, skipComments=True, skipEmptyLines=True)
	true = {pair[0]:pair[1] for pair in lines}

	data = list()
	for seq, classification in true.items():
		if classification is 'weak':
			class_type = 0
		elif classification is 'strong':
			class_type = 1
		else:
			continue
		if seq in results:
			data.append((class_type, scores[seq]))
		else:
			print 'Error: ' + seq + ' not found in results.'

	import pyroc
	roc_data = pyroc.ROCData(data)

	if auroc_path is not None:
		auroc = roc_data.auc()
		list2file([auroc], auroc_path)

	if curve_path is not None:
		roc.plot_and_save(curve_path, 'STATIUM-based Binding Prediction Model')

#from true classification key (text file of format, SEQ\tCLASSIFICATION of 'weak' or 'strong'),
#returns dictionary
#helper function
def get_true_class(in_res_path, true_class_file, class_results=None, in_pdb_orig=None):
	
	true_class = filelines2deeplist(true_class_file, skipComments = True, useTuples = False, skipEmptyLines = True)
	true_class_dict = dict()
		
	#read back from .res file where Chain B starts
	lines = filelines2list(in_res_path)
	residue_positions = [int(line.strip()) - 1 for line in lines]
			
	for pair in true_class:
		seq = pair[0] if (len(pair) == 2) else pair[0]+' '+pair[1]
		seq = fix_sequence_line(seq, len(residue_positions), in_pdb_orig)
		true_class_dict[seq] = pair[1] if (len(pair) == 2) else pair[2]
		
	return true_class_dict

#helper function
def get_sorted_results(results_file):
	results = filelines2deeplist(results_file, skipComments=True, useTuples=True, skipEmptyLines=True)
	filtered_results = []
	for result in results:
		#remember seq stored in results file is already 'fixed' so we don't need to fix sequence it again
		filtered_results.append((float(result[-1]), result[0]))
	
	sorted_results = sorted(filtered_results)
	
	return sorted_results

#function to combine sequence energy and true classification files so that we can easily read file in format: SEQ	   ENERGY  TRUE_CLASS
def summarize(in_energy_file_path, in_true_class_file_path, in_res_path, in_pdb_orig, out_file = '1ycr_mdm2_summarize.txt'):
	lines = open(in_energy_file_path).readlines()
	lines2 = open(in_true_class_file_path).readlines()
	lines = [line for line in lines if line[0] != '#' and line[0] != '\n']
	lines2 = [line for line in lines2 if line[0] != '#' and line[0] != '\n' and line[0] != '\r']
	lines = [line[:-1] for line in lines]
	table1 = [line.split('\t') for line in lines]
	lines2 = [line[:-2] for line in lines2]
	table2 = [line.split('\t') for line in lines2 if not line.isspace()]

	res_lines = filelines2list(in_res_path)
	residue_positions = [int(line.strip()) - 1 for line in res_lines]
	table2 = [[fix_sequence_line(pair[0], len(residue_positions), in_pdb_orig),pair[1]] for pair in table2]

	out = []
	for pair in table1:
		for i, pair2 in enumerate(table2):
			if pair2[0] == pair[0]:
				s = pair[0] + ' ' + pair[1] + ' ' +  pair2[1]
				out.append(s)
				break

	out_parts = [(triplet.split(' ')[0], float(triplet.split(' ')[1]), triplet.split(' ')[2]) for triplet in out]
	out_parts = sorted(out_parts, key=itemgetter(1), reverse=True) 
	if out_file is not None:
		list2file([part[0] + ' ' + str(part[1]) + ' ' + part[2] for part in out_parts], out_file)
	return out_parts

