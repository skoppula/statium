from util import filelines2deeplist

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
		if seq in scores:
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

def print_merged(score_path, true_path, out):
	#Read in scores data
	lines = filelines2deeplist(scores_path, skipComments=True, skipEmptyLines=True)
	scores = {pair[0]:pair[1] for pair in lines}

	#Read in true classification data
	lines = filelines2deeplist(true_path, skipComments=True, skipEmptyLines=True)
	true = {pair[0]:pair[1] for pair in lines}

	data = list()
	for seq, classification in true.items():
		if seq in scores:
			data.append((seq, scores[seq], classification))
		else:
			data.append((seq, float('inf'), classification))

	data.sort(key=lambda t: t[1])
	
	with open(out, 'w') as f:
		for datum in data:
			out.write(datum[0] + '\t' + datum[1] + '\t' + datum[2]+ '\n')	
