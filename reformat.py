from util import AA2char
from util import get_random_AA
import re
import random

#Function to create *.res file:
#	each line a residue position (in renumbered file) for STATIUM analysis
def create_res(pdb_orig_path, pdb_renum_path, out_res_path, start=None, end=None, chain='B'):
	
	orig =  open(pdb_orig_path, 'r')
	renum =  open(pdb_renum_path, 'r')
	out =  open(out_res_path, 'w')
	chain = 'B' if chain is None else chain
	
	lines_orig = orig.readlines()
	lines_renum = renum.readlines()
	
	chain_start_line = ""
	num_residues = 0
	init = 0
	
	#count the number of residues
	#	by identifying chain's start line (or line with start=residue number)
	#	and end line
	for line in lines_orig:
		if(line[21]==chain and chain_start_line == "" and line[0:4] =='ATOM'):
			if (not start) or int(line[22:27]) == start:
				chain_start_line = line

		elif(line[21]==chain and chain_start_line != "" and (line[0:3] =='TER' or int(line[22:27]) == end)):
			num_residues = int(line[22:27]) - int(chain_start_line[22:27]) + 1
			break

	#Get start position in renumbered file
	#	If AA identity and coordinates are the same, must be same atom and residue
	for line in lines_renum:
		if(line[17:21] == chain_start_line[17:21] and line[32:55] == chain_start_line[32:55]):
			init = int(line[22:27])
			break
	
	if(num_residues < 1):
		print('Could not find a valid' + chain + '-chain.')
		return
	
	#count from start position in renumbered file
	#	to start position + num_residues	
	for i in range(num_residues-1):
		out.write(str(i+init) + '\n')
		
	outfile.write(str(num_residues-1+init)) #to ensure no newline at end of file
	
	out.close()
	orig.close()
	renum.close()

#Take a .pdb file from pdb.org and strip away meta-data so that output PDB only contains atoms and coordinates
def renumber(start_res_num, start_atom_num, in_pdb, out_pdb):
   
	infile =  open(in_pdb_path, 'r')
	outfile =  open(out_pdb_path, 'w')
	lines = infile.readlines()
	
	RN = start_res_num #initialize residue number
	AN = start_atom_num #initialize atom number
	
	FRN = -1; 
	
	for line in lines:
		line = line.strip()
		#add whitespace so checks below don't error out
		line = (line + ' '*16 + '\n') if len(line) < 61 else (line + '\n') 
		#check if the line corresponds to an ATOM
		# or coordinated-molecule (specifically, MSE = selenomethionine)
		if(line[0:4] == 'ATOM' or (line[0:6] == 'HETATM' and line[17:20] == 'MSE')):
			
			#if first line, set FRN to first residue number
			FRN = int(line[23:26]) if FRN<0 else FRN  
			CRN = int(line[23:26])	#current residue number
			
			#checks if next residue has been reached
			#if so, update FRN and RN
			if CRN != FRN: 
				RN += 1
				FRN = CRN
			
			#Replace old number with new residue count
			num_digits = len(str(RN))
			line = line[:21] + ' '*(5-num_digits) + str(RN) + line[26:]
			
			#Replace old number with new atom count
			num_digits = len(str(AN))
			line = line[:6] + ' '*(5-num_digits) + str(AN) + line[11:]
			AN += 1
			
			#replace selenomethionine with MET
			if(line[17:20] == 'MSE'):
				line = line[:18] + 'ET' + line[20:]
				line = 'ATOM  ' + line[6:]
		
			#replace occupancy and temperature factor
			line = line[:56] + '1.00  0.00' + line[66:]
		
			outfile.write(line)
		
	outfile.write('TER\nEND')
	infile.close()
	outfile.close()

#Get the original AA sequence of chain B, along with stats like the length and position of that chain
def get_orig_seq(in_pdb_path_orig):
	
	infile =  open(in_pdb_path_orig, 'r')
	lines = infile.readlines()
	
	bchain_started = False
	start, end = 0, 0
	num_residues = 0
	residues = []
	sequence = ''
	
	for line in lines:
		line = line.strip()
		line = (line + ' '*16 + '\n') if (len(line) < 61) else (line + '\n')
	 
		if(line[0:4] == 'ATOM' or (line[0:6] == 'HETATM' and line[17:20] == 'MSE') or line[0:3] == 'TER'):
			
			if(line[21]=='B' and not bchain_started):
				bchain_started = True
				start = int(line[22:27])
				residues.append(start)
				sequence += AA2char(line[17:20])
			
			elif(line[0:3] == 'TER' and bchain_started):
				end = int(line[22:27])
				num_residues = end - start + 1
				return (sequence, num_residues, start, end)

			elif(line[21]=='B' and (int(line[22:27]) not in residues)):
				residues.append(int(line[22:27]))
				sequence += AA2char(line[17:20])

def generate_random_seqs(seq_length, num_seqs, library_path = 'data/all_protein_sequences.txt'):
	
	use_protein_library = False if(library_path == None) else True
	
	size = 0
	if(use_protein_library):
		try:
			with open (library_path, "r") as myfile:
				data = myfile.read().replace('\n', '')
				pattern = re.compile(r'\s+')
				data = re.sub(pattern, '', data)
				size = len(data)
		except IOError:
			print('Could not find file ' + library_path + '. Using random AA generation')
			use_protein_library = False
		
	sequences = []
	for _ in range(num_seqs):
		if(use_protein_library):
			start = random.randint(0, size - seq_length)
			sequence = data[start:(start+seq_length)]
			
		else:
			sequence = ''
			for j in range(seq_length):
				sequence += get_random_AA()
				
		sequences.append(sequence)
	
	return sequences
