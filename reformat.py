from util import AA2char
from util import get_random_AA
from util import get_pdb_info
from util import print_pdb
import re
import random
def renumber(start_res_num, start_atom_num, chains, in_pdb, out_pdb):
   
	RN = start_res_num 
	AN = start_atom_num
	
	residues = get_pdb_info(in_pdb)
	new_residues = list()

	#Putting non-mutating chains first
	for residue in residues:
		if residue.chainID in chains: continue
		residue.position = RN
		for atom in residue.atoms:
			atom.num = AN
			AN += 1
		RN += 1
		new_residues.append(residue)

	#Put mutating chains second
	for residue in residues:
		if residue.chainID not in chains: continue
		residue.position = RN
		for atom in residue.atoms:
			atom.num = AN
			AN += 1
		RN += 1
		new_residues.append(residue)

	print_pdb(new_residues, out_pdb)

def create_res(pdb_orig_path, pdb_renum_path, out_res_path, positions):
	
	orig =  get_pdb_info(pdb_orig_path)
	renum =  get_pdb_info(pdb_renum_path)

	#Adding all ligand residues we're interesting in mutating to this set()
	#	and converting them to residues in new_res
	res = set()
	new_res = set()

	for position in positions:
		if isinstance(position, tuple):
			(chain, start, end) = position
			i = 0
			while orig[i].chainID != chain and orig[i].position != start:
				i += 1
			while orig[i].chainID != chain and orig[i].position != end:
				res.add(orig[i])
				i += 1
			res.add(orig[i])

		elif len(position) == 1:
			inchain = False
			for residue in orig:
				if residue.chainID == position:
					inchain = True
					res.add(residue)	
				else:
					if inchain: break
		else:
			for residue in orig:
				if residue.chainID == position[0] and residue.position == position[1:]:
					res.add(residue)
					break

	for residue in res:
		i = renum.index(residue)
		new_res.add(i)

	out =  open(out_res_path, 'w')
	for residue in new_res:
		out.write(residue.position + '\n')
	
	out.close()


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
