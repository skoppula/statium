from util import AA2char
from util import get_random_AA
from util import get_pdb_info
from util import print_pdb
from util import filelines2list
import re
import random
import sys

def renumber(start_res_num, start_atom_num, chains, in_pdb, out_pdb):
   
	RN = start_res_num 
	AN = start_atom_num
	
	residues = get_pdb_info(in_pdb)
	new_residues = list()

	#Putting non-mutating chains first
	for residue in residues:
		if residue.chainID in chains: continue
		residue.position = str(RN)
		for atom in residue.atoms:
			atom.num = AN
			AN += 1
		RN += 1
		new_residues.append(residue)

	#Put mutating chains second
	for residue in residues:
		if residue.chainID not in chains: continue
		residue.position = str(RN)
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
			while not (orig[i].chainID == chain and orig[i].position == start):
				i += 1
			while not (orig[i].chainID == chain and orig[i].position == end):
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
		i = 0
		for r in renum:
			if r.atom_dict['CA'].coordinates == residue.atom_dict['CA'].coordinates and r.int_name == residue.int_name and r.chainID == residue.chainID:
				new_res.add(r)
				break

	out =  open(out_res_path, 'w')
	positions = [int(residue.position) for residue in new_res]	
	positions.sort()	
	for position in positions:
		out.write(str(position) + '\n')
	
	out.close()


def get_orig_seq(res_path, orig_pdb_path, renum_pdb_path):
	print 'Extracting residue position from ' + res_path + '...'
	res_lines = filelines2list(res_path)
	residues = [int(line.strip()) for line in res_lines]

	orig =  get_pdb_info(orig_pdb_path)
	renum =  get_pdb_info(renum_pdb_path)
	to_print = set()

	for residue in residues:
		renum_aa = renum[residue-1]
		for orig_aa in orig:
			if renum_aa == orig_aa:
				to_print.add(orig_aa)

	to_print_tuples = list()
	for res in to_print:
		try:
			to_print_tuples.append((res.chainID,int(res.position),res.string_name, res.char_name))
		except:
			to_print_tuples.append((res.chainID,res.position,res.string_name, res.char_name))
	to_print_tuples.sort()

	for t in to_print_tuples:
		print t[2] + ' ' + t[0] + ' ' + str(t.[1])

	for t in to_print_tuples:
		sys.stdout.write(t[3])

	print ''
				


