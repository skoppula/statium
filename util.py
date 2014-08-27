import bidict
import re
import random
import math
import warnings

def generate_random_seq(seq_length, library_path='data/all_protein_sequences.txt'): 
	
	size = 0
	try:
		with open (library_path, "r") as myfile:
			data = myfile.read()
			data = re.sub('>.+?\n', '', data)
			pattern = re.compile(r'\s+')
			data = re.sub(pattern, '', data)
			size = len(data)
	except IOError:
		print 'Could not find file ' + library_path
		
	start = random.randint(0, size - seq_length)
	sequence = data[start:(start+seq_length)]
	
	return sequence


def calc_seq_zscore(mean, std, energy):
	zscore = (energy - mean)/std
	return zscore


def list2file(t, infile):
	with open(infile, 'w') as f:
		for element in t:
			f.write(str(element) + "\n")


def filelines2deeplist(infile, skipComments=False, useTuples=False, skipEmptyLines=False):
	t = []
	file = open(infile, 'r')
	with open(infile, 'r') as f:
		lines = f.read().split('\n')
	
	for line in lines:
		items = line.strip().split()
		if (skipEmptyLines and not items) or (skipComments and line != '' and line[0] == '#'):
			continue
		if useTuples:
			items = tuple(items)
		t.append(items)

	return t


#Returns lines of file in a list
def filelines2list(infile):
	t = []
	file = open(infile, 'r')
	lines = file.readlines()
	
	for line in lines:
		t.append(line.strip())
		
	return t


def mean(s):
	return sum(s) * 1.0/len(s)


def std(s):
	avg = mean(s)
	variance = map(lambda x: (x - avg)**2, s)
	return math.sqrt(mean(variance))

	
def nCr(n,r):
	f = math.factorial
	return f(n) / f(r) / f(n-r)


def get_random_AA():
	library_AA_distro = {'A':1391008, 'C':379677, 'D':825255 , 'E':911171 , 'F':1415079 , 'G':947640 , 'H':522382 , 'I':1838459 , 'K':822707 , 'L':2756022 , 'M':491788 , 'N':708103 , 'P':808079 , 'Q':614885 , 'R':983182 , 'S':954331 , 'T':1090982 , 'V':2059360 , 'W':539073 , 'Y':1190549 , 'X':0}
	maxFreq = max(library_AA_distro.values())
	normalizer = 10 ** (len(str(maxFreq))-1)
	total = sum([x/normalizer for x in library_AA_distro.values()])

	num = random.uniform(0, total)
	s = 0
	for AA in library_AA_distro:
		s += library_AA_distro[AA]/normalizer
		if(num < s):
			return AA

	return 'X'


def AA2char(a):
	AA = bidict.bidict({'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F', 'GLY':'G', 'ILE':'I', 'LYS':'K', 'LEU':'L', 'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R', 'SER':'S', 'THR':'T', 'VAL':'V', 'TRP':'W', 'TYR':'Y', 'MSE':'X', 'HIS':'H'})
	return AA[a]
	
def char2AA(a):
	AA = bidict.bidict({'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F', 'GLY':'G', 'ILE':'I', 'LYS':'K', 'LEU':'L', 'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R', 'SER':'S', 'THR':'T', 'VAL':'V', 'TRP':'W', 'TYR':'Y', 'MSE':'X', 'HIS':'H'})
	return AA[:a]

def isAA(a):
	AA = {'ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR', 'MSE', 'HIS'}
	return (a in AA)

def AAchar2int(a):
	AA = bidict.bidict({'0':'A', '1':'C', '2':'D', '3':'E', '4':'F', '5':'G', '7':'I', '8':'K', '9':'L', '10':'M', '11':'N', '12':'P', '13':'Q', '14':'R', '15':'S', '16':'T', '17':'V', '18':'W', '19':'Y', '20':'X', '6':'H'})
	return int(AA[:a])

def AAint2char(a):
	AA = bidict.bidict({'0':'A', '1':'C', '2':'D', '3':'E', '4':'F', '5':'G', '7':'I', '8':'K', '9':'L', '10':'M', '11':'N', '12':'P', '13':'Q', '14':'R', '15':'S', '16':'T', '17':'V', '18':'W', '19':'Y', '20':'X', '6':'H'})
	return AA[str(a)]

def get_sidechain_atoms(res):
	
	sidechain_mappings = {'A':{'CA', 'CB'}, 'C':{'CA', 'CB', 'SG'},
							'D':{'CA', 'CB', 'CG', 'OD1', 'OD2'},
							'E':{'CA', 'CB', 'CG', 'CD', 'OE1', 'OE2'},
							'F':{'CA', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'},
							'G':{'CA', 'CB'},  
							'H':{'CA', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'},
							'I':{'CA', 'CB', 'CG1', 'CG2', 'CD1'},
							'K':{'CA', 'CB', 'CD', 'CG', 'CE', 'NZ'},
							'L':{'CA', 'CB', 'CG', 'CD1', 'CD2'},
							'M':{'CA', 'CB', 'CG', 'SD', 'CE'}, 
							'N':{'CA', 'CB', 'CG', 'OD1', 'ND2'},
							'P':{'CA', 'CB', 'CG', 'CD'},
							'Q':{'CA', 'CB', 'CG', 'CD', 'OE1', 'NE2'},
							'R':{'CA', 'CB', 'CD', 'CG', 'NE', 'CZ', 'NH1', 'NH2'},
							'S':{'CA', 'CB', 'OG'}, 
							'T':{'CA', 'CB', 'OG1', 'CG2'},
							'V':{'CA', 'CB', 'CG1', 'CG2'},
							'W':{'CA', 'CB', 'CG', 'CD1', 'CD2', 'CE2', 'NE1', 'CZ2', 'CE3', 'CZ3', 'CH2'}, 
							'Y':{'CA', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'}}
	
	return sidechain_mappings[res]


class Atom:
	def __init__(self, name, x, y, z, num):
		self.name = name
		self.x = x
		self.y = y
		self.z = z
		self.num = num
		self.coordinates = (self.x, self.y, self.z)


	def distanceTo(self, atom):
		return math.sqrt((self.x-atom.x)**2+(self.y-atom.y)**2+(self.z-atom.z)**2)
	
	def __hash__(self):
		return hash((self.x, self.y, self.z))

	def __eq__(self, other):
		return (self.name, self.x, self.y, self.z) == (other.name, other.x, other.y, other.z)

	def __str__(self):
		return '(' + self.name + ', ' + str(self.num) + ',' + str(self.coordinates) + ')'

#	def __repr__(self):
#		return self.name + '[' + str(self.num) + ']'

	def sameName(self, other):
		return True if self.name == other.name else False


def magnitude(x, y, z):
	return math.sqrt(x**2+y**2+z**2)

def quadsolve(a, b, c):
	disc = math.sqrt(b*b-4*a*c)
	return ((-b + disc)/(2*a), (-b+disc)/(2*a))

def cross(a, b):
	return (a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0])

def dot(a, b):
	return sum(p*q for p,q in zip(a, b))

class Residue:
	#string_name is the three letter amino acid identifier
	#position is the position along the chain in the renumbered PDB
	#chainID is the chain identifier (duh)
	#atoms is a list of Atom() objects in the order listed in the corresponding PDB
	def __init__(self, string_name, position, chainID, atoms):
		self.string_name = string_name
		self.char_name = AA2char(string_name)
		self.int_name = AAchar2int(self.char_name)
		
		self.position = position
		self.chainID = chainID
		self.atoms = atoms

		self.atom_names = list()
		self.atom_dict = dict()
		for a in atoms:
			self.atom_names.append(a.name)
			self.atom_dict[a.name] = a	
                self.atom_names.sort()

		self.stubIntact = True if 'CA' in self.atom_names and 'CB' in self.atom_names else False

	#Returns dictionary of pairs of atoms mapped to each pair's distance
	#DO NOT USE THIS. NOT THIS.
	def filteredDistancesTo(self, other, cutoff, output_dict=True):
                isIP = False

                if output_dict:
                    out = dict()
                    for atom_one in self.atoms:
                            for atom_two in other.atoms:
                                    dist = atom_one.distanceTo(atom_two)
                                    if dist < cutoff: isIP = True
                                    out[(atom_one, atom_two)] = dist
                                    out[(atom_two, atom_one)] = dist

                    return out if isIP else None

                else:
                    dists = list()
                    pairs = list()
                    for atom_one in self.atoms:
                            for atom_two in other.atoms:
                                    dist = atom_one.distanceTo(atom_two)
                                    if dist < cutoff: isIP = True
                                    dists.append(dist)
                                    pairs.append((atom_one, atom_two))

                    return [pairs, dists] if isIP else None

	def fastFilteredDistancesTo(self, other, cutoff, output_dict=True):
                isIP = False

                if output_dict:
		    out = dict()

                    for atom_one in self.atoms:
                            for atom_two in other.atoms:
                                    dist = atom_one.distanceTo(atom_two)
                                    if dist < cutoff: isIP = True
	         		    out[(atom_one.name, atom_two.name)] = dist
		    return out if isIP else None

                else:
		    dists = list()
		    pairs = list()

                    for atom_one in self.atoms:
                            for atom_two in other.atoms:
                                    dist = atom_one.distanceTo(atom_two)
                                    if dist < cutoff: isIP = True
				    dists.append(dist)
				    pairs.append((atom_one.name, atom_two.name))

		    return [pairs, dists] if isIP else None

	def distancesTo(self, other, output_dict=True):
                if output_dict:
                    out = dict()
                    for atom_one in self.atoms:
                            for atom_two in other.atoms:
                                    dist = atom_one.distanceTo(atom_two)
                                    out[(atom_one, atom_two)] = dist
                                    out[(atom_two, atom_one)] = dist
                    return out

                else:
                    dists = list()
                    pairs = list()
                    for atom_one in self.atoms:
                            for atom_two in other.atoms:
                                    dist = atom_one.distanceTo(atom_two)
                                    dists.append(dist)
                                    pairs.append((atom_one, atom_two))

                    return [pairs, dists]

	def fastDistancesTo(self, other, output_dict=True):
                if output_dict:
                    out = dict()
                    for atom_one in self.atoms:
                            for atom_two in other.atoms:
                                    dist = atom_one.distanceTo(atom_two)
                                    out[(atom_one.name, atom_two.name)] = dist
                    return out

                else:
                    dists = list()
                    pairs = list()
                    for atom_one in self.atoms:
                            for atom_two in other.atoms:
                                    dist = atom_one.distanceTo(atom_two)
                                    dists.append(dist)
                                    pairs.append((atom_one.name, atom_two.name))

                    return [pairs, dists]


	def __hash__(self):
		return hash(self.string_name + self.chainID + str(self.atom_dict['CA'].coordinates))

	def __eq__(self, other):
		return (self.int_name, self.chainID, self.atom_dict['CA'].coordinates) == \
			(other.int_name, other.chainID, self.atom_dict['CA'].coordinates)

	def __str__(self):
		return '(' + self.string_name + ', ' + str(self.position) + ', ' + self.chainID + ')'

#	def __repr__(self):
#		return self.string_name + '[' + str(self.position) + ']'

        def isStubIntact(self):
		self.stubIntact = True if 'CA' in self.atom_names and 'CB' in self.atom_names else False
                return self.stubIntact

	def correct(self):
		from scipy.optimize import fsolve
		if 'CA' not in self.atom_names or 'C' not in self.atom_names or 'N' not in self.atom_names:
			print 'Cannot correct residue %s with no alpha carbon, nitrogen, or secondary carbon.' % self.position
			return False
		else:
			#uses vector analysis to place CB in correct tetrahedral place
			#	according to coordinate geometry/SP3 CA hybridization
			(p1x, p1y, p1z) = self.atom_dict['N'].coordinates
			(p2x, p2y, p2z) = self.atom_dict['CA'].coordinates
			(p3x, p3y, p3z) = self.atom_dict['C'].coordinates

			ca_n = (p1x-p2x,p1y-p2y,p1z-p2z)
			ca_c = (p3x-p2x,p3y-p2y,p3z-p2z)
			d_can = magnitude(*ca_n)
			d_cac = magnitude(*ca_c)
	
			def equations(x):
				c = math.cos(math.radians(109.5))
				out = [ca_n[0]*x[0]+ca_n[1]*x[1]+ca_n[2]*x[2]-d_can*d_cac*c]
				out.append(ca_c[0]*x[0]+ca_c[1]*x[1]+ca_c[2]*x[2]-d_cac*d_cac*c)
				out.append(x[0]**2+x[1]**2+x[2]**2-d_can*d_cac)
				return out
	
			init1 = [p2x+1,p2y,p2z]
			init2 = [p2x-1,p2y,p2z]
			with warnings.catch_warnings():
				warnings.simplefilter("ignore")
				ca_x_1 = fsolve(equations, init1)
				ca_x_2 = fsolve(equations, init2)
	
			if all(ca_x_1 == ca_x_2):
				init1 = [p2x,p2y+1,p2z]
				init2 = [p2x,p2y-1,p2z]
				with warnings.catch_warnings():
					warnings.simplefilter("ignore")
					ca_x_1 = fsolve(equations, init1)
					ca_x_2 = fsolve(equations, init2)
	
			product = cross(ca_x_1, ca_n)
			product = dot(product, ca_c)
	
			x_vec = ca_x_1 if product > 1 else ca_x_2
		
			atom = Atom('CB', p2x+x_vec[0], p2y+x_vec[1], p2z+x_vec[2], 0)
			self.atom_dict['CB'] = atom
			self.atoms.append(atom)
                        self.atom_names.append('CB')
			self.stubIntact = True
			print '\tCorrected residue %s by adding CB' % self.position
			return True
        
        def strip_backbone(self):
            def is_backbone(atom):
                if atom.name == 'N' or atom.name == 'C':
                    del self.atom_dict[atom.name]
                    self.atom_names.remove(atom.name)
                    return False
                else:
                    return True

            self.atoms = filter(is_backbone, self.atoms)

		
		
def get_pdb_info(pdb_path, filter_sidechains = False, include_backbone = True):

	info = list()
	res = set()
 
	lines = filelines2list(pdb_path)
	
	curr_position = None
	curr_chainID = None
	curr_res_name = None
	curr_atoms = list()
	curr_possible_atoms = set()
	curr_found_atoms = set()

	atom_count = 1
	first_run = True

	for i, line in enumerate(lines):

		if line[0:4] == 'ATOM' or (line[0:6] == 'HETATM' and line[17:20] == 'MSE'):
			res_name = 'MET' if line[17:20] == 'MSE' else line[17:20]
			if not isAA(res_name): continue

			try:
				position = line[22:28].strip()
				chain = line[21:22].strip()
				atom_name = line[12:16].strip()
				x = float(line[30:38].strip())
				y = float(line[38:46].strip())
				z = float(line[46:54].strip())			
			except IndexError:
				continue

			if curr_position == position and chain == curr_chainID and curr_res_name == res_name:
				#print [atom.name for atom in curr_atoms]
				#print curr_found_atoms
				if atom_name not in curr_found_atoms:
					if not filter_sidechains or (atom_name in curr_possible_atoms):
						if include_backbone or (atom_name != 'N' and atom_name != 'C'):
							curr_atoms.append(Atom(atom_name, x, y, z, atom_count))
							curr_found_atoms.add(atom_name)
							atom_count += 1
				
			else:
				if not first_run and 'CA' in curr_found_atoms:
					info.append(Residue(curr_res_name, curr_position, curr_chainID, curr_atoms))
					res.add((curr_chainID, curr_position))
				first_run = False

				curr_position = position
				curr_chainID = chain
				curr_res_name = res_name
				curr_atoms = [Atom(atom_name, x, y, z, atom_count)]
				curr_found_atoms = {atom_name}
				atom_count += 1
				
				if isAA(curr_res_name):
					curr_possible_atoms = get_sidechain_atoms(AA2char(curr_res_name))

	info.append(Residue(curr_res_name, curr_position, curr_chainID, curr_atoms))
	
	return info

def print_pdb(residues, path):
	outfile = open(path, 'w')
	for residue in residues:
		rname = (3-len(residue.string_name))*' ' + residue.string_name
		rchainID = residue.chainID
		rnum = (5-len(str(residue.position)))*' ' + residue.position
		
		for atom in residue.atoms:
			anum = (5-len(str(atom.num)))*' ' + str(atom.num)
			aname =  atom.name + (4-len(atom.name))*' '
			ax = (8-len(str(atom.x)))*' ' + str(atom.x)
			ay = (8-len(str(atom.y)))*' ' + str(atom.y)
			az = (8-len(str(atom.z)))*' ' + str(atom.z)
			line = 'ATOM  ' + anum + ' ' + aname + ' ' + rname + ' ' + rchainID + rnum + '   ' + ax + ay + az + '\n'
			outfile.write(line)

	outfile.write('TER\nEND')
	outfile.close()
