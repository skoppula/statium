import bidict
import random
import math

def read_results(results_file, valueIsNum=True):
	results = filelines2deeplist(results_file, skipComments=True, useTuples=True, skipEmptyLines=True)
	results_dict = dict()
	for result in results:
		#remember seq stored in results file is already 'fixed' so we don't need to fix sequence it again
		if(valueIsNum):
			results_dict[result[0]] = float(result[-1])
		else:
			results_dict[result[0]] = result[-1]
	
	return results_dict


def list2file(t, infile):
	file = open(infile, 'w')
	for element in t:
		file.write(element + "\n")
	file.close()


def filelines2deeplist(infile, skipComments=False, useTuples=False, skipEmptyLines=False):
  
	t = []
	file = open(infile, 'r')
	lines = file.readlines()
	
	for line in lines:
		if(skipComments and line[0] == '#'):
			continue
		items = line.strip().split()
		if(skipEmptyLines and not items):
			continue
		if(useTuples):
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


#uses binary search to determine where element e ranks in a sorted list
def binary_placement(elements, e):
	left, right = 0, len(elements)-1
	mid = (left+right)/2
	
	if(elements[left] == e):
		return left
	elif(elements[right] == e):
		return right
	
	while(right-left > 1):
		if(e == elements[mid]):
			return mid
		elif(e > elements[mid]):
			left = mid
		else:
			right = mid
			
		mid = (left+right)/2
	
	return mid


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

	
def get_sidechain_atoms(a):
	
	sidechain_mappings = {'A':{'CB'}, 'C':{'CB', 'SG'}, 'D':{'CB', 'OD1', 'OD2'}, 'E':{'CB', 'OE1', 'OE2'}, 'F':{'CB', 'CE1', 'CE2', 'CZ'}, 'G':{'CA'},  
	'H':{'CB', 'CE1', 'NE2'}, 'I':{'CB', 'CG1', 'CG2', 'CD1'}, 'K':{'CB', 'NZ'}, 'L':{'CB', 'CD1', 'CD2'}, 'M':{'CB', 'CE'}, 
	'N':{'CB', 'OD1', 'ND2'}, 'P':{'CB', 'CG', 'CD'}, 'Q':{'CB', 'OE1', 'NE2'}, 'R':{'CB', 'NE', 'CZ'}, 'S':{'CB', 'OG'}, 
	'T':{'CB', 'OG1', 'CG2'}, 'V':{'CB', 'CG1', 'CG2'}, 'W':{'CB', 'CD1', 'CD2', 'NE1', 'CZ2', 'CZ3', 'CH2'}, 
	'Y':{'CB', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'}}
	
	return sidechain_mappings[a]


class Atom:
	def __init__(self, name, x, y, z):
		self.name = name
		self.x = x
		self.y = y
		self.z = z
		self.num = 0
		self.coordinates = (self.x, self.y, self.z)

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
		return hash((self.name, self.x, self.y, self.z))

	def __eq__(self, other):
		return (self.name, self.x, self.y, self.z) == (other.name, other.x, other.y, other.z)

	def __str__(self):
		return '(' + self.name + ', ' + str(self.num) + ',' + str(self.coordinates) + ')'


def magnitude(dx, dy, dz):
	return math.sqrt(dx**2+dy**2+dz**2)

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

		self.stubIntact = True if 'CA' in self.atom_names and 'CB' in self.atom_names else False

	#Returns dictionary of pairs of atoms mapped to each pair's distance
	def distancesTo(self, other):
		out = dict()
		for atom_one in self.atoms:
			for atom_two in other.atoms:
				out[(atom_one, atom_two)] = atom_one.distanceTo(atom_two)
				out[(atom_two, atom_one)] = atom_one.distanceTo(atom_two)
		return out

	def __hash__(self):
		return hash(self.string_name + self.chainID + str(self.atom_dict['CA'].coordinates))

	def __eq__(self, other):
		return (self.int_name, self.chainID, self.atom_dict['CA'].coordinates) == \
			(other.int_name, other.chainID, self.atom_dict['CA'].coordinates)

	def __str__(self):
		return '(' + self.string_name + ', ' + str(self.position) + ', ' + self.chainID + ')'

	def correct(self):
		if 'CA' not in self.atom_names:
			print 'Cannot correct residue %d with no alpha carbon.' % self.position
			return False
		else:
			ca = self.atom_dict['CA']
			n = self.atom_dict['N']
			c = self.atom_dict['C']
			cos_angle = math.cos(109.5)
			can_vector = (n.x-ca.x, n.y-ca.y, n.z-ca.z)
			(can.x, can.y, can.z) = can_vector
			cac_vector = (c.x-ca.x, c.y-ca.y, c.z-ca.z)
			(cac.x, cac.y, cac.z) = cac_vector

			n = magnitude(can_vector)
			m = magnitude(cac_vector)

			k1 = m*m - m*n*cac.x/can.x
			k2 = cac.y - cac.x*can.y/can.x
			k3 = cac.z - cac.x*can.z/can.x
			k13 = k1/k2
			k23 = k2/k3
			
			kappa = m**2*(1 - n**2/can.x**2)
			alpha = (1 + can.y**2/can.x**2)
			beta = (1 + can.z**2/can.x**2)
			gamma  = 2*can.y/can.x**2
			delta  = 2*can.y*m*n/can.x**2
			epsilon  = 2*can.z*m*n/can.x**2

			a = alpha + beta*k23**2 - gamma*k23
			b = -beta*k13*k23 + gamma*k13 + delta - epsilon*k23
			c = beta*k13**2 + epsilon*k13 - kappa

			y = solve(a,b,c)
			z = (k13-k23*y[0], k13-k23*y[1])
			x = (m*n/can.x - can.y/can.x*y[0] - can.z/can.x*z[0], m*n/can.x - can.y/can.x*y[1] - can.z/can.x*z[1])

			cb = (x[0], y[0], z[0])
			cacb_vector = (cb[0]-ca.x, cb[1]-ca.y, cb[2]-ca.z)
			cross = cross_product(ccab_vector, can_vector)
			dot = dot_product(cross, cac_vector)/(magnitude(cross)*m)
			
			atom = Atom('CB', x[0], y[0], z[0]) if dot > 0 else Atom('CB', x[1], y[1], z[1])
			self.atom_dict['CB'] = atom
			self.atoms.append(atom)
			print 'Corrected residue %d by adding CB' % self.position
			return True
		
		
#NEW VERSION:
#       Outputs list of Residues()
#OLD VERSION:
#       Creates a list with information for each residue:
#       e.g for each residue: [[1, ''], '16', [['CA', 'CB', 'OG1', 'CG2'], [[21.142, -19.229, 4.185], [21.957, -18.596, 5.322], [23.023, 17.818, 4.773], [22.547, -19.67, 6.206]]]]
def get_pdb_info(pdb_path):

	info = list()
	res = set()
 
	lines = filelines2list(pdb_path)
	
	curr_position = None
	curr_chainID = None
	curr_name = None
	curr_atoms = list()
	curr_possible_atoms = set()
	curr_found_atoms = set()

	atom_count = 1
	first_run = True

	for i, line in enumerate(lines):

		if line[0:4] == 'ATOM' or (line[0:6] == 'HETATM' and line[17:20] == 'MSE'):
			name = line[12:16].strip()
			if name ==  'N':
				if lines[i+1][12:16].strip() == 'CA':
					if not first_run: info.append(Residue(curr_name, curr_position, curr_chainID, curr_atoms))
					first_run = False
					curr_found_atoms.add(name)
					x = float(line[30:38].strip())
					y = float(line[38:46].strip())
					z = float(line[46:54].strip())			
					curr_atoms = [Atom(name, x, y, z, atom_count)]
					atom_count += 1
					continue
	
			if name == 'CA':
				try:
					curr_position = line[22:28].strip()
					curr_chainID = line[21:22].strip()
					curr_name = 'MET' if line[17:20] == 'MSE' else line[17:20]
					if (curr_chainID, curr_position) in res or not isAA(curr_name): continue
					res.add((curr_chainID, curr_position))
					
				except:
					continue

				x = float(line[30:38].strip())
				y = float(line[38:46].strip())
				z = float(line[46:54].strip())
				name = 'CA'
				curr_atoms .append(Atom(name, x, y, z, atom_count))
				atom_count += 1

				curr_possible_atoms = get_sidechain_atoms(AA2char(curr_name))
				curr_found_atoms = {'CA'}

			else:
				try:
					pos = line[22:28].strip()
					chain = line[21:22].strip()
				except: continue

				if chain != curr_chainID or pos != curr_position: continue
				
				if name not in curr_found_atoms or name in curr_found_atoms:# and name in curr_possible_atoms:
					curr_found_atoms.add(name)
					x = float(line[30:38].strip())
					y = float(line[38:46].strip())
					z = float(line[46:54].strip())
					curr_atoms.append(Atom(name, x, y, z, atom_count))
					atom_count += 1

	if curr_name:
		info.append(Residue(curr_name, curr_position, curr_chainID, curr_atoms))
	
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
