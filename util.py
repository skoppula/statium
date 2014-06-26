import bidict
import random
import math
import warnings

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

		self.stubIntact = True if 'CA' in self.atom_names and 'CB' in self.atom_names else False

	#Returns dictionary of pairs of atoms mapped to each pair's distance

	def filteredDistancesTo(self, other, cutoff):
		out = dict()
		isIP = False
		for atom_one in self.atoms:
			for atom_two in other.atoms:
				dist = atom_one.distanceTo(atom_two)
				if dist < cutoff: isIP = True
				out[(atom_one, atom_two)] = dist
				out[(atom_two, atom_one)] = dist

		return out if isIP else None

	def distancesTo(self, other):
		out = dict()
		for atom_one in self.atoms:
			for atom_two in other.atoms:
				dist = atom_one.distanceTo(atom_two)
				out[(atom_one, atom_two)] = dist
				out[(atom_two, atom_one)] = dist
		return out

	def __hash__(self):
		return hash(self.string_name + self.chainID + str(self.atom_dict['CA'].coordinates))

	def __eq__(self, other):
		return (self.int_name, self.chainID, self.atom_dict['CA'].coordinates) == \
			(other.int_name, other.chainID, self.atom_dict['CA'].coordinates)

	def __str__(self):
		return '(' + self.string_name + ', ' + str(self.position) + ', ' + self.chainID + ')'

#	def __repr__(self):
#		return self.string_name + '[' + str(self.position) + ']'

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
			self.stubIntact = True
			print '\tCorrected residue %s by adding CB' % self.position
			return True
		
		
#NEW VERSION:
#	   Outputs list of Residues()
#OLD VERSION:
#	   Creates a list with information for each residue:
#	   e.g for each residue: [[1, ''], '16', [['CA', 'CB', 'OG1', 'CG2'], [[21.142, -19.229, 4.185], [21.957, -18.596, 5.322], [23.023, 17.818, 4.773], [22.547, -19.67, 6.206]]]]
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
				try:
					ca_check1 = lines[i+1][12:16].strip() == 'CA' and line[22:28].strip()==lines[i+1][22:28].strip()
					ca_check2 = lines[i+3][12:16].strip() == 'CA' and line[22:28].strip()==lines[i+3][22:28].strip()
					if ca_check1 or ca_check2:
						#only add to residue list if there are actual values in curr vars
						#	and circumvent weird case where first residue has two N's
						#	for two different conformations
						if not first_run and isAA(curr_name) and (curr_chainID, curr_position) not in res:
							info.append(Residue(curr_name, curr_position, curr_chainID, curr_atoms))
							res.add((curr_chainID, curr_position))

						first_run = False
						curr_found_atoms.add(name)
						x = float(line[30:38].strip())
						y = float(line[38:46].strip())
						z = float(line[46:54].strip())			
						curr_position = line[22:28].strip()
						curr_chainID = line[21:22].strip()
						curr_atoms = [Atom(name, x, y, z, atom_count)]
						atom_count += 1

						continue

				except IndexError:
					continue	
		
			if name == 'CA':
				try:
					curr_position = line[22:28].strip()
					curr_chainID = line[21:22].strip()
					curr_name = 'MET' if line[17:20] == 'MSE' else line[17:20]
					
				except:
					continue

				x = float(line[30:38].strip())
				y = float(line[38:46].strip())
				z = float(line[46:54].strip())
				name = 'CA'
				curr_atoms.append(Atom(name, x, y, z, atom_count))
				atom_count += 1

				if isAA(curr_name):
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

	if isAA(curr_name):
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
