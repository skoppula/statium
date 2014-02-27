import bidict
import random
import math

def list2file(t, infile):
    file = open(infile, 'w')
    for element in t:
        file.write(element + "\n")
    file.close()

def filelines2deeplist(infile):
  
    t = []
    file = open(infile, 'r')
    lines = file.readlines()
    
    for line in lines:
        items = line.strip().split()
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

def AA_cutoff_dist(AA):
    
    if(AA == 'A' or AA == 'G'):
        return 0.2
    else:
        return 0.4

def get_random_AA():
    library_AA_distro = dict({'A':1391008, 'C':379677, 'D':825255 , 'E':911171 , 'F':1415079 , 'G':947640 , 'H':522382 , 'I':1838459 , 'K':822707 , 'L':2756022 , 'M':491788 , 'N':708103 , 'P':808079 , 'Q':614885 , 'R':983182 , 'S':954331 , 'T':1090982 , 'V':2059360 , 'W':539073 , 'Y':1190549 , 'X':0})
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
    AA = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR', 'MSE', 'HIS']
    return (a in AA)

def AAchar2int(a):
    AA = bidict.bidict({'0':'A', '1':'C', '2':'D', '3':'E', '4':'F', '5':'G', '7':'I', '8':'K', '9':'L', '10':'M', '11':'N', '12':'P', '13':'Q', '14':'R', '15':'S', '16':'T', '17':'V', '18':'W', '19':'Y', '20':'X', '6':'H'})
    return int(AA[:a])

def AAint2char(a):
    AA = bidict.bidict({'0':'A', '1':'C', '2':'D', '3':'E', '4':'F', '5':'G', '7':'I', '8':'K', '9':'L', '10':'M', '11':'N', '12':'P', '13':'Q', '14':'R', '15':'S', '16':'T', '17':'V', '18':'W', '19':'Y', '20':'X', '6':'H'})
    return AA[str(a)]

    
def get_sidechain_atoms(a):
    
    sidechain_mappings = {'A':['CB'], 'C':['CB', 'SG'], 'D':['CB', 'OD1', 'OD2'], 'E':['CB', 'OE1', 'OE2'], 'F':['CB', 'CE1', 'CE2', 'CZ'], 'G':['CA'],  
    'H':['CB', 'CE1', 'NE2'], 'I':['CB', 'CG1', 'CG2', 'CD1'], 'K':['CB', 'NZ'], 'L':['CB', 'CD1', 'CD2'], 'M':['CB', 'CE'], 
    'N':['CB', 'OD1', 'ND2'], 'P':['CB', 'CG', 'CD'], 'Q':['CB', 'OE1', 'NE2'], 'R':['CB', 'NE', 'CZ'], 'S':['CB', 'OG'], 
    'T':['CB', 'OG1', 'CG2'], 'V':['CB', 'CG1', 'CG2'], 'W':['CB', 'CD1', 'CD2', 'NE1', 'CZ2', 'CZ3', 'CH2'], 
    'Y':['CB', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH']}
    
    return(sidechain_mappings[a])