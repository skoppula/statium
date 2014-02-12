import bidict

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
        if(line[0] == '#'):
            continue;
        t.append(line.strip())
        
    return t

def AA_cutoff_dist(AA):
    
    if(AA == 'A' or AA == 'G'):
        return 0.2
    else:
        return 0.4

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