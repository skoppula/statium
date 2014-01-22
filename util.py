import sys

def list2file(l, infile):
    file = open(infile, 'w')
    for element in l:
        file.write(element + "\n")
    file.close()

'''    
def lines2list(infile):
  
    v1 = []
    lines = readlines(infile)
    for i in range(len(lines)):
        items = lines[i].strip().split()
        sub = []
        for j in range(len(items)):
        try: sub.append(float(items[j]))
        except: sub.append(items[j])
    v1.append(sub)

    return v1
'''
    
def readlines(file_path):
    
    file = open(file_path, 'r')
    lines = file.readlines()
    file.close()
    return lines