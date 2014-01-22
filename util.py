def write_list_to_file(l, infile):
    file = open(infile, 'w')
    for element in l:
        file.write(element + "\n")
    file.close()