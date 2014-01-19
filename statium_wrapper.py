import sys
from docopt import docopt

def renumber(start_res_num, start_atom_num, in_pdb_path, out_pdb_path):
    
    infile =  open(out_pdb_path, 'r')
    outfile =  open(in_pdb_path, 'w')
    lines = infile.readlines()
    
    res_num = int(start_res_num)
    atom_num = start_atom_num
    atom_count = 0
    
    l_1 = 'false'
    
    for line in lines:
        lines = line.strip()
        if len(line) < 61:
            line = line + '                \n'
        else:
            line = line + '\n'
        line_list = list(line)
    
    if (line_list[0] == 'A' and line_list[1] == 'T' and line_list[2] == 'O' and line_list[3] == 'M') or (line_list[0] == 'H' and line_list[1] == 'E' and line_list[2] == 'T' and line_list[3] == 'A' and line_list[4] == 'T' and line_list[5] == 'M' and line_list[17] == 'M' and line_list[18] == 'S' and line_list[19] == 'E'):
        if l_1 == 'false':
            hold_res_num = int(pdb_lines[i][23:26])
        l_1 = 'true'
        current_res_num = int(pdb_lines[i][23:26])
        if current_res_num == hold_res_num:
            res_size = len(str(int(res_num)))
            for j in range(5):
                line_list[25 - j] = ' '
            for j in range(res_size):
                line_list[25 - j] = str(res_num)[res_size - j - 1]
        else: 
            res_num = res_num + 1
            hold_res_num = current_res_num
            res_size = len(str(res_num))
            for j in range(5):
                line_list[25 - j] = ' '
            for j in range(res_size):
                line_list[25 - j] = str(res_num)[res_size - j - 1]
    
        atom_size = len(str(int(atom_num) + atom_count))
        for j in range(5):
            line_list[10 - j] = ' '
        for j in range(atom_size):
            line_list[10 - j] = str(int(atom_num) + atom_count)[atom_size - j - 1]
        
        if line_list[17] == 'M' and line_list[18] == 'S' and line_list[19] == 'E':
        line_list[18] = 'E'
        line_list[19] = 'T'
        
        line_list[0] = 'A'
        line_list[1] = 'T'
        line_list[2] = 'O'
        line_list[3] = 'M'
        line_list[4] = ' '
        line_list[5] = ' '

       # line_list[21] = 'A'
        
        line_list[56] = '1'
        line_list[57] = '.'
        line_list[58] = '0'
        line_list[59] = '0'
        line_list[60] = ' '
        line_list[61] = ' '        
        line_list[62] = '0'
        line_list[63] = '.'
        line_list[64] = '0'
        line_list[65] = '0'
        joined_list = string.join(line_list, '')
        renumber_file.write(joined_list)
        atom_count = atom_count + 1
    renumber_file.write('TER\n')
    renumber_file.write('END')

    pdb_file.close()
    renumber_file.close()

def renumber(in_pdb, out_pdb):
    #         Renumber(start_res_num, start_atom_num, pdb_path, renumber_path)

def main(argv):
    
    helpdoc =   """
                usage: statium_wrapper.py renumber (<IN_PDB>) [<OUT_PDB>] [-v | --verbose]
                       statium_wrapper.py create_res (<IN_PDB>) [<OUT_RES>] [-v | --verbose]
                       statium_wrapper.py run_statium (<IN_CFG> <IN_RES> <IN_PDB>) [-v | --verbose]
                       statium_wrapper.py calc_seq_energy (<DIR> <SEQUENCE>) [-v | --verbose]
                       statium_wrapper.py [-h | --help]
                """
    
    options = docopt(helpdoc, argv, help = True, version = "1.0.0", options_first=False)
    print(str(options))
    
    if(options['renumber']):
        if(options['OUT_PDB'] == False):
            renumber(1, 1, options['IN_PDB'], options['IN_PDB'][:-4]+'_renumbered.pdb')
        else:
            renumber(1, 1, options['IN_PDB'], options['OUT_PDB'])

    elif(options['create_res']):
        
    elif(options['run_statium']):
        
    elif(options['calc_seq_energy']):
        
        

if __name__ == "__main__":
    main(sys.argv[1:])