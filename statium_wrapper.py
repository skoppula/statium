import sys
from docopt import docopt

def renumber(start_res_num, start_atom_num, in_pdb_path, out_pdb_path):
    
    infile =  open(out_pdb_path, 'r')
    outfile =  open(in_pdb_path, 'w')
    lines = infile.readlines()
    
    RN = int(start_res_num) #residue number
    AN = int(start_atom_num) #atom number
    
    FRN = -1; #number of first residue
    
    for line in lines:
        line = line.strip()
        line = (line + ' '*16 + '\n') if (len(line) < 61) else (line + '\n') 
        
        if(line[0:4] == 'ATOM' or (line[0:6] == 'HETATM' and line[17:20] == 'MSE')):
            
            FRN = int(line[23:26]) if (FRN<0) else FRN  #if first run, set to non-zero first value
            CRN = int(line[23:26])                      #current residue number
            
            if CRN != FRN:  #checks if next residue has been reached
                RN += 1
                CRN = FRN
            
            num_digits = len(str(RN))
            line = line[:20] + ' '*5 + line[25:]
            line = line[:(25-num_digits)] + RN + line[25:]

            num_digits = len(str(AN))
            line = line[:5] + ' '*5 + line[10:]
            line = line[:(10-num_digits)] + (AN) + line[10:]
            AN += 1
            
            if(line[17:20] == 'MSE'):
                line = line[:18] + 'ET' + line[20:]
                line = 'ATOM  ' + line[6:]
        
            line = line[:56] + '1.00  0.00' + line[65:]
        
            outfile.write(line)
        
    outfile.write('TER\nEND')
    infile.close()
    outfile.close()

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