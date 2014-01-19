def create_res(in_pdb_path, out_res_path):
    
    infile =  open(in_pdb_path, 'r')
    outfile =  open(out_res_path, 'w')
    lines = infile.readlines()

def renumber(start_res_num, start_atom_num, in_pdb_path, out_pdb_path):
   
    infile =  open(in_pdb_path, 'r')
    outfile =  open(out_pdb_path, 'w')
    lines = infile.readlines()
    
    RN = int(start_res_num) #residue number
    AN = int(start_atom_num) #atom number
    
    FRN = -1; #first residue number
    
    for line in lines:
        line = line.strip()
        line = (line + ' '*16 + '\n') if (len(line) < 61) else (line + '\n') 
        
        if(line[0:4] == 'ATOM' or (line[0:6] == 'HETATM' and line[17:20] == 'MSE')):
            
            FRN = int(line[23:26]) if (FRN<0) else FRN  #if first run, set to non-zero first value
            CRN = int(line[23:26])                      #current residue number
            
            if CRN != FRN:  #checks if next residue has been reached
                RN += 1
                FRN = CRN
            
            #Replace old number with new residue count
            num_digits = len(str(RN))
            line = line[:21] + ' '*(5-num_digits) + str(RN) + line[26:]
            
            #Replace old number with new atom count
            num_digits = len(str(AN))
            line = line[:6] + ' '*(5-num_digits) + str(AN) + line[11:]
            AN += 1
            
            if(line[17:20] == 'MSE'):
                line = line[:18] + 'ET' + line[20:]
                line = 'ATOM  ' + line[6:]
        
            line = line[:56] + '1.00  0.00' + line[66:]
        
            outfile.write(line)
        
    outfile.write('TER\nEND')
    infile.close()
    outfile.close()
