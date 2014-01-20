import sys
from docopt import docopt
from statium_reformat import renumber
from statium_reformat import create_res
from statium_reformat import create_cfg 

def main(argv):
    
    helpdoc =   """
                usage: statium_wrapper.py renumber (IN_PDB) [OUT_PDB] [-v | --verbose]
                       statium_wrapper.py create_res (IN_PDB_ORIG IN_PDB_RENUMBERED) [OUT_RES] [-v | --verbose]
                       statium_wrapper.py create_cfg (IN_RES) [OUT_CFG] [-v | --verbose]
                       statium_wrapper.py run_statium (IN_CFG IN_RES IN_PDB) [-v | --verbose]
                       statium_wrapper.py calc_seq_energy (DIR SEQ) [-v | --verbose]
                       statium_wrapper.py [-h | --help]
                """
    
    options = docopt(helpdoc, argv, help = True, version = "2.0.0", options_first=False)
    
    if(options['renumber']):
        if(options['-v'] or options['--verbose']): print("Reformatting file: " + options['IN_PDB'])
        
        if(options['OUT_PDB'] == None):
            renumber(1, 1, options['IN_PDB'], options['IN_PDB'][:-4]+'_renumbered.pdb')
            if(options['-v'] or options['--verbose']): print("Done. Formatted file: " + options['IN_PDB'][:-4]+'_renumbered.pdb')
            
        else:
            renumber(1, 1, options['IN_PDB'], options['OUT_PDB'])
            if(options['-v'] or options['--verbose']): print("Done. Formatted file: " + options['OUT_PDB'])
    
        
    elif(options['create_res']):
        if(options['-v'] or options['--verbose']): print("Creating .res file using: " + options['IN_PDB_ORIG'] + " and " + options['IN_PDB_RENUMBERED'])
        
        if(options['OUT_RES'] == None):
            create_res(options['IN_PDB_ORIG'], options['IN_PDB_RENUMBERED'], options['IN_PDB_RENUMBERED'][:-4]+'_renumbered.pdb')
            if(options['-v'] or options['--verbose']): print("Done. .res file: " + options['IN_PDB_RENUMBERED'][:-4]+'.res')
            
        else:
            create_res(options['IN_PDB_ORIG'], options['IN_PDB_RENUMBERED'], options['OUT_RES'])
            if(options['-v'] or options['--verbose']): print("Done. .res file: " + options['OUT_RES'])
    
    
    elif(options['create_cfg']):
        if(options['-v'] or options['--verbose']): print("Creating .cfg file using: " + options['IN_RES'])
        
        if(options['OUT_CFG'] == None):   
            create_cfg(options['IN_RES'], options['IN_RES'][:-4]+'.cfg')
            if(options['-v'] or options['--verbose']): print("Done. .cfg file: " + options['IN_RES'][:-4]+'.cfg')
            
        else:
            create_cfg(options['IN_RES'], options['OUT_CFG'])
            if(options['-v'] or options['--verbose']): print("Done. .cfg file: " + options['OUT_CFG'])
        
        
    elif(options['run_statium']):
        if(options['-v'] or options['--verbose']): print("Running STATIUM with: " + options['IN_CFG'] + " " + options['IN_RES'] + " " + options['IN_PDB'])
        
    elif(options['calc_seq_energy']):
        if(options['-v'] or options['--verbose']): print("Calculating energy for sequence: " + options['SEQ'] + " with STATIUM output directory " + options['DIR'])
        
if __name__ == "__main__":
    main(sys.argv[1:])