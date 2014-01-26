import sys
from docopt import docopt
from statium_reformat import renumber
from statium_reformat import create_res
from statium_analysis import statium_pipeline

def main(argv):
    
    helpdoc =   """
                usage: statium_wrapper.py renumber (IN_PDB) [OUT_PDB] [-v | --verbose]
                       statium_wrapper.py create_res (IN_PDB_ORIG IN_PDB_RENUMBERED) [OUT_RES] [-v | --verbose]
                       statium_wrapper.py run_statium (IN_RES IN_PDB IN_PDB_LIB_DIR IN_IP_LIB_DIR) [OUT_DIR] [-v | --verbose]
                       statium_wrapper.py calc_seq_energy (DIR SEQ) [-v | --verbose]
                       statium_wrapper.py [-h | --help]
                """
    
    options = docopt(helpdoc, argv, help = True, version = "2.0.0", options_first=False)
    verbose = options['-v'] or options['--verbose']
    
    if(options['renumber']):
        if(verbose): print("Reformatting file: " + options['IN_PDB'])
        
        if(options['OUT_PDB'] == None):
            renumber(1, 1, options['IN_PDB'], options['IN_PDB'][:-4]+'_renumbered.pdb')
            if(verbose): print("Done. Formatted file: " + options['IN_PDB'][:-4]+'_renumbered.pdb')
            
        else:
            renumber(1, 1, options['IN_PDB'], options['OUT_PDB'])
            if(verbose): print("Done. Formatted file: " + options['OUT_PDB'])
    
        
    elif(options['create_res']):
        if(verbose): print("Creating .res file using: " + options['IN_PDB_ORIG'] + " and " + options['IN_PDB_RENUMBERED'])
        
        if(options['OUT_RES'] == None):
            create_res(options['IN_PDB_ORIG'], options['IN_PDB_RENUMBERED'], options['IN_PDB_RENUMBERED'][:-4]+'_renumbered.pdb')
            if(verbose): print("Done. .res file: " + options['IN_PDB_RENUMBERED'][:-4]+'.res')
            
        else:
            create_res(options['IN_PDB_ORIG'], options['IN_PDB_RENUMBERED'], options['OUT_RES'])
            if(verbose): print("Done. .res file: " + options['OUT_RES'])
    
    elif(options['run_statium']):
        
        if(verbose): print("Running STATIUM with: " + options['IN_RES'] + " " + options['IN_PDB'] + " " + options['IN_PDB_LIB_DIR'] + " " + options['IN_IP_LIB_DIR'])
        
        if(options['OUT_DIR'] == None):   
            statium_pipeline(options['IN_RES'], options['IN_PDB'], options['IN_PDB_LIB_DIR'], options['IN_IP_LIB_DIR'], options['IN_RES'][:-4], verbose)
            if(verbose): print("Done. STATIUM probabilities in output directory: " + options['IN_RES'][:-4]);
            
        else:
            statium_pipeline(options['IN_RES'], options['IN_PDB'], options['IN_PDB_LIB_DIR'], options['IN_IP_LIB_DIR'], options['OUT_DIR'], verbose)
            if(verbose): print("Done. STATIUM probabilities in output directory: " + options['OUT_DIR']);

    elif(options['calc_seq_energy']):
        if(verbose): print("Calculating energy for sequence: " + options['SEQ'] + " with STATIUM output directory " + options['DIR'])
        
if __name__ == "__main__":
    main(sys.argv[1:])