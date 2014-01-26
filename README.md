statium!
=======

Implementation of the STATIUM algorithm [protein-protein binding affinity scoring]

                usage: statium_wrapper.py renumber (IN_PDB) [OUT_PDB] [-v | --verbose]
                       statium_wrapper.py create_res (IN_PDB_ORIG IN_PDB_RENUMBERED) [OUT_RES] [-v | --verbose]
                       statium_wrapper.py run_statium (IN_RES IN_PDB IN_PDB_LIB_DIR IN_IP_LIB_DIR) [OUT_DIR] [-v | --verbose]
                       statium_wrapper.py calc_seq_energy (IN_RES IN_DIR PROBS_DIR SEQ) [-v | --verbose]
                       statium_wrapper.py [-h | --help]

An example of an analysis sequence: <br>
&nbsp;&nbsp;&nbsp;&nbsp;``python statium_wrapper.py renumber testing/4hfz_orig.pdb testing/4hfz.pdb -v``<br>
&nbsp;&nbsp;&nbsp;&nbsp;``python statium_wrapper.py create_res testing/4hfz_orig.pdb testing/4hfz.pdb testing/4hfz.res -v``<br>
&nbsp;&nbsp;&nbsp;&nbsp;``python statium_wrapper.py run_statium testing/4hfz.res testing/4hfz.pdb culled_90_reduced/ ip_90_wGLY/ testing/output -v`` <br>
&nbsp;&nbsp;&nbsp;&nbsp;``python statium_wrapper.py calc_seq_energy testing/4hfz.res testing/output/ testing/output_probs/ ETFSDLWKLLP``<br>
