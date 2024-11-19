from Bio.PDB import *

def load_plddt(pdb_file):

    parser = PDBParser(PERMISSIVE=True, QUIET=True)
    data = parser.get_structure("A", pdb_file)
    model = data[0]
    chainA = model["A"]

    sum = 0
    for res in chainA:
        res_plddt = res["CA"].get_bfactor()
        sum += res_plddt

    plddt = sum/len(chainA)

    return plddt