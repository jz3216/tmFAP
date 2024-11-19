import glob
import multiprocessing
import sys
from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.protocols.rosetta_scripts import *
from pyrosetta.rosetta.core.scoring import CA_rmsd
from sys import argv
from rmsd_libs import *
import os

list = argv[1]  # the query and ref pdb list
refpdb = argv[2]  # optional input of the refpdb file from commandline
pocket_res = argv[3]  # optional input of residues

# suppress output
null_device = open(os.devnull, 'w')
sys.stdout = null_device
pyrosetta.init( options="-beta -overwrite -mute all")  # -extra_res_fa /home/jingyi/Projects/HBC/params/599/HBC599_p.params
# Restore stdout
sys.stdout = sys.__stdout__

sele_res_list = pyrosetta.rosetta.std.list_unsigned_long_t()

if pocket_res:
    pocket_res_list = [int(res) for res in pocket_res.split(',')]
else:
    pocket_res_list = [15,16,19,23,51,54,55,58,62,65,100,103,106,107,110,114,135,138,139,142,145,146,149] #514-pocket

for res in pocket_res_list: #514-pocket
    sele_res_list.append(res)

def calc_rmsd(pdb, ref_pdb):
    pose = pose_from_pdb(pdb)
    ref_pose = pose_from_pdb(ref_pdb)
    ca_rmsd_value = CA_rmsd(pose, ref_pose, sele_res_list)
    return f'{pdb}  {f"{ca_rmsd_value:.3g}"}'


with open(list, 'r') as f:
    lines = f.readlines()
    pdbs = [x.split('\n')[0] for x in lines]
    ref_pdbs = [x.split()[1] if len(x.split()) == 2 else refpdb for x in lines ]

    with multiprocessing.Pool() as p:
        rmsds = p.starmap(calc_rmsd, zip(pdbs, ref_pdbs))

    for rmsd in rmsds:
        print(rmsd)