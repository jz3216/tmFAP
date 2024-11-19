import argparse
import numpy as np
from colabdesign import mk_afdesign_model, clear_mem

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

####################### i/o args ###########################
argparser.add_argument("--pdb_in", type=str, default='')
argparser.add_argument("--output", type=str, default='')
#argparser.add_argument("--design_config", type=str, default='soft_iters=40, hard_iters=20, num_models=4')
argparser.add_argument("--data_dir", type=str, default='')
argparser.add_argument("--motif", type=str, default='')
argparser.add_argument("--bias_np", type=str, default='')
argparser.add_argument("--wt_seq", type=str, default='')
argparser.add_argument("--soft_iters", type=int, default='100')
argparser.add_argument("--hard_iters", type=int, default='50')
argparser.add_argument("--rmsd", type=float, default='1.0')
argparser.add_argument("--WY_upweight", type=float, default='1.0')
argparser.add_argument("--RK_upweight", type=float, default='1.0')
argparser.add_argument("--num_models", type=int, default='4')
argparser.add_argument("--num_recycles", type=int, default='0')
argparser.add_argument("--fix_res_list", type=str, default='')
argparser.add_argument("--wy_ring_res_list", type=str, default='')
argparser.add_argument("--rk_ring_res_list", type=str, default='')
args = argparser.parse_args()

#design restrictions
aa_order="ARNDCQEGHILKMFPSTWYV"
seq_len=len(args.wt_seq)
fix_res_list=[int(x) for x in args.fix_res_list.split(',')]     #1-indexed    
fix_res=np.array(fix_res_list)-1    

wy_ring_res_list=[int(x) for x in args.wy_ring_res_list.split(',')]  #args.wy_ring_res
wy_ring_res=np.array(wy_ring_res_list)-1
rk_ring_res_list=[int(x) for x in args.rk_ring_res_list.split(',')] 
rk_ring_res=np.array(rk_ring_res_list)-1

wt_seq=args.wt_seq  

#lower M,G grad
def mod_grad(self):    
    for aa in ['M','G']:
        aa_i = aa_order.index(aa)
        self.aux['grad']['seq'][0,:,aa_i] *= 0.1

def add_WYring(self, wy_ring_res=wy_ring_res):    
    for aa in ['W','Y']:
        aa_i = aa_order.index(aa)
        self.opt['bias'][wy_ring_res, aa_i] = 0
        self.aux['grad']['seq'][0, wy_ring_res, aa_i] *= args.WY_upweight
        
def add_RKring(self, rk_ring_res=rk_ring_res):    
    for aa in ['W','Y']:
        aa_i = aa_order.index(aa)
        self.opt['bias'][rk_ring_res, aa_i] = 0
        self.aux['grad']['seq'][0, rk_ring_res, aa_i] *= args.RK_upweight

def main():
    clear_mem()
    model = mk_afdesign_model(protocol="partial", design_callback=[add_WYring, add_RKring], data_dir=args.data_dir)  #, add_RKring
    motif = args.motif
    model.set_opt(num_recycles=args.num_recycles)
    model.prep_inputs(pdb_filename=args.pdb_in, chain="A", pos=motif, fix_seq=True, length=seq_len)
    model.set_seq(rm_aa="C,P,D,E,N,Q,R,K,S,T,H,Y,W,M,G", bias=np.zeros((seq_len,20))) 
    model.set_weights(dgram_cce=1, rmsd=args.rmsd, con=0, plddt=1) 

    #run
    model.design_pssm_semigreedy(soft_iters=args.soft_iters, hard_iters=args.hard_iters, num_models=args.num_models, models=[0,1,2,3])  #ramp_models=False, 

    #output
    model.save_pdb(args.output)

if __name__ == "__main__":
    main()
