import argparse
import numpy as np
from colabdesign import mk_afdesign_model, clear_mem

#/storage/lulab/zhujy/599/hal_test/mono_skip.py

argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
####################### i/o args ###########################
argparser.add_argument("--pdb_in", type=str, default='')
argparser.add_argument("--output", type=str, default='')
#argparser.add_argument("--design_config", type=str, default='soft_iters=40, hard_iters=20, num_models=4')
argparser.add_argument("--data_dir", type=str, default='')
argparser.add_argument("--motif", type=str, default='')
#argparser.add_argument("--bias_np", type=str, default='')
argparser.add_argument("--wt_seq", type=str, default='')
argparser.add_argument("--soft_iters", type=int, default='100')
argparser.add_argument("--hard_iters", type=int, default='50')
argparser.add_argument("--rmsd", type=float, default='1')
argparser.add_argument("--num_models", type=int, default='5')
argparser.add_argument("--num_recycles", type=int, default='3')
argparser.add_argument("--non_TM_res_list", type=str, default='')
argparser.add_argument("--RK_ring_res_list", type=str, default='')
argparser.add_argument("--RK_upweight", type=float, default='1.0')
argparser.add_argument("--DE_ring_res_list", type=str, default='')
argparser.add_argument("--DE_upweight", type=float, default='1.0')
argparser.add_argument("--WY_ring_res_list", type=str, default='')
argparser.add_argument("--WY_upweight", type=float, default='1.0')
argparser.add_argument("--chains", type=str, default='A')
argparser.add_argument("--extra_loop", type=str, default='')
argparser.add_argument("--intra_loop", type=str, default='')
argparser.add_argument("--ramp_models", action='store_true')
argparser.add_argument("--use_multimer", action='store_true')
argparser.add_argument("--monomer", action='store_true')
argparser.add_argument("--apo_loop", type=str, default='')
#argparser.add_argument("--become_apo_res", type=str, default='')
args = argparser.parse_args()

# default allow input oligomers unless --monomer is set
if args.monomer:
    n_chain=1
    chains='A'
    args.use_multimer=False
    h_olig=False
else:
    h_olig=True

aa_order="ARNDCQEGHILKMFPSTWYV"
seq_len=len(args.wt_seq)    #161
n_chain=len(args.chains.split(','))
tot_len=seq_len*n_chain
wt_seq=args.wt_seq   
chains=args.chains

#all arrays are 0-indexed
#fix_res    # renamed to non_TM_res_list for clarity
if args.non_TM_res_list:
    fix_res=np.array([int(x) for x in args.non_TM_res_list.split(',') ])-1
else:
    fix_res=[]
design_res=[int(x) for x in range(seq_len) if x not in fix_res] 
#RK ring
if args.RK_ring_res_list:
    RK_ring_res=np.array([int(x) for x in args.RK_ring_res_list.split(',') ])-1
else:
    RK_ring_res=[]
if args.WY_ring_res_list:
    WY_ring_res=np.array([int(x) for x in args.WY_ring_res_list.split(',') ])-1
else:
    WY_ring_res=[]
if args.DE_ring_res_list:
    DE_ring_res=np.array([int(x) for x in args.DE_ring_res_list.split(',') ])-1
else:
    DE_ring_res=[]

#design loop seq
if args.extra_loop:
    extra_loop_res=np.array([int(x) for x in args.extra_loop.split(',') ])-1
else:
    extra_loop_res=[]

if args.intra_loop:
    intra_loop_res=np.array([int(x) for x in args.intra_loop.split(',') ])-1
else:
    intra_loop_res=[]
if args.apo_loop:
    apo_loop_res=np.array([int(x) for x in args.apo_loop.split(',') ])-1
else:    
    apo_loop_res=[]



#callbacks

def fixed_res(self, fix_res=fix_res, wt_seq=wt_seq):    
    #apply fix_res
    for resi in fix_res:
        aa = wt_seq[resi]
        aa_i = aa_order.index(aa)
        self._inputs['bias'][resi, :] = -1e8    #strict fix?
        self._inputs['bias'][resi, aa_i] = 1e8

def rings_(self, RK_ring_res=RK_ring_res, WY_ring_res=WY_ring_res, DE_ring_res=DE_ring_res):

    for aa in 'RK':
        aa_i = aa_order.index(aa)
        self._inputs['bias'][RK_ring_res, aa_i] = 0 #allow design
        #upweight RK at RK ring 
        self.aux['grad']['seq'][0, RK_ring_res, aa_i] *= args.RK_upweight
        #disbale RK at DE ring
        self._inputs['bias'][DE_ring_res, aa_i] = -1e8
    
    if len(WY_ring_res) > 0:
        for aa in 'WY':
            aa_i = aa_order.index(aa)
            self._inputs['bias'][WY_ring_res, aa_i] = 0 #allow design
            #upweight WY at WY ring
            self.aux['grad']['seq'][0, WY_ring_res, aa_i] *= args.WY_upweight
            #disbale WY at RK ring
            #self._inputs['bias'][RK_ring_res, aa_i] = -1e8

    for aa in 'DE':
        aa_i = aa_order.index(aa)
        self._inputs['bias'][DE_ring_res, aa_i] = 0 #allow design
        #upweight DE at DE ring
        self.aux['grad']['seq'][0, DE_ring_res, aa_i] *= args.DE_upweight
        #disbale DE at RK ring
        self._inputs['bias'][RK_ring_res, aa_i] = -1e8

    #design polar res at rings
    for aa in 'ASTNQ':
        aa_i = aa_order.index(aa)
        self._inputs['bias'][RK_ring_res, aa_i] = 0
        self._inputs['bias'][DE_ring_res, aa_i] = 0
        if len(WY_ring_res) > 0:
            self._inputs['bias'][WY_ring_res, aa_i] = 0


'''    #disable AST at rings
    for aa in 'AST':
        aa_i = aa_order.index(aa)
        self._inputs['bias'][RK_ring_res, aa_i] = -1e8
        self._inputs['bias'][WY_ring_res, aa_i] = -1e8'''

        
def design_TMspan(self, design_res=design_res):
    if len(design_res) > 0:
        #clear bias
        for resi in design_res:
            self._inputs['bias'][resi, :] = 0
        #disable polar
        for aa in 'CPSTNQHDERKWY':
            aa_i = aa_order.index(aa)
            self._inputs['bias'][design_res, aa_i] = -1e8
        #reduce MG
        for aa in 'MG':
            aa_i = aa_order.index(aa)
            self.aux['grad']['seq'][0, design_res, aa_i] *= 0.1

def extra_loop(self, extra_loop_res=extra_loop_res):
    if len(extra_loop_res) > 0:
        for resi in extra_loop_res:
            self._inputs['bias'][resi, :] = 0
        for aa in aa_order:
            if aa not in 'ADEFGHILNPQSTV':     #classic loop aas but no RK
                aa_i = aa_order.index(aa)
                self._inputs['bias'][extra_loop_res, aa_i] = -1e8

def intra_loop(self, intra_loop_res=intra_loop_res):
    if len(intra_loop_res) > 0:
        for resi in intra_loop_res:
            self._inputs['bias'][resi, :] = 0
        for aa in aa_order:
            if aa not in 'ARKFGHILNPQSTV':     #classic loop aas but no DE
                aa_i = aa_order.index(aa)
                self._inputs['bias'][intra_loop_res, aa_i] = -1e8

def apo_loop(self, apo_loop_res=apo_loop_res):
    if len(apo_loop_res) > 0:
        for resi in apo_loop_res:
            self._inputs['bias'][resi, :] = 0
        for aa in aa_order:
            if aa not in 'AFGILPVWY':     #tm res + GPWY
                aa_i = aa_order.index(aa)
                self._inputs['bias'][apo_loop_res, aa_i] = -1e8
        

post_design_callback=[] #[fixed_res, design_TMspan, rings_, extra_loop, intra_loop, apo_loop]
post_design_callback.append(fixed_res)
post_design_callback.append(design_TMspan)
post_design_callback.append(rings_)
if args.extra_loop:
    post_design_callback.append(extra_loop)
if args.intra_loop:
    post_design_callback.append(intra_loop)
if args.apo_loop:
    post_design_callback.append(apo_loop)



def main():
    #init af2model
    clear_mem()
    #print("fix_res: ", fix_res)
    model = mk_afdesign_model(protocol="partial", post_design_callback=post_design_callback, data_dir=args.data_dir, use_multimer=args.use_multimer)   #  , rings, become_TM

    motif=args.motif    #   15-16,19,23,26,30,51,54-55,58,62,65,100,103,106-107,110,114,131,135,138-139,142,145-146,149
    model.set_opt(num_recycles=args.num_recycles)
    model.prep_inputs(pdb_filename=args.pdb_in, chain=args.chains, homooligomer=h_olig, pos=motif, fix_seq=True, length=tot_len)  #   chains="A,B"

    #print(model.opt["pos"])

    #allow AST
    #disable MG 
    #model.set_seq(rm_aa="C,P", bias=np.zeros((seq_len,20)))
    model.set_seq(rm_aa="C,P,D,E,N,Q,R,K,S,T,H,Y,W", bias=np.zeros((seq_len,20))) #init bias with blank array. 
    # rm_aa is applied after this.
    model.set_weights(dgram_cce=1, rmsd=args.rmsd, con=0, plddt=1)   #dgram_cce=1, rmsd=2

    #print settings
    #print(f'chains: {chains}, n_chain: {n_chain}, monomer? {args.monomer},  use_multimer_model? {args.use_multimer}, \npost_design_callback: {post_design_callback} \n fix_res: {fix_res} \n motif res: {motif} \n RK_ring_res: {RK_ring_res} \n WY_ring_res: {WY_ring_res} \n DE_ring_res: {DE_ring_res} \n design_res: {design_res} \n extra_loop_res: {extra_loop_res} \n intra_loop_res: {intra_loop_res} \n apo_loop_res: {apo_loop_res}')

    #run
    model.design_pssm_semigreedy(soft_iters=args.soft_iters, hard_iters=args.hard_iters, ramp_models=args.ramp_models, num_models=args.num_models, models=[0,1,2,3,4])

    #output
    model.save_pdb(args.output)

if __name__ == "__main__":

    main()
