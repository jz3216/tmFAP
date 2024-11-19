# transmembrane span hallucination and tmFAPs

scripts used for tmFAP design and transmembrane span hallucination in general

## note

For making new designs, it is recommended to use the newer script `hal_tmspan.py` in the section `transmembrane span halluciantion in general`. The script `hal_tmFAPs.py` in the `design of tmFAPs` section was used with ColabDesign v1.0.8 to generate tmFAPs, but is not compatible with more recent ColabDesign versions. It's functions have been covered by the updated scripts, and we only put it here for archiving purpose.

## Required softwares 

[ColabDesign](https://github.com/sokrypton/ColabDesign)  v1.1.1

## transmembrane span hallucination

This script converts water-soluble proteins to transmembrane proteins with minimum perturbance to the functional motif structure.

- in an environment with colabdesign installed, run:  

```
python hal_tmspan.py \
        --pdb_in            $pdb_in \
        --soft_iters        $soft_iters \
        --hard_iters        $hard_iters \
        --wt_seq            $wt_seq \
        --output            $output \
        --data_dir          $data_dir \
        --num_recycles      $num_recycles \
        --non_TM_res_list   $non_TM_res_list \
        --motif             $motif \
        --chains            $chains    \
        --RK_ring_res_list  $RK_ring_res_list \
        --RK_upweight       $RK_upweight      \
        --DE_ring_res_list  $DE_ring_res_list \
        --DE_upweight       $DE_upweight      \
        --extra_loop        $extra_loop \
        --intra_loop        $intra_loop \
        --use_multimer 
```
### input and initialization

This script is compatible with both monomer and homo-oligomer input. You may input a model with or without transmembrane sequence initialization. For better success rate, start with a good input i.e. the structure agree with AF2 single seqeunce predictions, at least for the motif region.

### output 

In our test cases, successful designs converge to plddt > 85 and rmsd < 1. You may additionally validate the output using other strucuture prediction methods.

### the flags

- `$soft_iters`         - number of steps for relaxed sequence-space design  
- `$soft_iters`         - number of steps for relaxed sequence-space design  
        - The steps required are highly dependent on the input target. Start with `soft_iters=100` and `hard_iters=50` and increase them if not converging.  
- `wt_seq`              - typically fasta sequence of $pdb_in  e.g. "DEERLKEILFFLLLIIIFVVFLLIVHYKFLEEFKEKNV"  Used to set amino acid identity at the fixed positions. 
- `output`              - the output pdb file  
- `data_dir`            - path to AF2 weights  
- `non_TM_res_list`     - string of residue indices of non-TMspan residues. e.g. '1,2,3'  
- `motif`               - string of residue indices of pocket/functional site residues.  e.g. "15-16,19,23,51,54-55,58,62,65,100,103,106-107,110,114,135,138-139,142,145-146,149"  
- `chains`              - string of chain IDs. e.g. 'A,B,C,D'  
- `RK_ring_res_list`    - string of residue indices of RK ring residues. Remove this flag and `RK_upweight` if you don't want the RK ring.  e.g. '1,2,3'  
- `RK_upweight`         - increase RK. 1 is allow design of RK but no additional bias. Start with 10 and increase if you need more RK.  
- `DE_ring_res_list`    - Similar to `RK_ring_res_list`  
- `DE_upweight`         - Similar to `RK_upweight`  
- `extra_loop`          - allow sequence redesign of extracellualr polar loops  
- `intra_loop`          - allow sequence redesign of intracellualr polar loops  
- `monomer`             - use this flag if designing monomer  
- `use_multimer`        - use AF2 multimer weights. Use this flag only if designing multimer.  


---

## design of tmFAPs

Scripts and input used for generating the tmFAPs in the manuscript

### 1. transmembrane span definition and sequence initilization

Scripts and usage for this step is the same as described in the previous work  

> 'Xu, C., Lu, P., Gamal El-Din, T.M. et al. Computational design of transmembrane pores. Nature 585, 129–134 (2020).'

### 2. transmembrane span sequence hallucination

```
python hal_tmFAPs.py \
        --pdb_in            $pdb_in --output $output\
        --soft_iters        $soft_iters --hard_iters $hard_iters \
        --wt_seq            $wt_seq \
        --rmsd              $rmsd \
        --data_dir          $path_to_AF2_params \
        --num_recycles      $num_recycles \
        --fix_res_list      $fix_res_list \
        --motif             $motif \
        --wy_ring_res_list  $wy_ring_res_list\
        --rk_ring_res_list  $rk_ring_res_list
```

Runnning this script requires the output pdb and the TMspan definition (residues indices for each region) from the initialization step.

- `pdb_in`                 - the initialized pdb  
- `soft_iters`             - number of steps for relaxed sequence-space design  
- `hard_iters`             - number of steps for one-hot sequence optimization  
- `wt_seq`                 - fasta sequence of $pdb_in e.g. "DEERLKEILFFLLLIIIFVVFLLIVHYKFLEEFKEKNV"  
- `rmsd`                   - weights for rmsd loss  
- `motif`                  - residue indices of pocket residues (used for pocket deviation loss calculation) e.g. "15-16,19,23,51,54-55,58,62,65,100,103,106-107,110,114,135,138-139,142,145-146,149"  
- `wy_ring_res_list`       - string of residue indices of WY ring residues. e.g. '27,28,50,52,111,112'  
- `rk_ring_res_list`       - string of residue indices of RK ring residues  
- `fix_res_list`           - string of residue indices of non-TMspan residues  


