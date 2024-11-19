# wFAP

scripts used for generating the reported wFAPs

## Required softwares 

- [Rosetta](https://www.rosettacommons.org)
    - Rosetta Modeling Suit version 2019.47.61047 (this work) or other verions of this software are available from https://www.rosettacommons.org
- [RifDock](https://github.com/rifdock/rifdock)
- [ColabFold](https://github.com/sokrypton/ColabFold)

## Command line examples

### 1. parametric scaffold generation

input parameters for Crick equation for each helices:

- `r0_1, r0_2`      - r0 (major helix radius)   
- `dw1_1, dw1_2`    - Δω1 (rotation about the helix's own axis) 
- `zoff`            - z0 offset (displacement along the helical bundle's z-axis) 

```
$ROSETTA @bb.flag -parser:protocol c2.xml -parser:script_vars r0_1=${r0_1} r0_2=${r0_2} dw1_1=${dw1_1} dw1_2=${dw1_2} zoff=${zoff}
```



### 2. rifgen

rifgen takes the following as input:  
- `input/HBC599_O_conf${conf}.pdb`      - HBC599 pdb, with -CN replaced by -C=O for compatibility to rifdock (only used for rifgen/rifdock).  
- `input/599_O_p.params`                - HBC599 params file, with the above mentioned substitution  
- `input/tuning.x`                      - controls requirements for hbond search. Here we required negative charges (D/E) h-bond to the -OH in HBC599  

```
$rifgen @rifgen_599_O.flag -database $database_path -rifgen:outdir $outpath -rifgen:target input/HBC599_O_conf${conf}.pdb -extra_res_fa input/599_O_p.params -tuning_file input/tuning.x >> $log
```

- Each conformation of the target ligand requires a seperate run of rifgen, and is suffixed with ${conf}. 

### 3.rifdock

- To run rifdock, add the path of rifgen output to `rifdock_.flag` (available at the end of rifdock output).  
- `$scaffold_pdb` is the parametric scaffold generation output in the case of initial search. For iterative optimization, the protein chain of Rosetta designed poses are used as the input scaffold.

```
$rif_dock_test  @rifdock_.flag  -rif_dock:scaffolds $scaffold_pdb -rif_dock:outdir ${outpath} -rif_dock:output_tag conf${conf_i} >> $rifdock_log
```

### 4. sequence initilization by Rosetta fix-backbone design

This step initializes sequences of parametric helical bundles after rifdock. Both RIF residues and pre-defined hydrogen bond network residues are labeled as 'HBNet' residues, to indicate they are fixed during Rosetta design.

```
$ROSETTA @design.flag -parser:protocol pack.xml  -s $pdb_in -in:file:extra_res_fa HBC599_p.params -out:path:all $outpath
```

### 5. Rosetta sequence design ('update sequence' in iteration)

Required input:  
- `input/hbond.cst`     - constraint files defining hbonds found by rifdock  
- `input/all.comp`      - amino acid composition adjustment for Rosetta sequence design

```
$ROSETTA @design.flag -parser:protocol design_bbcst.xml -s $pdb_in -in:file:extra_res_fa HBC599_p.params c -parser:script_vars cst_file='input/hbond.cst' aacom_file='input/all.comp' -out:path:all $outpath
```

### 6. AF2 inference using ColabFold ('update scaffold' in iteration)

```
colabfold_batch --rank plddt  --msa-mode single_sequence --overwrite-existing-results  --num-recycle 3  --use-gpu-relax  --amber $fasta $outpath
```

### 7. filter by CA rmsd and plddt


```
python ca_rmsd_p.py $pdb_list $ref.pdb $pocket_res_list >> rmsd
```
pdb files (i.e. AF2 prediction models) in `$pdb_list` are superimposed to `$ref.pdb` (design model) by aligning residues in the `$pocket_res_list`

