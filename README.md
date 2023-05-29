# LearnMemTol
Genetic dissection of learning, memory, and thermal tolerance using the DSPR.

## Instructions for reproducing analyses

### Functions
1. `Functions/mappingfunctions.R` holds some basic functions used for QTL mapping
2. `Functions/h2_functs.R` holds functions for calculating heritability and genetic correlations with a linear model and jackknife method
3. `Functions/ggplot_theme.R` holds the plotting theme for publication format

### Obtain raw data from Zenodo

-Fetch raw data from repository (zenodo)- when done, update any relevant file paths and process raw scripts

### Phenotype analysis steps for the DSPR RILs and Founders: in the `LearnMemTolPheno` folder

1. Process raw data for learning, memory, and thermal tolerance
  - Learning and Memory
      - Processed data from HeatCalc is parsed using a perl script (`heatcalcRead_LearnMem.pl`) and outputs: `ProcessedData/HeatProc_Learn.txt`
        -use `perl myscript.pl > myoutput.txt 2> myerror.txt &` in terminal 
      - Raw `.asc` files are processed using an R script (`process_raw_learnmem.R`) to identify any other issues (see Step 2) and outputs `/ProcessedData/Learn_raw.rda`
  - Thermal tolerance
      - Raw `.asc` files are processed using an R script (`process_raw_therm.R`) to identify issues and calculate thermal tolerance for each individual and outputs `ProcessedData/Incapacitation.rda`. There are matching scripts and output for the data for the founders.

2. Raw data and processed data are checked using `check_raw_data.Rmd`. Notes within this script for inputs. Outputs the following:
  - `ThermalTol_processed.txt`
  - `LearnMem_processed.txt`
These are cleaned datasets ready for next steps. 

3. `Validate_incapacitation.R` compares human scored and automated incapacitation scores from the heatbox. Input is `ProcessedData/ThermalTol_processed.txt` and `ProcessedData/Thermotolerance_test.csv`. and Makes a plot: `Plots/HeatBox_Human_Val.pdf` and outputs the paired data: `ProcessedData/Validation.txt`.  `Compare_heatbox_heatplate.R` compares the two methods for measuring thermal tolerance (box and plate).  

4. `Processing_Pheno_data.Rmd` gets line means and does basic visualizations. Outputs: `ProcessedData/L_MDATA.rda` and `ProcessedData/T_TDATA.rda`. `Processing_ThermTol_Pheno_Founders.Rmd` does this for the founder data and outputs Figure 3b. 

5. `visualize_pheno_data.Rmd` makes several phenotype plots for the RIL data. See script for additional output plots. Outputs Figure 2 and Figure 3c. 

6. Additional script `correlation_between_traits` analyzes the relationship between thermal tolerance, learning, and memory and makes some scatterplots. These are not included in the publication. 

### QTL mapping and peak analysis: in `LearnMemTolQTL` folder

1. `h2_rg.Rmd` calculates heritabilities and genetic correlations
2. `mapping_perms` performs the QTL mapping and permutations and calculation of FDR and FWER, generating several rda files holding the results
2. `identify_peaks` finds each QTL peak and saves it in a list. `IndividualPeaks.Rmd` incorporates additional datasets to make a large summary list. There are notes within the file. The list contains our three phenotypes containing a data.frame with the following columns with one QTL per row:
  - "chr" - the chromosome where the QTL is located
  - "Ppos" - the physical position of the QTL in R5 coordinates
  - "Gpos" - the genetic position of the QTL in cM
  - "Gaxis" - the genetic position on a long scale (chromosomes stacked)
  - "LL" - the LOD score at the peak
  - "lp" - the lower bound physical position in R5 coordinates for the BCI
  - "up" - the upper bound physical position in R5 coordinates for the BCI
  - "lpchr" - the lower bound chromosome for the BCI
  - "upchr" - the upper bound chromosome for the BCI
  - "lg" - the lower bound genetic position for the BCI
  - "ug" - the uppper bound genetic position for the BCI
  - "ulod" - the LOD score at the upper bound
  - "llod" - the LOD score at the lower bound 
  - "chrR6" - the chromosome in R6 coordinates
  - "PposR6" - the physical position of the QTL in R6 coordinates
- "lpR6" the lower bound physical position in R6 coordinates for the BCI
  - "upR6" - the upper bound physical position in R6 coordinates for the BCI
  - "chrN" - the chromosome number without L & R 
"A1" - "A8" - the haplotype effect estimates (haplotype means at the QTL)
  - "A1se" - "A8se" - the standard errors of the haplotype effect estimates
  - "A1N" - "A8N" - the number of RILs with each haplotype genotype at the QTL location   
- "L1" - "L8" - the number of RILs in the "low pool" with each haplotype genotype
- "H1" - "H8" - the number of RILs in the "high pool" with each haplotype genotype 
3.  `QTL_effects_visualization_LM.Rmd` makes the plots of the haplotype effects at each QTL for learning and memory phenotypes. `TT-founder-effect.Rmd` constructs these plots for thermal tolerance data. Makes Figure 3c and Figure S4.
4. `visualize_genome_scan.Rmd` plots the genome scan and makes a composite figure - Figure 3. 
5. `individualpeaks.Rmd` plots each peak separately
6. `Visualize_Peak_Intervals.Rmd` plots each peak interval with the DE genes. Makes figures S6 and S7.
7. `ThermTol_additional_mapping.R` performs a genome scan after correcting for the large Q5 peak and one including subpopulation in the model to compare to our genome scan. Creates figure S2.

### RNAseq analysis: in `LearnMemTolRnaseq` folder

#### The alignment and assembly steps are in `LearnMemTolRnaseq/RNAseq_Processing_LearnMemTol`

1.Set_up_commands_masterpipeline.R is the master script that outputs 8 separate scripts for processing RNAseq data. It sets up scripts to run the RNA-Seq data on a cluster using SLURM. It will require adaptation for your own system and file paths.

Output sripts are named: 

 - S01_Trim_LearnMemRNA.txt: trims data
 - S02_QC_LearnMemRNA.txt: performs quality control
 - S03_Align_LearnMemRNA.txt: aligns data to the fly genome 
 - S04_SamtoBam_LearnMemRNA.txt: converts file from .sam to .bam
 - SO5_MergeBam_LearnMemRNA.txt: merges bam files
 - SO6_Assemble_LearnmemRNA.txt: assembles data

 - S08_Abundances_LearnMemRNA_eb.txt - finds abundances (Stringtie shorter version. This step only identifies genes that are already known. No novel genes are detected.)

or

 - S08_Abundances_LearnMemRNA.txt - finds abundances including novel transcripts


#### Differential Expression Analysis is in `LearnMemTolRnaseq/DEseq_LearnMemTol`

1. `prepDEseq.py` takes the output form S08_Abundances_LearnMemRNA_pl.txt, stringtie (shorter version) and converts it into a usual version for the DEseq pipeline.
2. `LearnMem_RNAproc_DESeq_eb.Rmd` and `ThermTolRNAproc_DEseq_eb.Rmd` runs the DEseq analysis pipeline on the different high and low cohorts for the RNAseq data (learning and memory and thermal tolerance respectively). Additional notes in this script. Output data is in `Data`
3. `Visualize_LMRNA_Data.Rmd` (for learning and memory) and `Visualize_TTRNA_Data.Rmd` (for thermal tolerance) create MA and volcano plots. This script creates Figure 4 and Figure S5.
4. `Visualize_ind_genes.R` creates a series of plots of the normalized counts in each group for individual genes in the Q5 interval. Creates Figure S9.


### RNAi analysis: in `LearnMemTolPheno` folder

1. The `HeatPlate_PreProcess` folder holds the python scripts for all pre processing. This includes: 
a) The script for calibrating the Raspberry Pi camera (`camera_calibration.py`).  
b) The script for collecting the image data during assays (`thermocouple_5.py`). 
c) The `fly_tracker_batch.py` script uses deeplabcut to auto track the flies from the recordings. 
d) The `Flymove.R` script processes the deeplabcut output to make a combined `.csv` file with fly positions over time for each fly. 

2. `process_raw_therm_heatplate.R` and `process_raw_therm_heatplate_val.R` processes the raw output from deeplabcut to score incapacitation for the RNAi lines and the RIL validation set. 

3. `RNAi_post_raw_process.R` checks for and corrects errors and combines the two rounds of phenotyping. 

4. `Models_RNAi.R` performs ANOVAs for each gene and the specific post-doc test for the cross vs the control cross. Outputs: `ProcessedData/RNAi_model_output.csv`

5. `RNAi_viz_individualgenes.R` makes a plot of RNAi results for every focal gene, which is Figure S8. `RNAi_viz.R` makes a plot of several genes together. 

6. `RNAi_DEseq_compare.R` collates the RNAi results and DEseq results to make Figure 5.
