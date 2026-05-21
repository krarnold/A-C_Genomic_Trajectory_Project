## Beta-Binomial GLM and Permutation FDR
Core scripts for the genome-wide SNP analysis using a beta-binomial 
generalized linear model with permutation-based FDR correction.
- GLM_Betabi_AC_Trajectory.R: Main beta-binomial GLM for treatment, generation, and interaction terms
- GLM_Betabi_Permutation.R: Permutation of population labels to generate null distribution for FDR
- GLM_Betabi_Perm_FDR_Plots.R: FDR calculation, effect size filtering, Manhattan plot, and antiparallel hexbin figure
- Reformat_Full_Trajectory_SNP_Table.R: Reformats raw SNP table for GLM input
- Reformat_ACO1_DeepSeq_SNP_table.R: Reformats deep sequencing SNP table for downstream analysis
- Top_SNP_Annotation.R: Annotates top 100 candidate SNPs with genomic features