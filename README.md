# A-C Genomic Trajectory Project

Scripts used for the manuscript "Beyond Fixation: Persistent Genetic 
Variation Under Intense Selection" (10.64898/2026.03.02.706684).

## Repository Structure

- **AFS/**: Allele frequency spectrum analysis and deep sequencing scripts
- **BetaBinomial_ME/**: Beta-binomial GLM and permutation-based FDR analysis
- **GLMM/**: Old vs new population comparisons using permutation-based GLM
- **Heterozygosity/**: Genome-wide heterozygosity computation and visualization
- **PCA/**: Principal components analysis and effective population size estimation
- **Phenotypic/**: Larval development and age-specific mortality analyses
- **LOO_parallelism_test/**: Deprecated — see DEPRECATED.txt

## Data Availability

Major input and output files are available through Dryad 
(DOI: 10.5061/dryad.tdz08kqdc).

Raw sequencing data are available through NCBI SRA 
(BioProject: PRJNA1426013).

## Updates

Scripts updated for resubmission to Molecular Biology and Evolution 
(May 2026). Major changes include: permutation-based FDR framework 
replacing Benjamini-Hochberg and Bonferroni corrections, revised PCA 
with mean centering and no scaling, and updated deep sequencing 
analysis with minimum minor allele count threshold of 3.

## Contact

Kenneth Arnold (Kenneth.Arnold@oregonstate.edu)
Mark Phillips (philmark@oregonstate.edu)
