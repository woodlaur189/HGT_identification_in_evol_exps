# HGT_identification_in_evol_exps
Scripts for HGT identification workflow.

This repository includes Python and R scripts used in the analysis of an evolution experiment with HGT. The project is in active development, and at the moment, the workflow is not fully automated. For that reason, though the entire workflow makes use of the Breseq pipeline for analysis of microbial populations (Deatherage and Barrick 2014), Breseq is not explicity called here. 

Read files available upon request or at BioProject ID PRJNA720176.

1. Optional: Trim and filter all reads using bbduk () to desired quality and minimum read length. Eg. 
Include the 'ref' option to remove reads aligning to a specific contaminant.

2. Install and run breseq () in polymorphic mode (-p flag) for all samples, including the recipient and donor(s), with the appropriate recipient reference sequence (Figure 1a). Users have reported issues running Breseq on Linux CentOS, so if installing with Anaconda, it is recommended that a seperate environment is created for breseq which is activated and deactivated as needed. To detect all variants, ensure that the polymorphism frequency threshold is decreased to the user-preferred minimum. When selecting a cutoff, note that the minimum threshold is applied to both the reference and alternate versions of a polymorphism, so if the threshold is 0.05, variants that occur at 96% will be called as fixed. Note also that variants can be filtered by conditions including frequency at a later stage using gdtools FILTER.

3. Retrive the set of variants that correspond to HGT events with the gdtools suite of modules that accompany breseq using gdtools SUBTRACT and COMPARE (see Figure 1b). Use of gdtools NORMALIZE can allow better comparisons across samples.

4. If there is little divergence between the donor and recipient strain, and the donor is expected to map well to the recipient reference genome, the analysis can end here. If required for downstream analyses, HGT, or class of variants can be applied to the reference genome using gdtools APPLY (Figure 1c). Note that only fixed variants will be applied in this case, so it is possible to apply all variants (regardless of frequency) if frequency values are replaced with 1.

![image](https://user-images.githubusercontent.com/46662443/123540531-cf62a180-d782-11eb-95fb-625ad4f79450.png)
 
5. If divergence between the donor and receipient is high and especially if the user notices that a high proportion of HGT-evolved reads are not mapping to the recipient, one may want to check for HGT-derived variants that are not able to be mapped to the recipient genome. One solution is to assemble the discarded donor reads from the reference-guided assembly into genomic fragments or scaffolds against which evolved reads can be mapped. Here, I used SPAdes ().

6. If a full genomic assembly for the donor is available, the donor scaffolds can be checked against it using Blast. Scaffolds may also be filtered for length to eliminate noise. Use the makeblastdb function from the BLAST+ executables () to make a database from recipient and donor(s) genomes. 

7. Running the BLASTN algorithm is relatively simple but a script using the Biopython Blast wrapper is given as an example (See Step 7).

8. If a full genome is not available, the user may have to prune sequences with matches to the recipient or to contaminants. Contaminants can often be identified by their low coverage compared to scaffold length in a _de novo_ assembly, though high levels of contamination may make this more difficult.

9. Concatenate the donor scaffolds formed from the discarded donor reads and the recipient reference sequence into a single reference. Run breseq for each of the evolved samples against this new concatenated reference.

10. 

