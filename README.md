# HGT_identification_in_evol_exps
Scripts for HGT identification workflow.

This repository includes Python and R scripts used in the analysis of an evolution experiment with HGT. The project is in active development, and at the moment, the workflow is not fully automated. For that reason, though the entire workflow makes use of the Breseq pipeline for analysis of microbial populations (Deatherage and Barrick 2014), Breseq is not explicity called here. 

Read files available upon request or at BioProject ID PRJNA720176.

1. Optional: Trim and filter all reads using bbduk (https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/) to desired quality and minimum read length. Include the 'ref' option to remove reads aligning to a specific contaminant.

2. Install and run breseq (Deatherage and Barrick 2014) in polymorphic mode (-p flag) for all samples, including the recipient and donor(s), with the appropriate recipient reference sequence (Figure 1a). Users have reported issues running Breseq on Linux CentOS, so if installing with Anaconda, it is recommended that a seperate environment is created for breseq which is activated and deactivated as needed. To detect all variants, ensure that the polymorphism frequency threshold is decreased to the user-preferred minimum. When selecting a cutoff, note that the minimum threshold is applied to both the reference and alternate versions of a polymorphism, so if the threshold is 0.05, variants that occur at 96% will be called as fixed. Note also that variants can be filtered by conditions including frequency at a later stage using gdtools FILTER.

3. Retrive the set of variants that correspond to HGT events with the gdtools suite of modules that accompany breseq using gdtools SUBTRACT and COMPARE (see Figure 1b). Use of gdtools NORMALIZE can allow better comparisons across samples.

4. If there is little divergence between the donor and recipient strain, and the donor is expected to map well to the recipient reference genome, the analysis can end here. If required for downstream analyses, HGT, or class of variants can be applied to the reference genome using gdtools APPLY (Figure 1c). Note that only fixed variants will be applied in this case, so it is possible to apply all variants (regardless of frequency) if frequency values are replaced with 1 (Eg. using a simple sed command). The reference genome with variants applied for each population will be known as that population's "foundational HGT genome."

![image](https://user-images.githubusercontent.com/46662443/123540531-cf62a180-d782-11eb-95fb-625ad4f79450.png)
 
5. If divergence between the donor and receipient is high and especially if the user notices that a high proportion of HGT-evolved reads are not mapping to the recipient, one may want to check for HGT-derived variants that are not able to be mapped to the recipient genome. One solution is to assemble the discarded donor reads from the reference-guided assembly into genomic fragments or scaffolds against which evolved reads can be mapped. Here, I used SPAdes (Nurk, Bankevich et al., 2013).

6. If a full genomic assembly for the donor is available, the donor scaffolds can be checked against it using Blast. Scaffolds may also be filtered for length to eliminate noise. Use the makeblastdb function from the BLAST+ executables (Camacho, 2009) to make a database from recipient and donor(s) genomes. 

7. Running the BLASTN algorithm is relatively simple but a script using the Biopython (Cock et al. 2009) Blast wrapper is given as an example (See Step 7).

8. If a full genome is not available, the user may have to prune sequences with matches to the recipient or to contaminants. Contaminants can often be identified by their low coverage compared to scaffold length in a _de novo_ assembly, though high levels of contamination may make this more difficult.

9. Concatenate the donor scaffolds formed from the discarded donor reads and the recipient reference sequence into a single reference. Run breseq for each of the evolved samples against this new concatenated reference.

10. Apply the variants to the concatenated reference using gdtools APPLY. The aim here is to try and "claw back" some of reference sequence on the donor scaffolds where there is short-read coverage. In some cases, this will effectively transform the donor scaffold back into the original reference, meaning that the reads mapping there were mismapped and are now accounted for. In other cases, this may help to identify recombination breakpoints and make it easier to determine where the donor DNA has recombined.

11. Now, excise the regions of the donor scaffolds that have short read coverage using the script indicated as Step 11.

12. Here, we filter the donor coverage regions for all populations. Blast all coverage regions against the reference and donor genome sequences to remove any scaffolds that simply match the reference (and are therefore the result of read mismapping to the donor coverage regions). See Step 12 for an example Blast script using the BioPython Blast wrapper. Most of the donor coverage regions for the non-HGT controls will be removed at this stage, as will some of the donor coverage regions in the HGT populations, but there should be marked difference between the two. If not, consider whether DNA contamination of the non-HGT controls with donor DNA is a possibility.

13. Now, filter any donor coverage regions remaining in the non-HGT control populations from the HGT populations. Again, this should remove a very small proportion of the coverage regions identified in the HGT populations.

14. For each population, use a long-read mapper to map the donor coverage regions to the foundational HGT genomes and apply the read varinats to produce a consensus genome. Using the foundational HGT genome may improve the ability to map the donor coverager regions through overlapping HGT variants between the two methods. It also merges the two sets of variants identified into a single genome known as the "maxHGT genome" for each population. This genome is likely to have HGT variants that ultimately do not have short read coverage, but it represents a putative sequence against which to map the population short reads.

In the next few steps, the aim is to "bin" short reads to either the reference genome or the maxHGT genome. Because reads representing recombination breakpoints may map ambiguously, we want to merge the reference vs HGT variants for both bins. However, these positions may differ between the two bins due to the application of indels to the maxHGT genome sequence. The next few steps shows how to solve this.

15. Use breseq (or any aligner) to map short reads to the reference and maxHGT genome sequences simultaneously.

16. Use a pairwise genome alignment tool to get recipricol variants calls--that is, each variants is called in reference to both genome sequences. I use progressiveMauve and exproted SNPs and gaps, and the scripts here assume progressiveMauve outputs.

17. If progressiveMauve was used, use the script indicated at Step 17 to obtain the allele frequencies at each site using the information from both bins. Input the bam alignments from the alignment in Step 15. This script will output gd-formatted variant calls and a csv file with absolute read cocverage valuess for each putative HGT variant. The former is appropriate for use with gdtools while the latter is useful for input to the R glm script for determining selection coefficients by comparing relative variant proportions over time.
![image](https://user-images.githubusercontent.com/46662443/124859007-6d136780-dff2-11eb-8c15-83cb3151c76d.png)

18. Many low-frequency variants may be identified, so the user can filter the gd file for variants above specific thresholds. Eg. 1%+, 50%+ and 80%+. If these variants are then applied to the reference sequence, the user can examine, for example, only the effects of high frequency HGT alleles.
