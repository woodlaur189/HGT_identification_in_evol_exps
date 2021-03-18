#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 23:49:40 2020

@author: lwoo0005
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 13:55:01 2020

@author: lwoo0005
"""

#Script must be run through conda environment with pysam installed
#(interferes with breseq)
#use conda activate pysam_env

import pysam
from pysam import bcftools
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from string import Template
from numpy import format_float_scientific
from decimal import *
import csv

#pop='MCE3T0'
pop_in_fasta='MCH1'
#Be sure that HGT and reference strain IDs are consistent in bam, mauve, and genbank files
#for pop in ['T1CM','T1MTZ','T2Ab-','T2CLR','T2CM','T2MTZ','T4Ab-','T4CLR','T4CM','T4MTZ','T5Ab-','T5CLR','T5CM','T5MTZ']:
#for pop in ['428F1abT1','428F1abT3','428F1abT5','428F1noT1','428F1noT3','428F1noT5']:
for pop in [
    '428F1',
    '428G1',
    '428G4',
    '428G5',
    '428G6',
    '428H1',
    '428H2',
    '428H4',
    '428H5',
    '428H6',
    'S428E2',
    'S428E3',
    ]:
    #pop_in_fasta=pop
    #t='5'
    ref_strain='NC_011498.1'
    #Use Mauve files where sequence 0 (a) is the reference
    #In Mauve alignments, seq1 is the real reference (P12, for example)
    #mauve_SNPs = '/Users/lwoo0005/Documents/An_H_pylori/Mauve_P12_v_maxHGTs/Mauve_P12_'+pop_in_fasta+'_maxHGT_genome_SNPs'
    #mauve_SNPs='/Users/lwoo0005/Documents/An_H_pylori/Mauve_P12_v_maxHGTs/Mauve_P12_426F4_maxHGT_genome_SNPs'
    #mauve_gaps = '/Users/lwoo0005/Documents/An_H_pylori/Mauve_P12_v_maxHGTs/Mauve_P12_'+pop_in_fasta+'_maxHGT_genome_gaps'
    #mauve_gaps='/Users/lwoo0005/ab_mauve_maxHGTs/Mauve_P12_ShamiFA_MCE3T0_maxHGT_gaps'
    #Open the BAM file
    mauve_SNPs='/Users/lwoo0005/Documents/An_H_pylori/mauve_snps_n_gaps/Mauve_P12_v_maxHGTs/Mauve_P12_428_rdxA_from_MCH1_SNPs'
    #mauve_SNPs = '/Users/lwoo0005/Documents/An_H_pylori/mauve_snps_n_gaps/Mauve_P12_'+pop_in_fasta+'_misc_SNPs'
    #mauve_SNPs='/Users/lwoo0005/Documents/An_H_pylori/Mauve_P12_v_maxHGTs/Mauve_P12_426F4_maxHGT_genome_SNPs'
    mauve_gaps='/Users/lwoo0005/Documents/An_H_pylori/mauve_snps_n_gaps/Mauve_P12_v_maxHGTs/Mauve_P12_428_frxA_rdxA_HGT_genome_gaps'
    #bam_file='/Users/lwoo0005/Documents/An_H_pylori/evolved_strains/'+pop+'_vs_P12_w_HPP12_w_'+pop+'_maxHGT_genome_breseq_p0_bam_info/'+pop+'_vs_P12_w_HPP12_w_'+pop+'_maxHGT_genome_breseq_p0_alignment.bam'
    #bam_file='/Users/lwoo0005/Documents/An_H_pylori/Shami_stuff/bam_data/'+pop+'_L8_breseq_v_P12_w_HPP12_w_'+pop_in_fasta+'_maxHGT_genome_p0_bam_info/reference.bam'
    #fasta_file='/Users/lwoo0005/Documents/An_H_pylori/evolved_strains/'+pop+'_vs_P12_w_HPP12_w_'+pop+'_maxHGT_genome_breseq_p0_bam_info/'+pop+'_vs_P12_w_HPP12_w_'+pop+'_maxHGT_genome_breseq_p0_ref.fasta'
    #fasta_file='/Users/lwoo0005/Documents/An_H_pylori/Shami_stuff/bam_data/'+pop+'_L8_breseq_v_P12_w_HPP12_w_'+pop_in_fasta+'_maxHGT_genome_p0_bam_info/reference.fasta'
    #bam_file= '/Users/lwoo0005/Documents/Laura_stuff/H_py_An/Concatanated_maxHGT_p12_genomes/bam_data/'+pop+'_vs_P12_w_HPP12_w_'+pop+'_maxHGT_genome_breseq_p0_bam_info/'+pop+'_vs_P12_w_HPP12_w_'+pop+'_maxHGT_genome_breseq_p0_alignment.bam'
    #bam_file='/Users/lwoo0005/Documents/An_H_pylori/Clr_Crtl_Fitness_assay/v_maxHGT_bam_data/'+pop+'T'+t+'_v_'+pop+'_maxHGT_genome_w_P12_w_HPP12_bam_info/reference.bam'
    #bam_data = pysam.AlignmentFile(bam_file, "rb")
    #Pathway to the reference file associated with the bam file
    #fasta_file='/Users/lwoo0005/Documents/Laura_stuff/H_py_An/Concatanated_maxHGT_p12_genomes/bam_data/'+pop+'_vs_P12_w_HPP12_w_'+pop+'_maxHGT_genome_breseq_p0_bam_info/'+pop+'_vs_P12_w_HPP12_w_'+pop+'_maxHGT_genome_breseq_p0_ref.fasta'
    #fasta_file='/Users/lwoo0005/Documents/An_H_pylori/Clr_Crtl_Fitness_assay/v_maxHGT_bam_data/'+pop+'T'+t+'_v_'+pop+'_maxHGT_genome_w_P12_w_HPP12_bam_info/reference.fasta'
    bam_file='/Users/lwoo0005/Documents/An_H_pylori/rdxA_frxA_HGT_search/no_Mtz_428_frxA_rdxA_HGT_search_bam_data/'+pop+'_vs_P12_w_HPP12_'+pop_in_fasta+'_maxHGT_genome_breseq_p0_bam_info/'+pop+'_vs_P12_w_HPP12_'+pop_in_fasta+'_maxHGT_genome_breseq_p0_alignment.bam'
    fasta_file='/Users/lwoo0005/Documents/An_H_pylori/rdxA_frxA_HGT_search/no_Mtz_428_frxA_rdxA_HGT_search_bam_data/'+pop+'_vs_P12_w_HPP12_'+pop_in_fasta+'_maxHGT_genome_breseq_p0_bam_info/'+pop+'_vs_P12_w_HPP12_'+pop_in_fasta+'_maxHGT_genome_breseq_p0_ref.fasta'
    
    gb_inp=''
    
    #/Users/lwoo0005/Documents/An_H_pylori/Clr_Crtl_Fitness_assay/All_428G5T5_breseq_v_P12_w_HPP12_p0_HGT_events_on_P12_as_counts.csv
    real_counts_file='/Users/lwoo0005/Documents/An_H_pylori/HGT_binned_for_HGT_only_All_'+pop+'_428_rdxA_HGT_from_MCH1_on_P12_as_counts.csv'
    #Output location for the resultant vcf file
    gd_op_path='/Users/lwoo0005/Documents/An_H_pylori/HGT_binned_for_HGT_only_All_'+pop+'_428_rdxA_HGT_from_MCH1_on_P12.gd'    
    #real_counts_file='/Users/lwoo0005/Documents/An_H_pylori/All_'+pop+'_breseq_v_P12_w_HPP12_p0_HGT_events_on_P12_as_counts.csv'
    #real_counts_file='/Users/lwoo0005/Documents/An_H_pylori/Clr_Crtl_Fitness_assay/All_'+pop+'T'+t+'_breseq_v_P12_w_HPP12_p0_HGT_events_on_P12_as_counts.csv'
    #Output location for the resultant vcf file
    #gd_op_path='/Users/lwoo0005/Documents/Laura_stuff/H_py_An/S426E6_investigation/All_'+pop+'_HGT_events_on_P12.gd'
    #gd_op_path='/Users/lwoo0005/Documents/An_H_pylori/All_'+pop+'_breseq_v_P12_w_HPP12_p0_HGT_events_on_P12.gd'
    #Open gd file for P12+max breseq run as dictionary, where key=position,
    #value[0]=change, value[1]=allele frequency, value[2]=depth (we will use this
    #wherever possible, as I expect breseq to have better cutoffs than me :))
    #Use value[1]=allele_frequency for P12, but value[1]=1-allele_frequency
    #for maxHGT
    #Also, change multy SNPs to single SNPs
    """
    breseq_gd = '/Users/lwoo0005/Documents/Laura_stuff/H_py_An/Concatanated_maxHGT_p12_genomes/output.gd'
    gd_data = [line.strip() for line in open(breseq_gd, 'r') if line.split('\t')[0] not in ['RA','MC','JC','UN']][13:]
    gd_ref_dic={}
    gd_alt_dic={}
    for line in gd_data:
        line_info=line.split('\t')
        if line_info[0]=='SNP':
            gd_chro, gd_pos, gd_HGT_type, gd_HGT_change, gd_freq=str(line_info[3]),int(line_info[4]),str(line_info[0]),str(line_info[5]),float(line_info[6])
            if pop in gd_chro:
                gd_HGT_freq=1-gd_freq
                gd_alt_dic[gd_pos]=[gd_HGT_change, gd_HGT_freq]
            elif ref_strain in gd_chro:
                gd_HGT_freq=gd_freq
                gd_ref_dic[gd_pos]=[gd_HGT_change, gd_HGT_freq]
        elif line_info[0]=='SUB':
            gd_chro, gd_pos, gd_HGT_type, gd_change_length, gd_HGT_change, gd_freq=str(line_info[3]),int(line_info[4]),str(line_info[0]),int(line_info[5]),str(line_info[6]),float(line_info[7])
            if pop in gd_chro:
                gd_HGT_freq=1-gd_freq
                for i in range(len(gd_change_length)):
                    gd_pos=gd_pos+i
                    gd_change=gd_HGT_change[i]
                    gd_alt_dic[gd_pos]=[gd_HGT_change, gd_HGT_freq]
            elif ref_strain in gd_chro:
                gd_HGT_freq=gd_freq
                for i in range(len(gd_change_length)):
                    gd_pos=gd_pos+i
                    gd_change=gd_HGT_change[i]
                    gd_ref_dic[gd_pos]=[gd_HGT_change, gd_HGT_freq]
            elif line_info[0]=='DEL':
            #Same as SNP but change is a number (int)
                gd_chro, gd_pos, gd_HGT_type, gd_HGT_change, gd_freq=str(line_info[3]),int(line_info[4]),str(line_info[0]),int(line_info[5]),float(line_info[6])
                if pop in gd_chro:
                    gd_HGT_freq=1-gd_freq
                    gd_alt_dic[gd_pos]=[gd_HGT_change, gd_HGT_freq]
                elif ref_strain in gd_chro:
                    gd_HGT_freq=gd_freq
                    gd_ref_dic[gd_pos]=[gd_HGT_change, gd_HGT_freq]  
            elif line_info[0]=='INS':                                                                                                                           for i in range(gd_change_length):
            #More complex, expect issues here
                                                                                                                                             
        if pop in gd_chro:
            gd_HGT_freq=1-gd_freq
            gd_alt_dic[gd_pos]=[gd_HGT_change, gd_HGT_freq]
        elif ref_strain in gd_chro:
            gd_HGT_freq=gd_freq
            gd_ref_dic[gd_pos]=[gd_HGT_change, gd_HGT_freq]
    """
    
    #Mauve SNPs and gaps are seperate files, so the best plan of action would be
    # to combine the SNPs and gaps into lists for each comparison and then add the
    # tags.
    
    
    #ab means a with reference to b
    with open(mauve_SNPs, 'r') as SNP_handle:
        SNPs=[line.strip('\n') for line in SNP_handle][1:]
        SNP_handle.close()
    
    SNP_ab=[]
    SNP_ba=[]
    SNP_pos_ab=[]
    SNP_pos_ba=[]
    SNP_lengths = []
    for SNP in SNPs:
        # Mauve shows actual base at location, so swapped to show base changed to
        SNP_ba.append(SNP.split('\t')[0][1])
        SNP_ab.append(SNP.split('\t')[0][0])
        # Positions don't need to be swapped
        SNP_pos_ba.append(int(SNP.split('\t')[3]))
        SNP_pos_ab.append(int(SNP.split('\t')[6]))
        SNP_lengths.append(1)
    with open(mauve_gaps, 'r') as gap_handle:
        gaps=[line.strip('\n') for line in gap_handle][1:]
        gap_handle.close()
    
    all_recs=[]
    rec_ids=[]
    del_ab=[]
    del_ba=[]
    ins_ab=[]
    ins_ba=[]
    ins_len_ab=[]
    ins_len_ba=[]
    del_pos_ab=[]
    del_pos_ba=[]
    ins_pos_ab=[]
    ins_pos_ba=[]
    
    records = [rec for rec in SeqIO.parse(fasta_file, "fasta")]
    for rec in records:
        all_recs.append(rec.id)
        if pop_in_fasta in rec.id:
            for gap in gaps:
                gap_cols=gap.split('\t')
                if gap_cols[0] == 'sequence_0':
                    #col5 always seq0
                    #col8 always seq1
                    del_pos_ab.append(int(gap_cols[8]))
                    in_pos_opposite=int(gap_cols[8])
                    ins_pos_ba.append(int(gap_cols[3]))
                    del_ab.append(int(gap_cols[4]))
                    in_len=int(gap_cols[4])
                    ins_len_ba.append(in_len)
                    #ins_ba.append(in_len)
                    ins_ba.append(str(rec.seq[in_pos_opposite:in_pos_opposite+in_len]))
                if gap_cols[0] == 'sequence_1':
                    del_pos_ba.append(int(gap_cols[5]))
                    ins_pos_ab.append(int(gap_cols[3]))
                    del_ba.append(int(gap_cols[4]))
                    in_len=int(gap_cols[4])
                    ins_len_ab.append(in_len)
                    ins_ab.append(in_len)
    
    #Added this way so that pop is second
    for rec in all_recs:
        if rec in ref_strain:
            rec_ids.append(rec)
    for rec in all_recs:
        if pop_in_fasta in rec:
            rec_ids.append(rec)        
                  
    mut_lengths_ba = [SNP_lengths, del_ba, ins_len_ba]                
    mut_lengths_ab = [SNP_lengths, ins_len_ab, del_ab]
    
    HGT_lists = [[SNP_ba, del_ba, ins_ba],[SNP_ab, ins_ab, del_ab]]
    pos_lists = [[SNP_pos_ba, del_pos_ba, ins_pos_ba],[SNP_pos_ab, ins_pos_ab, del_pos_ab]]
    mut_length_lists = [mut_lengths_ba, mut_lengths_ab]
    type_lists = [['SNP', 'DEL', 'INS'],['SNP', 'INS', 'DEL']]
    
    #write text file with each line showing chromosome and position to check
    #But how do I tag the positions properly if there are just regions??
    #Write the lines ref1, alt1, ref2, alt2, etc.
    #
    
    n=0
    sam_tags=[]
    ref_pos=[]
    change_from_ref=[]
    sam_covs=[]
    for n in [0,1]:
    #for n in [0]:
        print("New strain")
        #For ref (n=0) and then for maxy (n=1)
        sam_cov=[]
        n_rec_id=rec_ids[n]
        n_pos_list=pos_lists[n]
        n_HGT_list = HGT_lists[n]
        n_mut_length_list = mut_length_lists[n]
        n_type_list = type_lists[n]
        for m in range(len(n_pos_list)):
        #for m in [1,2]:
            #Going through mutation types: SNPs, ins, dels
            pos_list=n_pos_list[m]
            HGT_list=n_HGT_list[m]
            mut_length_list=n_mut_length_list[m]
            mut_type=n_type_list[m]
            for HGT, pos, mut_length, idx in zip(HGT_list, pos_list, mut_length_list, range(len(pos_list))):
                if mut_type == "SNP":
                    if n==0:
                        sam_tags.append("SNP")
                        ref_pos.append(pos)
                        change_from_ref.append(HGT)
                    sam_pos=str(n_rec_id)+':'+str(pos)+'-'+str(pos)
                    pysam_pileup_ref=(pysam.mpileup("-r", sam_pos, "-f", fasta_file, bam_file)).split('\t')
                    if len(pysam_pileup_ref) <= 1:
                        ref_cov=0
                        alt_cov=0
                    elif len(pysam_pileup_ref) > 1:
                        if n==1:
                            if pysam_pileup_ref[2]==HGT:
                                new_var=str(HGT_lists[0][0][idx])
                                ref_cov=int(pysam_pileup_ref[4].count(new_var))+int(pysam_pileup_ref[4].count(new_var.lower()))
                                alt_cov=int(pysam_pileup_ref[4].count('.'))+int(pysam_pileup_ref[4].count(','))
                            else:
                                ref_cov=int(pysam_pileup_ref[4].count('.'))+int(pysam_pileup_ref[4].count(','))
                                sam_HGT=str(HGT)
                                alt_cov=int(pysam_pileup_ref[4].count(sam_HGT))+int(pysam_pileup_ref[4].count(sam_HGT.lower()))
                        else:
                            ref_cov=int(pysam_pileup_ref[4].count('.'))+int(pysam_pileup_ref[4].count(','))
                            sam_HGT=str(HGT)
                            alt_cov=int(pysam_pileup_ref[4].count(sam_HGT))+int(pysam_pileup_ref[4].count(sam_HGT.lower()))      
                    sam_cov.append([ref_cov,alt_cov])
                if mut_type in ("INS","DEL"):
                    if n==0:
                        sam_tags.append(mut_type)
                        ref_pos.append(pos)
                        change_from_ref.append(HGT)
                    if mut_type=="INS":
                        sam_pos=str(n_rec_id)+':'+str(pos)+'-'+str(pos)
                            #sam_pos="MCA1:1205268-1205268"
                            #Insertions in sam format are very tricky
                            #Can't seem to bring them up here, so will likely
                            #be ignored if not aligned properly
                        pysam_pileup_ref=(pysam.mpileup("-r", sam_pos, "-f", fasta_file, bam_file)).split('\t')
                        if len(pysam_pileup_ref) <= 1:
                            ref_cov=0
                            alt_cov=0
                        else:
                            ref_cov=(int(pysam_pileup_ref[4].count('.'))+int(pysam_pileup_ref[4].count(',')))
                            alt_cov=0
                        sam_cov.append([ref_cov,alt_cov])
                    elif mut_type=="DEL":
                        sam_pos=str(n_rec_id)+':'+str(pos)+'-'+str(pos)
                        sam_HGT='*'
                        pysam_pileup_ref=(pysam.mpileup("-r", sam_pos, "-f", fasta_file, bam_file)).split('\t')
                        if len(pysam_pileup_ref) <= 1:
                            ref_cov=0
                            alt_cov=0
                        else:
                            ref_cov=(int(pysam_pileup_ref[4].count('.'))+int(pysam_pileup_ref[4].count(',')))
                            alt_cov=int(pysam_pileup_ref[4].count(sam_HGT))
                        sam_cov.append([ref_cov,alt_cov])
                #print(sam_pos)
        sam_covs.append(sam_cov)
    
    refer_cov=[]
    HGT_cov=[]
    total_cov=[]
    sam_cov_a=sam_covs[0]
    sam_cov_b=sam_covs[1]
    
          
    
    #SNP 535
#Use this for normal ref+HGT consideration; use the other for just HGT bin 
# (inlcuding refer only for reference coverage)
    """   
    for refer, HGT in zip(sam_cov_a, sam_cov_b):
        refer_cov.append(int(refer[0])+int(HGT[1]))
        HGT_cov.append(int(refer[1]+int(HGT[0])))
        total_cov.append(int(refer[0])+int(HGT[1])+int(refer[1])+int(HGT[0]))
        if int(refer[0])+int(HGT[1])+int(refer[1])+int(HGT[0]) == 0:
            print("This one:" +str(n))
 """ 
    for refer, HGT in zip(sam_cov_a, sam_cov_b):
        refer_cov.append(int(refer[0])+int(HGT[1]))
        HGT_cov.append(int(HGT[0]))
        total_cov.append(int(refer[0])+int(HGT[1])+int(HGT[0]))
        if int(refer[0])+int(HGT[1])+int(refer[1])+int(HGT[0]) == 0:
            print("This one:" +str(n))          
    n=0
    HGT_freq=[]
    ref_freq=[]
    for HGT,refer,total in zip(HGT_cov, refer_cov, total_cov):
        if total==0:
            print("Here")
            print(n)
            HGT_freq.append(0)
            ref_freq.append(0)
        else:
            HGT_freq.append(Decimal(HGT)/Decimal(total))
            ref_freq.append(Decimal(refer)/Decimal(total))
        n+=1   
            
            
    #HGT_freq=[Decimal(HGT/total) for HGT, total in zip(HGT_cov, total_cov)]
    #ref_freq=[Decimal(refer)/Decimal(total) for refer, total in zip(refer_cov, total_cov)]
      
    num_HGTs=len(total_cov)
    chroms=[ref_strain for i in range(num_HGTs)]    
    ids=["." for chrom in chroms]
    
    #Stuff for gd  
    count_col=[i+1 for i in range(len(total_cov))]
    HGT_freq_col=[]
    for freq in HGT_freq:
        if freq==1.0:
            HGT_freq_col.append("frequency=1")
        elif freq=="No coverage":
            HGT_freq_col.append(freq)
        elif freq==0:
            HGT_freq_col.append("frequency=0")
        else:
            if int(str(freq)[-1])<=8:
                freq=float(freq)+float(0.000000000000001)
            elif int(str(freq)[-1])==9:
                freq=float(freq)-float(0.000000000000001)
            gd_freq_a=str(format_float_scientific(freq, exp_digits=2))[0:9]
            gd_freq_b=str(format_float_scientific(freq, exp_digits=2))[-5:]
            gd_freq=gd_freq_a+gd_freq_b
            HGT_freq_col.append("frequency="+str(gd_freq))
    
    parent_col=ids
    
    
    #Just need to figure out how to report insertions properly
    #For writing a gd file:
    #SNP:
    #SNP	17	485	MCA1	119017	A	frequency=1
    #DEL:
    #DEL	36	513	MCA1	205036	1	frequency=1
    #INS:
    #INS	35	512	MCA1	202955	T	frequency=1	insert_position=1
    
    #Can place "." for third column (Parent IDs)
    #col1=sam_tags, col2 is just a count, col4=chroms, col5=ref_pos, col6=HGT_freq
    #but with frequency= and col7 should be '' for SNPs and DELs and "insert_position=1" for INS
    
    #Change this if possible
    
    gd_data=[]
    gd_header=str("#=GENOME_DIFF	1.0\n#=CREATED	14:32:21 10 Jun 2020\n#=PROGRAM	breseq/HGTquantifier\n#=COMMAND	\n#=READSEQ\n#=CONVERTED-BASES	\n#=CONVERTED-READS	\n#=INPUT-BASES	\n#=INPUT-READS	\n#=MAPPED-BASES\n#=MAPPED-READS")
    gd_data.append(gd_header)
    for mut, mut_id, parent_id, chrom, pos, var, freq in zip(sam_tags, count_col, ids, chroms, ref_pos, change_from_ref, HGT_freq_col):
        if mut=="DEL":
            pos+=1
        line="\n"+("\t").join([str(mut), str(mut_id), str(parent_id), str(chrom), str(pos), str(var), str(freq)])
        if mut=="INS":
            line=line+"\tinsert_position=1"
        if mut=="SNP":
            if var in ['Y','R','W','S','K','M','K','H','B','D','V','X']:
                var='N'
        gd_data.append(line)

    with open(gd_op_path, 'w') as gd:
        for line in gd_data:
            gd.write(line)
        gd.close()

    """
    vcf_data=[]
    vcf_header = "#Writing the VCF file\n##fileformat=VCFv4.1\n##fileDate\n##source=Laura\'s cool program\n##contig=<ID=NC_011498.1,length=1673813>\n##contig=<ID=NC_011499.1,length=10225>\n##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n##INFO=<ID=AD,Number=1,Type=Float,Description=\"Allele Depth (avg read count)\">\n##INFO=<ID=DP,Number=1,Type=Float,Description=\"Total Depth (avg read count)\">\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n"
    with open(vcf_op, 'w') as vcf_out:
        vcf_out.write(vcf_header)
        vcf_out.close()
    """    
    with open(real_counts_file, 'w') as real_counts:
        writer=csv.writer(real_counts)
        writer.writerow(['Position','Mutation','Reference coverage','Alternate coverage'])
        for mut, pos, var, r, a in zip(sam_tags, ref_pos, change_from_ref, refer_cov,HGT_cov):
            writer.writerow([pos, mut, var, r,a])