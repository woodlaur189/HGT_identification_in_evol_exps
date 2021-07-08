#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 16:22:28 2020

@author: lwoo0005
"""

from Bio import SeqIO
import numpy as np
import glob as glob

#Requires as input:
#fasta sequences for the concatenated genome (recipient + donor nodes) with mutations applied (see Deatherage and Berrick 2014, gdtools).
#gd files should be in the same place (see below)
singles_folder='/Users/lwoo0005/Documents/Laura_stuff/H_py_An/PHun2628p300_breseq_gd_files/*MUTSON.fasta'


#output folder
op_path='/Users/lwoo0005/Documents/Laura_stuff/H_py_An/'

#Parse the gd_file for missing coverage start and end to get donor node regions with coverage (the spaces between)
def seq_cutter(op_path, fasta_file, pop, donors, gd_file):
    seq_records = [seq_record for seq_record in SeqIO.parse(fasta_file, "fasta")]
    gd_list_a=[line for line in open(gd_file, 'r') if 'MC' == line.split('\t')[0]]
    gd_list=[]
    for d in donors:
        name_file = glob.glob('/Users/lwoo0005/Documents/Laura_stuff/H_py_An/Unmatched_donor_scaffs/Filtering/un2[68]p300_v_full_donors_P12_HPP12/Scaffold_blaster/Blast_results/filtered_Trimmed_' + str(d) +'*_v_P12*_unmatched_spades_scaffolds_300+.fasta')[0]
        names = [line.strip('\n') for line in open(name_file, 'r') if d in line]
        for record in seq_records:  
            for name in names:
                if d in name:
                    if record.id in name:
                        record.id=name.strip('>')
    gd_list_a=[line.strip('/n') for line in open(gd_file, 'r') if 'MC' == line.split('\t')[0]]
    gd_list=[]
    for d in donors:
        name_file = glob.glob('/Users/lwoo0005/Documents/Laura_stuff/H_py_An/Unmatched_donor_scaffs/Filtering/un2[68]p300_v_full_donors_P12_HPP12/Scaffold_blaster/Blast_results/filtered_Trimmed_' + str(d) +'*_v_P12*_unmatched_spades_scaffolds_300+.fasta')[0]
        names = [line.strip('\n') for line in open(name_file, 'r') if d in line]
        for gd in gd_list_a:
            for name in names:
                if d in name:
                    if gd.split('\t')[3] in name:
                        gd_elements=gd.split('\t')
                        gd_elements[3]=name.strip('>')
                        gd_list.append('\t'.join(gd_elements))
    uni_ids = list(np.unique([line.split('\t')[3] for line in gd_list]))
    main=[]
    i=0
    for iden in uni_ids:
        uni_id_dic={}
        for line in gd_list:
            if iden == line.split('\t')[3]:
                i+=1
                uni_id_dic.update({i:[line.split('\t')[3],line.split('\t')[4],line.split('\t')[5]]})
        main.append(uni_id_dic)
    out_data=[]
    cut_lens=[]
    for dic in main[1:]:
        node_length=0
        cut_num=0
        match_seq=''
        match_name=''
        false_start=1
        for key, value in sorted(dic.items()):
            iden=value[0]
            start=int(value[1])
            end=int(value[2])
            for record in seq_records:
                if iden in record.id:
                    node_length = len(record.seq)
                    cut= record.seq[false_start:start-1]
                    if len(cut) >= 1:
                        match_seq = record.seq
                        match_name = record.id
                        cut_num+=1
                        out_data.append('>%s_coverage_region_%i\n' % (match_name, cut_num))
                        out_data.append(str(cut)+'\n')
                        cut_lens.append(len(cut))
            false_start=end             
        final_cut=match_seq[false_start-1:node_length-1]
        if len(final_cut) >= 1:
            cut_num+=1
            out_data.append('>%s_coverage_region_%i\n' % (match_name, cut_num))
            out_data.append(str(final_cut)+'\n')
            cut_lens.append(len(final_cut))
    with open(str(op_path)+'/'+str(pop)+'_donor_coverage_regions.fasta', 'w') as op_file:
        for line in out_data:
            op_file.write(str(line))
        op_file.close()
    print "All done!"
    if len(cut_lens) >= 1:
        print "And the longest node is %i bp long--wow!" % (max(cut_lens)) 

for g in glob.iglob(singles_folder):
    fasta_file=g
    pop=g.split('/')[-1].split('_')[0]
    gd_file=glob.glob('/Users/lwoo0005/Documents/Laura_stuff/H_py_An/PHun2628p300_breseq_gd_files/'+str(pop)+'*op.gd')[0]
    if '426' in pop:
        donors=['DMJM26']
    elif '428' in pop:
        donors=['DMJM28']
    elif 'MCA' in pop:
        donors=['DMJM26']
    elif 'MCB' in pop:
        donors=['DMJM26']
    elif 'MCE4' in pop:
        donors=['DMJM26']
    elif 'MCE5' in pop:
        donors=['DMJM26']
    elif 'MCE6' in pop:
        donors=['DMJM26']
    elif 'MCF4' in pop:
        donors=['DMJM26']
    elif 'MCF5' in pop:
        donors=['DMJM26']
    elif 'MCF6' in pop:
        donors=['DMJM26']
    elif 'MCC2' in pop:
        donors=['DMJM26']
    elif 'MCD3' in pop:
        donors=['DMJM26']
    elif '428' in pop:
        donors=['DMJM28']
    elif 'MCF1' in pop:
        donors=['DMJM28']
    elif 'MCE2' in pop:
        donors=['DMJM28']
    elif 'MCE3' in pop:
        donors=['DMJM28']
    elif 'MCG' in pop:
        donors=['DMJM28']
    elif 'MCH' in pop:
        donors=['DMJM28']
    #Controls
    elif '12' in pop:
        donors=['DMJM26','DMJM28']
    elif 'MCC4' in pop:
        donors=['DMJM26','DMJM28']
    elif 'MCD4' in pop:
        donors=['DMJM26','DMJM28']    
    else:
        break
    seq_cutter(op_path, fasta_file, pop, donors, gd_file)
