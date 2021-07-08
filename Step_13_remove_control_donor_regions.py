#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 17:39:59 2020

@author: lwoo0005
"""

from Bio import SeqIO
from glob import glob

#open all control results and collect the names of all nodes to make a uniques
# -only list

control_folder='/Users/lwoo0005/Documents/Laura_stuff/H_py_An/dd_donor_coverage_regions_unfiltered/Controls'
c_files = glob(control_folder+'/*.fasta')

HGT_folder='/Users/lwoo0005/Documents/Laura_stuff/H_py_An/dd_donor_coverage_regions_unfiltered/HGT'
h_files = glob(HGT_folder+'/*.fasta')

results_dir='/Users/lwoo0005/Documents/Laura_stuff/H_py_An/PHun2628p300_donor_coverage_regions/HGT_pops/blast_filtering_against_controls'

c_seq_list = []
for c_file in c_files:
    c_ids = [record.id.split('_coverage_region')[0] for record in SeqIO.parse(c_file, 'fasta')]
    for c_id in c_ids:
        if c_id not in c_seq_list:
            c_seq_list.append(c_id)
    
#remove any nodes from the HGT pop donor coverage regions that share the
# record.id in the list
            
for h_file in h_files:
    filtered_h_recs=[]
    new_h_file = results_dir+'/'+str(((str(h_file).split('/')[-1])).split('.fasta')[0])+'_filtered_against_controls.fasta'
    h_recs = [record for record in SeqIO.parse(h_file, 'fasta')]
    for h_rec in h_recs:
        if h_rec.id.split('_coverage_region')[0] not in c_seq_list:
            filtered_h_recs.append(h_rec)
    SeqIO.write(filtered_h_recs, new_h_file, 'fasta')
