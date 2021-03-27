#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 01:23:19 2020

@author: lwoo0005
"""

#Remove non-alignments from evolved populations
#If the coverage region doesn't align with any of the three original strains,
#it is tossed (but collected in a folder)

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from glob import glob
import os, errno

def scaf_blaster(qfile, output_path, label, db):
    checked_coverage=[]
    op = str(output_path)+"/Scaffold_blaster/Blast_XMLs/"
    out= op+str(label)+".xml"
    if not os.path.exists(op):
        try:
            os.makedirs(op)
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise
    blastn_cline = NcbiblastnCommandline(query = qfile,
    db = db, outfmt = "5", out= out, evalue=1e-100)
    stdout, stderr = blastn_cline()
    result_handle = open(out)
    blast_records = NCBIXML.parse(result_handle)
    for record in blast_records:
        if len(record.alignments) >= 1:
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                        results_dir = output_path+"/Scaffold_blaster/Blast_results/"
                        if not os.path.exists(results_dir):
                            try:
                                os.makedirs(results_dir)
                            except OSError as exc:
                                if exc.errno != errno.EEXIST:
                                    raise
                        checked_coverage.append(record.query)
                        #Just need to collect the name of the nodes that have blast matches, then filter through the fasta sequences using the list
    return checked_coverage

q_folder='/Users/lwoo0005/Documents/Laura_stuff/H_py_An/dd_donor_coverage_control_filtered/'
db = "/Users/lwoo0005/Documents/Laura_stuff/H_py_An/Unmatched_donor_scaffs/Filtering/P12_HPP12_CH426hybrid_CH428hybrid.fasta"
output_path ="/Users/lwoo0005/Documents/Laura_stuff/H_py_An/dd_donor_coverage_all_filtered"

qs = glob(q_folder+'/*.fasta')

all_checked_coverage=[]
for q in qs:
    qfile=q
    label=qfile.split("/")[-1].split(".fasta")[0]
    checked_cov=scaf_blaster(qfile, output_path, label, db)
    all_checked_coverage+=checked_cov

duds=[]
for q in qs:
    qfile=q
    label=qfile.split("/")[-1].split(".fasta")[0]
    new_cov_file=output_path+'/'+label+'_filtered_through_P12_and_donors.fasta'
    binned_cov_file=output_path+'/'+label+'_binned_through_P12_and_donors.fasta'
    seq_records = [record for record in SeqIO.parse(q, 'fasta')]
    filtered_seq_records=[]
    binned_seq_records=[]
    i=0
    for record in seq_records:
        if record.id in all_checked_coverage:
            i+=1
            filtered_seq_records.append(record)
        else:
            print "In "+label+":"
            print "\nThis one didn't make it:"
            print record.id
            binned_seq_records.append(record)
    SeqIO.write(filtered_seq_records, new_cov_file, 'fasta')
    SeqIO.write(binned_seq_records, binned_cov_file, 'fasta')       
    if i==0:
        duds.append(label)
print "There were "+str(len(duds))+" duds. Sad."
    
