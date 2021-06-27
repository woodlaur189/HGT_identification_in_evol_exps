#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 11:07:12 2020

@author: lwoo0005
"""
#Note: currently written in Python2
#Install the relevant dependencies if running Blast in Python

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from glob import glob
import os, errno

#For each donor scaffold, blast the scaffold against the recipient and any other donor(s) and filter for length to retrieve "pruned" donor scaffolds
def scaf_blaster(qfile, output_path, label, db):
    op = str(output_path)+"/Scaffold_blaster/Blast_XMLs/"
    print op
    out= op+str(label)+".xml"
    print out
    if not os.path.exists(op):
        try:
            os.makedirs(op)
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise
    blastn_cline = NcbiblastnCommandline(query = qfile,
    db = db, outfmt = "5", out= out, evalue=1e-100,
    perc_identity=100)
    stdout, stderr = blastn_cline()
    result_handle = open(out)
    blast_records = NCBIXML.parse(result_handle)
    for record in blast_records:
        if len(record.alignments) >= 1:
            for alignment in record.alignments:
                i=0
                sequence = str(alignment.title)
                for hsp in alignment.hsps:
                    query_length = len(hsp.query)
                    if "NC_01149" not in sequence and "CH426" not in sequence and query_length>=300:
                        i+=1
                        template = """
Query name: %s
Start position of query: %s
End position of query: %s
Subject name: %s
Start position of subject: %s
End position of subject: %s
Alignment length: %s
E-value: %s
Bit score: %s
Score: %s
Identities: %s of %s
Positives: %s of %s
Gaps: %s
Query:      %s %s... %s
Alignment:  %s...
Subject:    %s %s... %s
\n
"""                  
                        fasta_template = """
>%s_cut_to_length_%s
%s
"""             
                        query_name = record.query
                        quer_all = list(hsp.query)
                        quer_gaps = quer_all.count("-")
                        quer_start = int(hsp.query_start)
                        quer_end = quer_start+int(len(hsp.query))-quer_gaps
                        
                        sub_all = list(hsp.sbjct)
                        sub_gaps = sub_all.count("-")
                        sub_start = int(hsp.sbjct_start)
                        sub_end = sub_start+int(len(hsp.sbjct))-sub_gaps
                        length = str(alignment.length)
                        exp = str(hsp.expect)
                        bit = str(hsp.bits)
                        sco = str(hsp.score)
                        ids = str(hsp.identities)
                        pos = str(hsp.positives)
                        gaps = str(hsp.gaps)
                        quer = hsp.query[0:100]
                        query_full = hsp.query
                        match = hsp.match[0:100]
                        subj = hsp.sbjct[0:100]
                        results_dir = output_path+"/Scaffold_blaster/Blast_results/"
                        if not os.path.exists(results_dir):
                            try:
                                os.makedirs(results_dir)
                            except OSError as exc:
                                if exc.errno != errno.EEXIST:
                                    raise
                        blast_results=open(results_dir+"blast_summary_"+str(label)+".txt","a")
                        blast_results.write(str(template %(query_name, quer_start, quer_end, sequence, sub_start, sub_end, length, exp, bit, sco, ids, length, pos, length, gaps, quer_start, quer, quer_end, match, sub_start, subj, sub_end)))
                        blast_results.close()
                        fasta_file=open(results_dir+label+".fasta","a")
                        fasta_file.write(str(fasta_template %(query_name, query_length, query_full)))
                        fasta_file.close()
    return

qfile = "/Users/lwoo0005/Documents/Laura_stuff/H_py_An/Unmatched_donor_scaffs/CH428/Trimmed_DMJM28_L7_v_P12_w_plasmid_breseq_unmatched_spades_scaffolds_300+.fasta"
db = "/Users/lwoo0005/Documents/Laura_stuff/H_py_An/Unmatched_donor_scaffs/Filtering/P12_HPP12_CH426hybrid_CH428hybrid.fasta"
output_path ="/Users/lwoo0005/Documents/Laura_stuff/H_py_An/Unmatched_donor_scaffs/Filtering/un28p300_v_full_donors_P12_HPP12"
label = qfile.split("/")[-1].split(".fasta")[0]
scaf_blaster(qfile, output_path, label, db)
