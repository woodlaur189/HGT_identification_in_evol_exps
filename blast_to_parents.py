#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 17:13:01 2020

@author: lwoo0005
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 18:51:44 2018

@author: lwoo0005
"""

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from glob import glob
import os, errno

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
    db = db, outfmt = "5", out= out, evalue=1e-100)
    stdout, stderr = blastn_cline()
    result_handle = open(out)
    blast_records = NCBIXML.parse(result_handle)
    for record in blast_records:
        if len(record.alignments) >= 1:
            for alignment in record.alignments:
                for hsp in alignment.hsps:
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
%s...
%s...
%s...
\n
"""                  
                    query_name = record.query
                    quer_all = list(hsp.query)
                    quer_gaps = quer_all.count("-")
                    quer_start = int(hsp.query_start)
                    quer_end = quer_start+int(len(hsp.query))-quer_gaps
                    sequence = str(alignment.title)
                    sub_all = list(hsp.sbjct)
                    sub_gaps = sub_all.count("-")
                    sub_start = int(hsp.sbjct_start)
                    sub_end = sub_start+int(len(hsp.sbjct))-sub_gaps
                    length = str(hsp.align_length)
                    exp = str(hsp.expect)
                    bit = str(hsp.bits)
                    sco = str(hsp.score)
                    ids = str(hsp.identities)
                    pos = str(hsp.positives)
                    gaps = str(hsp.gaps)
                    quer = hsp.query
                    match = hsp.match
                    subj = hsp.sbjct
                    results_dir = output_path+"/Scaffold_blaster/Blast_Results/"
                    if not os.path.exists(results_dir):
                        try:
                            os.makedirs(results_dir)
                        except OSError as exc:
                            if exc.errno != errno.EEXIST:
                                raise
                    blast_results=open(results_dir+label+".txt","a")
                    blast_results.write(str(template %(query_name, quer_start, quer_end, sequence, sub_start, sub_end, length, exp, bit, sco, ids, length, pos, length, gaps, quer, match, subj)))
                    blast_results.close()            
    return

qfile = "/Users/lwoo0005/Documents/Laura_stuff/H_py_An/local_roary/roary_90_all_evolved_parents/_1592894899/evolve_uniques/evolved_unique_genes.fasta"
db = "/Users/lwoo0005/Documents/Laura_stuff/H_py_An/Unmatched_donor_scaffs/Filtering/P12_HPP12_CH426hybrid_CH428hybrid.fasta"
output_path ="/Users/lwoo0005/Documents/Laura_stuff/H_py_An/local_roary/roary_90_all_evolved_parents/_1592894899/Scaffold_blast_evolved_uniques"
label = qfile.split("/")[-1].split(".fasta")[0]
scaf_blaster(qfile, output_path, label, db)