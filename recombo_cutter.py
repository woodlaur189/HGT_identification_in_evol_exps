#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 14:50:10 2020

@author: lwoo0005
"""

#Cutting out the donor sequences based on the coverage regions.
from Bio import SeqIO

op_path='/Users/lwoo0005/Documents/Laura_stuff/H_py_An/'
fasta_file = '/Users/lwoo0005/Documents/Laura_stuff/H_py_An/426A3_max_HGT_consensus_genome.fasta'

#Should be extracted from file name in later updates
pop='426A3'
seq_record = SeqIO.read(fasta_file, "fasta")

#read each line of the MC gd file into a list
block_file = '/Users/lwoo0005/Documents/Laura_stuff/H_py_An/426A3_recombo_blocks.csv'
block_list=[line for line in open(block_file, 'r')][1:]

main=[]
i=1
for line in block_list:
    if 'End' in line.split(',')[2]:
        recombo_dic={i:[line.split(',')[0],line.split(',')[1],len(seq_record.seq)]}
    else:
        recombo_dic={i:[line.split(',')[0],line.split(',')[1],line.split(',')[2]]}
    main.append(recombo_dic)
    i+=1

out_data=[]
cut_num=0
for dic in main:
    for key, value in dic.items():
        block=value[0]
        start=int(value[1])
        end=int(value[2])
        cut= str(seq_record.seq[start-1:end])
        cut_num+=1
        out_data.append('>%s_block_%s\n' % (seq_record.id, block))
        out_data.append('%s\n' % (cut))           

with open(str(op_path)+'/'+str(pop)+'_blocks.fasta', 'w') as op_file:
    for line in out_data:
        op_file.write(str(line))
    op_file.close()
