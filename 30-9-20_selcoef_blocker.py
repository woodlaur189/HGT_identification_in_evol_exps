#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 06:24:33 2020

@author: lwoo0005
"""


#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 17:01:51 2020

@author: lwoo0005
"""
import numpy as np
import pandas as pd

"""
pop='426F4_all_minD0_for_image'
inp_path='/Users/lwoo0005/Documents/An_H_pylori/Clr_Crtl_Fitness_assay/All_selection_coefs_426F4.xlsx'
sheet='426F4_all_minD0_for_image'
op_path='/Users/lwoo0005/Documents/An_H_pylori/Clr_Crtl_Fitness_assay/'

"""
inp_path='/Users/lwoo0005/Documents/An_H_pylori/Shami_stuff/partMCE3_allvars_selection_coefs.xlsx'
for pop in ['partAb-']:
    #pop='MCE3_allvars_Clr_rm1tp_sig'
    sheet='MCE3_'+pop+'_allvars_rm1tp_sig'
    op_path='/Users/lwoo0005/Documents/An_H_pylori/Shami_stuff/'
    
    
    data=pd.read_excel(open(inp_path,'rb'), sheet_name=sheet)
    #data=data[:600]
    data=data.sort_values('Position')
    
    #Should also have DNp5, HGT, ancestral, other as metadata--add to selcoef sheet
    pos_list=data['Position'].tolist()
    selcoef_list=data['Coefficient'].tolist()
    lowers_list=data['Lower_CI_bound'].tolist()
    uppers_list=data['Upper_CI_bound'].tolist()
    gene_list=data['Gene'].tolist()
    description_list=data['Description'].tolist()
    mutation_list=data['S/NS'].tolist()
    annotation_list=data['Mutation_annotation'].tolist()
    type_list=data['Type'].tolist()
    pval_list=data['p-value'].tolist() 
    
    """
    pos_list=[10,12,14,23,25]
    selcoef_list=[-0.2,-0.22,-0.18,0.1,0.14]
    lowers_list=[-0.23,-0.25,-0.2,0,0.02]
    uppers_list=[-0.17,-0.19,-0.16,0.2,0.26]
    """
    
    i=0
    block_dic={}
    block_number=1
    old_lower=0
    old_upper=0
    block_positions=[]
    block_selcoefs=[]
    block_genes=[]
    block_descriptions=[]
    block_mutations=[]
    block_annotations=[]
    block_types=[]
    block_pvals=[]
    
    for p, s, lower, upper, g, d, m, a, t, pv in zip(pos_list, selcoef_list, lowers_list, 
                                             uppers_list, gene_list, description_list,
                                             mutation_list, annotation_list, type_list,
                                             pval_list):
       #print s
       #print old_lower
       #print old_upper
       i+=1
       if i==1:
           block_positions.append(p)
           block_selcoefs.append(s)
           block_genes.append(g)
           block_descriptions.append(d)
           block_mutations.append(m)
           block_annotations.append(a)
           block_types.append(t)
           block_pvals.append(pv)
           block_dic[block_number]=[block_positions,block_selcoefs,block_genes,
                    block_descriptions,block_mutations, block_annotations,
                    block_types, block_pvals]
           old_lower=lower
           old_upper=upper
       else:
           if s >= old_lower and s <= old_upper:
               #print "Within bounds"
               block_positions.append(p)
               block_selcoefs.append(s)
               block_genes.append(g)
               block_descriptions.append(d)
               block_mutations.append(m)
               block_annotations.append(a)
               block_types.append(t)
               block_pvals.append(pv)
               block_dic[block_number]=[block_positions,block_selcoefs,block_genes,
                    block_descriptions,block_mutations, block_annotations,
                    block_types, block_pvals]
               old_lower=lower
               old_upper=upper
           else:
               #print "outside bounds"
               block_positions=[]
               block_selcoefs=[]
               block_genes=[]
               block_descriptions=[]
               block_mutations=[]
               block_annotations=[]
               block_types=[]
               block_pvals=[]
               block_number+=1
               block_positions.append(p)
               block_selcoefs.append(s)
               block_genes.append(g)
               block_descriptions.append(d)
               block_mutations.append(m)
               block_annotations.append(a)
               block_types.append(t)
               block_pvals.append(pv)
               block_dic[block_number]=[block_positions,block_selcoefs,block_genes,
                    block_descriptions,block_mutations, block_annotations,
                    block_types, block_pvals]
               old_lower=lower
               old_upper=upper
    
    block_positions={}
    for k,v in block_dic.items():
       pos=v[0]
       selcoefs=v[1]
       genes=v[2]
       descriptions=v[3]
       mutations=v[4]
       annotations=v[5]
       types=v[6]
       vals=v[7]
       num_vars=len(types)
       num_HGT=types.count('HGT')
       freq_HGT=float(num_HGT)/float(len(types))
       num_syn=mutations.count("S")
       freq_syn=float(num_syn)/float(len(mutations))
       if 'Coding' in mutations:
           coding_flag='Yes'
           coding_indices = [a for a, x in enumerate(mutations) if x == "Coding"]
           coding_selcoef=', '.join([str(selcoefs[b]) for b in coding_indices])
           coding_pval=', '.join([str(vals[b]) for b in coding_indices])
       else:
           coding_flag='No'
           coding_selcoef=', '.join(['NA'])
           coding_pval=', '.join(['NA'])
       if 'Nonsense' in mutations:
           nonsense_flag='Yes'
           nonsense_indices = [a for a, x in enumerate(mutations) if x == "Nonsense"]
           nonsense_selcoef=', '.join([str(selcoefs[b]) for b in nonsense_indices])
           nonsense_pval=', '.join([str(vals[b]) for b in nonsense_indices])
       else:
           nonsense_flag='No'
           nonsense_selcoef=', '.join(['NA'])
           nonsense_pval=', '.join(['NA'])
       if 'Nonstop' in mutations:
           nonstop_flag='Yes'
           nonstop_indices = [a for a, x in enumerate(mutations) if x == "Nonstop"]
           nonstop_selcoef=', '.join([str(selcoefs[b]) for b in nonstop_indices])
           nonstop_pval=', '.join([str(vals[b]) for b in nonstop_indices])
       else:
           nonstop_flag='No'
           nonstop_selcoef=', '.join(['NA'])
           nonstop_pval=', '.join(['NA'])
       if 'NS' in mutations:
           missense_flag='Yes'
           missense_indices = [a for a, x in enumerate(mutations) if x == "NS"]
           missense_selcoef=', '.join([str(selcoefs[b]) for b in missense_indices])
           missense_pval=', '.join([str(vals[b]) for b in missense_indices])
       else:
           missense_flag='No'
           missense_selcoef=', '.join(['NA'])
           missense_pval=', '.join(['NA'])
       min_pos=min(pos)
       max_pos=max(pos)
       min_gene=genes[pos.index(min_pos)].encode('ascii', 'ignore').decode('ascii')
       max_gene=genes[pos.index(max_pos)].encode('ascii', 'ignore').decode('ascii')
       avg_sc=np.average(selcoefs)
       min_sc=min(selcoefs)
       max_sc=max(selcoefs)
       driver_index=np.argmax(np.abs(selcoefs))
       extreme_sc=selcoefs[driver_index]
       driver=pos[driver_index]
       driver_gene=genes[driver_index].encode('ascii', 'ignore').decode('ascii')
       driver_description=descriptions[driver_index].encode('ascii', 'ignore').decode('ascii')
       driver_mutation=mutations[driver_index].encode('ascii', 'ignore').decode('ascii')
       driver_annotation=annotations[driver_index].encode('ascii', 'ignore').decode('ascii')
       driver_type=types[driver_index].encode('ascii', 'ignore').decode('ascii')
       driver_pval=vals[driver_index]
       block_positions[k]=[str(min_pos)+'-'+str(max_pos),
                           str(abs(max_pos-min_pos)+1),
                           min_gene, max_gene,
                           str(num_vars), str(freq_HGT),
                           str(freq_syn),
                           str(avg_sc), str(min_sc), str(max_sc),
                           str(driver), str(extreme_sc), 
                           driver_gene, driver_description,
                           driver_mutation, driver_annotation,
                           driver_type, driver_pval,
                           coding_flag, coding_selcoef, coding_pval,
                           nonsense_flag, nonsense_selcoef, nonsense_pval,
                           nonstop_flag, nonstop_selcoef, nonstop_pval,
                           missense_flag, missense_selcoef, missense_pval
                           ]
    pandas_blocks=pd.DataFrame(block_positions).T
    csv_columns=['Position_range','Block_length',
                 'Start_gene', 'End_gene',
                 'Number_of_vars', 'Frequency_HGT',
                 'Frequency_synonymous',
              'Average_SC', 'Min_SC', 'Max_SC',
              'Putative_driver_pos','Putative_driver_SC', 
              'Putative_driver_gene','Putative_driver_description',
              'Putative_driver_mutation','Putative_driver_annotation',
              'Putative_driver_type', 'Putative_driver_p_value',
              'Coding_in_block?', 'Coding_selcef(s)', 'Coding_p_value(s)',
              'Nonsense_in_block?', 'Nonsense_selcoef(s)', 'Nonsense_p_value(s)',
              'Nonstop_in_block?', 'Nonstop_selcoef(s)', 'Nonstop_p_value(s)',
              'Missense_in_block?', 'Missense_selcoef(s)', 'Missense_p_values']
    pandas_blocks.columns=csv_columns
    pandas_blocks.to_csv(op_path+str('/')+sheet+"_selcoef_blocker_more_info.csv",encoding='utf-8')
"""
with open(op_path+str('/')+pop+"_selcoef_blocker_var_type_info.csv",'w') as op:
   writer=csv.writer(op)
   writer.writerow(csv_columns)
   for k,v in block_positions.items():
       writer.writerow([k,v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10]])
""" 