#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 16:25:35 2020

@author: lwoo0005
"""

"""
Credit to:

@author: kerry h.

for continuous colour cmaps for use with matplotlib (Lines 22-62).
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def hex_to_rgb(value):
    '''
    Converts hex to rgb colours
    value: string of 6 characters representing a hex colour.
    Returns: list length 3 of RGB values'''
    value = value.strip("#") # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

def rgb_to_dec(value):
    '''
    Converts rgb to decimal colours (i.e. divides each value by 256)
    value: list (length 3) of RGB values
    Returns: list (length 3) of decimal values'''
    return [v/256 for v in value]

def get_continuous_cmap(hex_list, float_list=None):
    ''' creates and returns a color map that can be used in heat map figures.
        If float_list is not provided, colour map graduates linearly between each color in hex_list.
        If float_list is provided, each color in hex_list is mapped to the respective location in float_list. 
        
        Parameters
        ----------
        hex_list: list of hex code strings
        float_list: list of floats between 0 and 1, same length as hex_list. Must start with 0 and end with 1.
        
        Returns
        ----------
        colour map'''
    rgb_list = [rgb_to_dec(hex_to_rgb(i)) for i in hex_list]
    if float_list:
        pass
    else:
        float_list = list(np.linspace(0,1,len(rgb_list)))
        
    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] for i in range(len(float_list))]
        cdict[col] = col_list
    cmp = cm.colors.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=256)
    return cmp

#Colour blind palette
okabe_tio = ["#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]

#Colour maps to colour by different driver mutation types
misc_hex_list=['#F6FFFF','#D3FFFD', '#9ED4D1', '#7BB7B4', '#579491', '#386E6C']
misc_cmap=get_continuous_cmap(misc_hex_list)

DNp5_hex_list=['#FDF2E9', '#FAE5D3', '#F5CBA7', '#F0B27A', '#EB984E',
               '#E67E22', '#CA6F1E', '#AF601A', '#935116', '#784212']
DNp5_cmap=get_continuous_cmap(DNp5_hex_list)

anc_hex_list=['#f9fbe7', '#f0f4c3', '#e6ee9c', '#dce775', '#d4e157', '#cddc39', '#c0ca33',
              '#afb42b', '#9e9d24', '#827717']
anc_cmap=get_continuous_cmap(anc_hex_list)

HGT_hex_list=['#FBF6FF','#E7D3FF', '#B69ED4', '#967BB7', '#725794', '#51386E']

HGT_cmap=get_continuous_cmap(HGT_hex_list)

#polar_hex_list=["#CC79A7","#0072B2"]
polar_hex_list=["#5697AF","#6EB09F","#A6C099","#E5CC50","#DEA73F","#DA7D2E","#E1392E"]
polar_cmap=get_continuous_cmap(polar_hex_list)

import pandas as pd
import matplotlib.gridspec as gridspec

npop='MCE3_partAb-_sig_rm1tp'
inp_path='/Users/lwoo0005/Documents/An_H_pylori/Shami_stuff/partMCE3_allvars_selection_coefs.xlsx'
sheet='MCE3_partAb-_allvars_rm1tp_bloc'
op_path='/Users/lwoo0005/Documents/An_H_pylori/Shami_stuff/'

data=pd.read_excel(open(inp_path,'rb'), sheet_name=sheet)
data=data.sort_values('Block')

# Only including blocks greater than a codon length

data = data.drop(data[data.Block_length < 4].index)

#Length of H. pylori P12 genome
genome_length=1673813

HGT_freq_list=data["Frequency_HGT"].tolist()

a = data['Average_SC'].to_numpy()

#Input minimum and maximum values that correlate to spread of selection coefficients
scale_min=-0.3
scale_max=0.3
sc_list_scaled=np.interp(a, (scale_min, scale_max), (0, 1))

type_list=data['Putative_driver_type'].tolist()
positions=[int(a.split('-')[0]) for a in data['Position_range'].tolist()]
widths=data['Block_length'].tolist()
centers=[p+(0.5*w) for p,w in zip(positions, widths)]

fig = plt.figure(figsize=(20,4))
gs = gridspec.GridSpec(ncols=1, nrows=5, figure=fig)

#ax = fig.add_subplot(212, aspect='equal')
ax = fig.add_subplot(gs[0:,:])
#ax0 = fig.add_subplot(gs[0:2,:], sharex=ax)
ax.set_xlim([0, genome_length])
ax.set_ylim([0, 100000])
base_rec = plt.Rectangle((0,0), width=genome_length, height=100000, color='white')
ax.add_artist(base_rec)
for p,w,t,s,h in zip(positions,widths,type_list, sc_list_scaled, HGT_freq_list):
    x=p
    y = 0
    """
    #For colouring by driver type option.
    if t=='HGT':
        c=HGT_cmap(0.5)
        #c=my_cmap(s)
    elif t== 'DNp5':
        c=DNp5_cmap(0.5)
    elif t=='Ancestral':
        c=anc_cmap(0.5)
    elif '-NS' in t:
        c=cm.Greys(0.5)
    else:
        #c=cm.Purples(s)
        c=misc_cmap(0.5)
    """
    #For colouring by average selection coefficient of block
    c=cm.coolwarm(s)
    new_rec=plt.Rectangle((x, y),width=w,height=h*100000,color=c)
    ax.add_artist(new_rec)
    #ax0.bar(x=p, height=h, width=w, color='silver')
avenir_font = {'fontname':'Avenir'}
ax.set_yticks([0,50000,100000])
ax.set_yticklabels(['0%','50%','100%'],  **avenir_font)
ax.set_ylabel('Proportion HGT', labelpad=10, fontsize=20, **avenir_font)
#ax.get_yaxis().set_visible(False)
ax.set_xlabel('Genomic position', labelpad=5, fontsize=40, **avenir_font)
ax.set_xticklabels(['0 Kb','200 Kb', '400 Kb', '600 Kb', '800 Kb', '1000 Kb', '1200 Kb', '1400 Kb', '1600 Kb'],  **avenir_font)
ax.tick_params(axis='x', which='major', labelsize=20, labelrotation=45, pad=7)
plt.xticks(rotation=40, ha='right')
ax.tick_params(axis='y', which='major', labelsize=20)
ax.spines["top"].set_color("#66FFCB")
ax.spines["top"].set_linewidth(4)
ax.spines["left"].set_color("#66FFCB")
ax.spines["left"].set_linewidth(4)
ax.spines["right"].set_color("#66FFCB")
ax.spines["right"].set_linewidth(4)
ax.spines["bottom"].set_color("#66FFCB")
ax.spines["bottom"].set_linewidth(4)

"""
# For colouring by driver mutation type
#ax0.bar(x=positions, height=HGT_freq_list, width=widths, color='silver')
ax0.spines['top'].set_visible(False)
ax0.spines['right'].set_visible(False)
ax0.spines['left'].set_visible(False)
ax0.spines['bottom'].set_visible(False) 
ax0.tick_params(axis='both', bottom=False, left=False)

plt.setp(ax0.get_xticklabels(), visible=False)


fig.subplots_adjust(hspace=0, wspace=0)

#x, y = np.mgrid[-5:5:0.05, -5:5:0.05]
#z = (np.sqrt(x**2 + y**2) + np.sin(x**2 + y**2))
#im = ax.imshow(z, cmap=HGT_cmap)
data = np.arange(100, 0, -1).reshape(10, 10)

im = ax.imshow(data, cmap=HGT_cmap)
im2 = ax.imshow(data, cmap=DNp5_cmap)
im3 = ax.imshow(data, cmap=anc_cmap)
im4 = ax.imshow(data, cmap=misc_cmap)

fig.colorbar(im)
fig.colorbar(im2)
fig.colorbar(im3)
fig.colorbar(im4)
"""
plt.gcf().subplots_adjust(bottom=0.4)

#plt.savefig(op_path+npop+'_coloured_by_SC.png', format='png', dpi=600)
plt.show()
