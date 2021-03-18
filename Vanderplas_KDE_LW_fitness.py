#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 16:04:30 2021

"""

# Authors: Jake Vanderplas <jakevdp@cs.washington.edu>
#          Mehreen Saeed https://stackabuse.com/author/mehreen/
#
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import norm
from sklearn.neighbors import KernelDensity
from sklearn.utils.fixes import parse_version
from sklearn.model_selection import GridSearchCV


import pandas
import openpyxl
#Open data file
"""
sel_coefs_file='/Users/lwoo0005/Documents/An_H_pylori/All_MCE3_minD0_selcoefs_combined.xlsx'
sheet1='MCE3_sig_Ab-_minD0_selcoefs'
sheet2='partHGT_MCE3_sig_Ab-_selcoefs'

X_full=(pandas.read_excel(sel_coefs_file, sheet_name=sheet1, header = 0, names=['Coefficient'],usecols=[7])).dropna()
X_partial=(pandas.read_excel(sel_coefs_file, sheet_name=sheet2, header = 0, names=['Coefficient'],usecols=[7])).dropna()
"""
sel_coefs_file_1='/Users/lwoo0005/Documents/An_H_pylori/Shami_stuff/MCE3_allvars_rm1tp_selcoefs_and_blocks.xlsx'
sel_coefs_file_2='/Users/lwoo0005/Documents/An_H_pylori/Shami_stuff/partMCE3_allvars_selection_coefs.xlsx'

sheet1='MCE3_allvars_Ab-_rm1tp_blocks'
sheet2='MCE3_partAb-_allvars_rm1tp_bloc'

X_full=(pandas.read_excel(sel_coefs_file_1, sheet_name=sheet1, header = 0, names=['Average_SC'],usecols=[8])).dropna()
X_partial=(pandas.read_excel(sel_coefs_file_2, sheet_name=sheet2, header = 0, names=['Average_SC'],usecols=[8])).dropna()


def my_scores(estimator, X):
    scores = estimator.score_samples(X)
    # Remove -inf
    scores = scores[scores != float('-inf')]
    # Return the mean values
    return np.mean(scores)

h_vals=np.arange(0.05, 0.11, .01)
kernels= ['gaussian', 'linear', 'tophat']
X_plot = np.linspace(-1, 1, 5000)[:, np.newaxis]

grid = GridSearchCV(KernelDensity(),
                    {'bandwidth': h_vals, 'kernel': kernels},
                    scoring=my_scores
                    )
log_dens_list=[]
best_kdes=[]
for X in [X_partial, X_full]:
    grid.fit(X)
    best_kde = grid.best_estimator_
    best_kdes.append(best_kde)
    log_dens = best_kde.score_samples(X_plot)
    log_dens_list.append(log_dens)

from matplotlib import rcParams
from matplotlib.ticker import FormatStrFormatter
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Avenir']

fig = plt.figure()
ax = fig.add_subplot(111)
rgbs = [(102/255,255/255,203/255,0.5),(252/255,102/255,254/255,0.5)]
for name,n,X, log_dens, c, best_kde, r in zip(['Foundational','Expanded'],[5,3],[X_partial, X_full],log_dens_list,['#66FFCB','#FC66FE'],best_kdes, rgbs):
    ax.fill(X_plot[:,0], np.exp(log_dens),fc=r, ec=c, lw=2)
    #ax.text(-1,n,"Best Kernel for "+name+":\n" + best_kde.kernel+", h="+"{:.3f}".format(best_kde.bandwidth),fontsize=10)
    #plt.plot(X, np.full_like(X, -0.5), c, markeredgewidth=0.01,alpha=0.25)
#ax.legend(['Foundational','Expanded'])
ax.set_xlim(-0.3,0.3)
ax.set_ylim(0,10)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.tick_params(axis='both', labelsize=15)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.xticks(rotation=45)
plt.savefig('/Users/lwoo0005/Documents/An_H_pylori/MCE3_Ab-_sig_rm1tp_BLOCK_SCs_partial2full_KDE.png',format='png',dpi=600)
plt.show()

"""
#kde = KernelDensity(kernel='linear', bandwidth=0.02).fit(X)
fig, ax = plt.subplots(1, 1, sharex=True, sharey=True)

#log_dens = kde.score_samples(X_plot)
ax.fill(X_plot[:,0], np.exp(log_dens), fc='#AAAAFF')
ax.text(0.2, 0.8, "Linear Kernel Density")
"""
