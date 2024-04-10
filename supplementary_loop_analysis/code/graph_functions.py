import matplotlib.image as mpimg
import matplotlib.cm as cm
import matplotlib.colors as colors
import pyBigWig
import subprocess
import numpy as np
import cooler
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3

from analysis_functions import *


chromo_keys = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrM', 'chrX', 'chrY']
chromos = {}
for chromo in chromo_keys:
    with open(f'./dm6/{chromo}.txt', 'r') as f:
        chromos[chromo] = f.readline()

# Plotting a specific network as a heatmap
def plot_network_heatmap(network, ind2anc, final_loops, cool):
    ancsoi = [ind2anc[n] for n in network.nodes]
    sort_ancs = sorted(ancsoi)

    start, end = sort_ancs[0][1], sort_ancs[-1][1]
    chrom = sort_ancs[0][0]

    plot(final_loops, chrom, start, end, cool, show = True)


import os
import shutil
def plot_networks(test_loopset, subgraph,  node_colors, adjacency, ind2anc, anc2ind, cool, vmin = 0, vmax = 17500, tads=None, save = False, save_dir = '../standard_plots', dpi = 100, show=True):
    fig, axs = plt.subplots(1, 2, figsize=(15, 6), dpi=dpi)
    cmap = plt.cm.bwr
    num = 5

    ax1 = plt.subplot2grid((100,100), (50, 50), colspan = 50, rowspan = 50, fig = fig)
    
    ax2 = plt.subplot2grid((100,100), (0, 50), colspan = 50, rowspan = 2, fig = fig)
    ax3 = plt.subplot2grid((100,100), (3, 50), colspan = 50, rowspan = 2, fig = fig)
    ax4 = plt.subplot2grid((100,100), (6, 50), colspan = 50, rowspan = 2, fig = fig)
    ax5 = plt.subplot2grid((100,100), (9, 50), colspan = 50, rowspan = 2, fig = fig)
    ax6 = plt.subplot2grid((100,100), (12, 50), colspan = 50, rowspan = 2, fig = fig)
    ax7 = plt.subplot2grid((100,100), (15, 50), colspan = 50, rowspan = 2, fig = fig)
    ax8 = plt.subplot2grid((100,100), (18, 50), colspan = 50, rowspan = 2, fig = fig)
    ax9 = plt.subplot2grid((100,100), (21, 50), colspan = 50, rowspan = 2, fig = fig)
    ax10 = plt.subplot2grid((100,100), (25, 50), colspan = 50, rowspan = 2, fig = fig)
    ax11 = plt.subplot2grid((100,100), (28, 50), colspan = 50, rowspan = 2, fig = fig)
    ax12 = plt.subplot2grid((100,100), (31, 50), colspan = 50, rowspan = 2, fig = fig)
    ax13 = plt.subplot2grid((100,100), (34, 50), colspan = 50, rowspan = 2, fig = fig)
    ax14 = plt.subplot2grid((100,100), (37, 50), colspan = 50, rowspan = 2, fig = fig)
    ax15 = plt.subplot2grid((100,100), (40, 50), colspan = 50, rowspan = 2, fig = fig)
    ax16 = plt.subplot2grid((100,100), (43, 50), colspan = 50, rowspan = 2, fig = fig)
    ax17 = plt.subplot2grid((100,100), (46, 50), colspan = 50, rowspan = 2, fig = fig)


    axoos = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12, ax13, ax14, ax15, ax16, ax17, ]

    nx.draw(subgraph, node_color = node_colors, ax=axs[0], cmap = cmap, vmin=vmin, vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin = vmin, vmax=vmax))
    fig.colorbar(sm, ax=axs[0])

    min1, max2 = 0, 0
    ancsoi = [ind2anc[l] for l in list(subgraph.nodes)]

    all_loops = {}
    loopsoi = {}
    for anc in ancsoi:
        ind = anc2ind[anc]
        connec = np.where(adjacency[ind, :] > 0)[0]
        for c in connec:
            anc2 = ind2anc[c]
            l1, l2 = sorted((anc, anc2))
            all_loops[l1, l2] = 1
            if anc2ind[l1] in subgraph.nodes and anc2ind[l2] in subgraph.nodes:
                loopsoi[l1, l2] = 1

    chrom = l1[0]
    otherloops = {}
    for loop in loopsoi:
        l1, l2 = loop
        if min1 > l1[1] or min1 == 0:
            min1 = l1[1]
        if max2 < l2[1] or max2 == 0:
            max2 = l2[1]
        chrom = l1[0]
        otherloops[loop]=1



    axs[0].set_title(f"{chrom}, {min1-40*100}, {max2+40*100}")
    print(f"{chrom}, {min1-40*100}, {max2+40*100}")
    anc_dict = {}
    anc_dict[chrom] = []
    for anc in ancsoi:
        anc_dict[chrom].append(anc)
    anc_dict[chrom].sort()
    t = np.argsort(list(zip(*ancsoi))[1])
    sorted_node_colors = np.asarray(node_colors)[t]

    looplist = [otherloop_to_loop(loopsoi), otherloop_to_loop(all_loops)]

    if tads is None:
        plot_looplist(looplist, chrom, min1-40*100, max2+40*100, cool, axs = axoos, 
            tads=anc_dict, tad_colors = sorted_node_colors, vmax=vmax, save = save, save_dir = save_dir, show=show)
    else:
        plot_looplist(looplist, chrom, min1-40*100, max2+40*100, cool, axs = axoos, 
            tads=tads,  vmax=vmax, save = save, save_dir = save_dir, show=show)


    for ax in axoos[2:]:
        bottom, top = ax.get_ylim()
        ntop = min(50, top)
        ax.set_ylim([0, ntop])
        ax.set_yticks([ntop])

    if save == True:
        fig.savefig(f"{save_dir}{chrom}_{min1-40*100}_{max2+40*100}.png")


def LINEAR_plot_networks(test_loopset, subgraph,  node_colors, adjacency, ind2anc, anc2ind, cool, vmin = 0, vmax = 17500, tads=None, save = False, save_dir = '../standard_plots', dpi = 100, show=True):
    fig, axs = plt.subplots(1, 2, figsize=(15, 6), dpi=dpi)
    cmap = plt.cm.bwr
    num = 5

    ax1 = plt.subplot2grid((100,100), (50, 50), colspan = 50, rowspan = 50, fig = fig)
    
    ax2 = plt.subplot2grid((100,100), (0, 50), colspan = 50, rowspan = 2, fig = fig)
    ax3 = plt.subplot2grid((100,100), (3, 50), colspan = 50, rowspan = 2, fig = fig)
    ax4 = plt.subplot2grid((100,100), (6, 50), colspan = 50, rowspan = 2, fig = fig)
    ax5 = plt.subplot2grid((100,100), (9, 50), colspan = 50, rowspan = 2, fig = fig)
    ax6 = plt.subplot2grid((100,100), (12, 50), colspan = 50, rowspan = 2, fig = fig)
    ax7 = plt.subplot2grid((100,100), (15, 50), colspan = 50, rowspan = 2, fig = fig)
    ax8 = plt.subplot2grid((100,100), (18, 50), colspan = 50, rowspan = 2, fig = fig)
    ax9 = plt.subplot2grid((100,100), (21, 50), colspan = 50, rowspan = 2, fig = fig)
    ax10 = plt.subplot2grid((100,100), (25, 50), colspan = 50, rowspan = 2, fig = fig)
    ax11 = plt.subplot2grid((100,100), (28, 50), colspan = 50, rowspan = 2, fig = fig)
    ax12 = plt.subplot2grid((100,100), (31, 50), colspan = 50, rowspan = 2, fig = fig)
    ax13 = plt.subplot2grid((100,100), (34, 50), colspan = 50, rowspan = 2, fig = fig)
    ax14 = plt.subplot2grid((100,100), (37, 50), colspan = 50, rowspan = 2, fig = fig)
    ax15 = plt.subplot2grid((100,100), (40, 50), colspan = 50, rowspan = 2, fig = fig)
    ax16 = plt.subplot2grid((100,100), (43, 50), colspan = 50, rowspan = 2, fig = fig)
    ax17 = plt.subplot2grid((100,100), (46, 50), colspan = 50, rowspan = 2, fig = fig)


    axoos = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12, ax13, ax14, ax15, ax16, ax17, ]
    n_nodes = len(subgraph.nodes)
    positions = {}
    for i, node in enumerate(subgraph.nodes):
        positions[node] = [0, i/n_nodes]
    nx.draw(subgraph, node_color = node_colors, ax=axs[0], cmap = cmap, vmin=vmin, vmax=vmax, pos=positions)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin = vmin, vmax=vmax))
    fig.colorbar(sm, ax=axs[0])

    min1, max2 = 0, 0
    ancsoi = [ind2anc[l] for l in list(subgraph.nodes)]

    all_loops = {}
    loopsoi = {}
    for anc in ancsoi:
        ind = anc2ind[anc]
        connec = np.where(adjacency[ind, :] > 0)[0]
        for c in connec:
            anc2 = ind2anc[c]
            l1, l2 = sorted((anc, anc2))
            all_loops[l1, l2] = 1
            if anc2ind[l1] in subgraph.nodes and anc2ind[l2] in subgraph.nodes:
                loopsoi[l1, l2] = 1

    chrom = l1[0]
    otherloops = {}
    for loop in loopsoi:
        l1, l2 = loop
        if min1 > l1[1] or min1 == 0:
            min1 = l1[1]
        if max2 < l2[1] or max2 == 0:
            max2 = l2[1]
        chrom = l1[0]
        otherloops[loop]=1



    axs[0].set_title(f"{chrom}, {min1-40*100}, {max2+40*100}")
    print(f"{chrom}, {min1-40*100}, {max2+40*100}")
    anc_dict = {}
    anc_dict[chrom] = []
    for anc in ancsoi:
        anc_dict[chrom].append(anc)
    anc_dict[chrom].sort()
    t = np.argsort(list(zip(*ancsoi))[1])
    sorted_node_colors = np.asarray(node_colors)[t]

    looplist = [otherloop_to_loop(loopsoi), otherloop_to_loop(all_loops)]

    if tads is None:
        plot_looplist(looplist, chrom, min1-40*100, max2+40*100, cool, axs = axoos, 
            tads=anc_dict, tad_colors = sorted_node_colors, vmax=vmax, save = save, save_dir = save_dir, show=show)
    else:
        plot_looplist(looplist, chrom, min1-40*100, max2+40*100, cool, axs = axoos, 
            tads=tads,  vmax=vmax, save = save, save_dir = save_dir, show=show)


    for ax in axoos[2:]:
        bottom, top = ax.get_ylim()
        ntop = min(50, top)
        ax.set_ylim([0, ntop])
        ax.set_yticks([ntop])

    if save == True:
        fig.savefig(f"{save_dir}{chrom}_{min1-40*100}_{max2+40*100}.png")





def plot_heatmap_with_network(loopset, chrom, start, end, big_graph, ind2anc, anc2ind, cool, vmin = 0, vmax = 17500, save = False, save_dir = '../standard_plots', dpi = 100, show=True):
    fig, axs = plt.subplots(1, 2, figsize=(15, 6), dpi=dpi)
    cmap = plt.cm.bwr
    num = 5
    ax1 = plt.subplot2grid((100,100), (35, 50), colspan = 50, rowspan = 70, fig = fig)
    ax2 = plt.subplot2grid((100,100), (0, 50), colspan = 50, rowspan = 8, fig = fig)
    ax3 = plt.subplot2grid((100,100), (10, 50), colspan = 50, rowspan = 4, fig = fig)
    ax4 = plt.subplot2grid((100,100), (15, 50), colspan = 50, rowspan = 4, fig = fig)
    ax5 = plt.subplot2grid((100,100), (20, 50), colspan = 50, rowspan = 4, fig = fig)
    ax6 = plt.subplot2grid((100,100), (25, 50), colspan = 50, rowspan = 4, fig = fig)
    ax7 = plt.subplot2grid((100,100), (30, 50), colspan = 50, rowspan = 4, fig = fig)

    axoos = [ax1, ax2, ax3, ax4, ax5, ax6, ax7]

    loopsoi = {}
    ancsoi = {}
    for key in loopset.other_loops.keys():
        l1, l2 = key
        if l1[0] != chrom:
            continue
        ss = l1[1]
        ee = l2[1]
        if ss > start and ee < end:
            loopsoi[l1, l2] = 1
            ancsoi[l1] = 1
            ancsoi[l2] = 1


    indsoi = [anc2ind[anc] for anc in ancsoi]
    print(indsoi)
    subgraph = big_graph.subgraph(indsoi)


    node_colors = np.arange(0, len(indsoi))

    nx.draw(subgraph, node_color = node_colors, ax=axs[0], cmap = cmap, vmin=vmin, vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin = vmin, vmax=vmax))
    fig.colorbar(sm, ax=axs[0])

    axs[0].set_title(f"{chrom}, {start}, {end}")
    anc_dict = {}
    anc_dict[chrom] = []
    for anc in ancsoi:
        anc_dict[chrom].append(anc)
    anc_dict[chrom].sort()
    t = np.argsort(list(zip(*ancsoi))[1])
    sorted_node_colors = np.asarray(node_colors)[t]

    loopsoi = otherloop_to_loop(loopsoi)
    plot(loopsoi, chrom, start, end, cool, axs = axoos, 
        tads=ancsoi, tad_colors = sorted_node_colors, vmax=vmax, save = save, save_dir = save_dir, show=show)

    for ax in axoos[2:]:
        bottom, top = ax.get_ylim()
        ntop = min(50, top)
        ax.set_ylim([0, ntop])
        ax.set_yticks([ntop])

    if save == True:
        fig.savefig(f"{save_dir}{chrom}_{min1-40*100}_{max2+40*100}.png")        

def plot_venn(dic1, dic2, l1, l2, title, ax):
    overlap = len(get_intersection(dic1, dic2))
#     neither = len(get_neither(dic1, dic2))
    left = len(get_a_not_b(dic1, dic2))
    right = len(get_a_not_b(dic2, dic1))
    venn2((left, right, overlap), set_labels = (l1, l2), ax=ax)
    ax.set_title(title)
    
def plot_venn3(dic1, dic2, dic3, l1, l2, l3, title, ax):
    set1, set2, set3 = set(), set(), set()
    
    for key in dic1.keys():
        v = dic1.get(key)
        if v:
            set1.add(key)
    for key in dic2.keys():
        v = dic2.get(key)
        if v:
            set2.add(key)
    for key in dic3.keys():
        v = dic3.get(key)
        if v:
            set3.add(key)
    venn3((set1, set2, set3), set_labels = (l1, l2, l3), ax=ax)
    ax.set_title(title)

def venn_diagrams(loopset, rand_loopset):
    name = loopset.name
    randname = rand_loopset.name
    fig, axs = plt.subplots(4, 4, figsize=(15, 10))
    plot_venn(loopset.peak_cp190, loopset.anchors, 'CP190', f'Loop anchors', '', ax=axs[0][0])
    plot_venn(loopset.peak_trl2, loopset.anchors, 'GAF ', f'Loop anchors', '', ax=axs[0][1])
    plot_venn(loopset.peak_pc2, loopset.anchors, 'PC ', f'Loop anchors', '', ax=axs[0][2])
    plot_venn(loopset.peak_pser, loopset.anchors, 'pol2 ', f'Loop anchors', '', ax=axs[0][3])


    plot_venn(rand_loopset.peak_cp190, rand_loopset.anchors, 'CP190 ', f'Shifted loop anchors', '', ax=axs[1][0])
    plot_venn(rand_loopset.peak_trl2, rand_loopset.anchors, 'GAF ', f'Shifted loop anchors', '', ax=axs[1][1])
    plot_venn(rand_loopset.peak_pc2, rand_loopset.anchors, 'PC', f'Shifted loop anchors', '', ax=axs[1][2])
    plot_venn(rand_loopset.peak_pol22, rand_loopset.anchors, 'pol2 ', f'Shifted loop anchors', '', ax=axs[1][3])


    plot_venn3(loopset.peak_pol22, loopset.peak_beaf,  loopset.peak_cp190, 'pol2 ', 'beaf32 ', 'cp190', '', ax=axs[2][0])
    plot_venn3(loopset.peak_pc2, loopset.peak_trl2,  loopset.peak_cp190, 'pc ', 'GAF ', 'cp190', '', ax=axs[2][1])
    plot_venn3(loopset.peak_pol22, loopset.peak_trl2, loopset.peak_pc2, 'pol2 ', 'GAF ', 'pc', '', ax=axs[2][2])
    plot_venn3(loopset.peak_pol22, loopset.peak_trl2,  loopset.peak_beaf, 'pol2 ', 'GAF ', 'beaf32', '', ax=axs[2][3])

    
    plot_venn3(rand_loopset.peak_pol22, rand_loopset.peak_beaf, rand_loopset.peak_cp190, 'pol2 ', 'beaf32 ', 'cp190', '', ax=axs[3][0])
    plot_venn3(rand_loopset.peak_pc2, rand_loopset.peak_trl2, rand_loopset.peak_cp190, 'pc', 'GAF ', 'cp190', '', ax=axs[3][1])
    plot_venn3(rand_loopset.peak_pol22, rand_loopset.peak_trl2, rand_loopset.peak_pc2, 'pol2 ', 'GAF ', 'pc', '', ax=axs[3][2])
    plot_venn3(rand_loopset.peak_pol22, rand_loopset.peak_trl2, rand_loopset.peak_beaf, 'pol2 ', 'GAF ', 'beaf32', '', ax=axs[3][3])

    
def loop_diagrams(loopset, rand_loopset):
    name = loopset.name
    fig, axs = plt.subplots(6, 4, figsize=(20, 20))    
    
    plot_venn3(loopset.one_trl_peak, loopset.both_trl_peak, loopset.other_loops, '1 TRL PEAK', '2 TRL PEAK', f'LOOPS', '', ax=axs[0][0])
    plot_venn3(loopset.one_pser_peak, loopset.both_pser_peak, loopset.other_loops, '1 PSER', '2 PSER ', f'LOOPS', '', ax=axs[0][1])
    plot_venn3(loopset.one_ctcf_peak, loopset.both_ctcf_peak, loopset.other_loops, '1 CTCF', '2 CTCF', f'LOOPS', '', ax=axs[0][2])

    plot_venn3(rand_loopset.one_trl_peak, rand_loopset.both_trl_peak, rand_loopset.other_loops, '1 TRL PEAK', '2 TRL PEAK', f'LOOPS', '', ax=axs[1][0])
    plot_venn3(rand_loopset.one_pser_peak, rand_loopset.both_pser_peak, rand_loopset.other_loops, '1 PSER PEAK', '2 PSER PEAK', f'LOOPS', '', ax=axs[1][1])
    plot_venn3(rand_loopset.one_ctcf_peak, rand_loopset.both_ctcf_peak, rand_loopset.other_loops, '1 CTCF peak', '2 CTCF peak', f'LOOPS', '', ax=axs[1][2])
    
    plot_venn3(loopset.one_trl_peak, loopset.both_trl_peak, loopset.both_pser_peak, '1 TRL PEAK', '2 TRL PEAK', f'2 PSER', '', ax=axs[2][0])
    plot_venn3(loopset.one_trl_peak, loopset.both_trl_peak, loopset.one_pser_peak, '1 TRL PEAK', '2 TRL PEAK', f'1 PSER', '', ax=axs[2][1])
    plot_venn3(loopset.one_trl_peak, loopset.both_pser_peak, loopset.one_pser_peak, '1 TRL PEAK', '2 PSER PEAK', f'1 PSER', '', ax=axs[2][2])    
    plot_venn3(loopset.both_trl_peak, loopset.both_pser_peak, loopset.one_pser_peak, '2 TRL PEAK', '2 PSER PEAK', f'1 PSER', '', ax=axs[2][3])        

    plot_venn3(rand_loopset.one_trl_peak, rand_loopset.both_trl_peak, rand_loopset.both_pser_peak, '1 TRL PEAK', '2 TRL PEAK', f'2 PSER', '', ax=axs[3][0])
    plot_venn3(rand_loopset.one_trl_peak, rand_loopset.both_trl_peak, rand_loopset.one_pser_peak, '1 TRL PEAK', '2 TRL PEAK', f'1 PSER', '', ax=axs[3][1])
    plot_venn3(rand_loopset.one_trl_peak, rand_loopset.both_pser_peak, rand_loopset.one_pser_peak, '1 TRL PEAK', '2 PSER PEAK', f'1 PSER', '', ax=axs[3][2])    
    plot_venn3(rand_loopset.both_trl_peak, rand_loopset.both_pser_peak, rand_loopset.one_pser_peak, '2 TRL PEAK', '2 PSER PEAK', f'1 PSER', '', ax=axs[3][3])        
    
    plot_venn3(loopset.both_trl_peak, loopset.both_pser_peak, loopset.both_ctcf_peak, '2 TRL PEAK', '2 PSER PEAK', f'2 CTCF peak', '', ax=axs[4][1])
    plot_venn3(loopset.one_ctcf_peak, loopset.both_ctcf_peak, loopset.one_pser_peak, '1 CTCF peak', '2 CTCF peak', f'1 PSER', '', ax=axs[4][2])

    plot_venn3(rand_loopset.both_trl_peak, rand_loopset.both_pser_peak, rand_loopset.both_ctcf_peak, '2 TRL PEAK', '2 PSER PEAK', f'2 CTCF peak', '', ax=axs[5][1])
    plot_venn3(rand_loopset.one_ctcf_peak, rand_loopset.both_ctcf_peak, rand_loopset.one_pser_peak, '1 CTCF peak', '2 CTCF peak', f'1 PSER', '', ax=axs[5][2])


def plot_anchor_pileup(bw_vals, real_vals, rand_vals):
    fig, axs = plt.subplots(figsize = (5, 3), dpi=200)

    ax2 = plt.subplot2grid((100,100), (0, 0), colspan = 40, rowspan = 4, fig = fig)
    ax3 = plt.subplot2grid((100,100), (10, 0), colspan = 40, rowspan = 4, fig = fig)
    ax4 = plt.subplot2grid((100,100), (15, 0), colspan = 40, rowspan = 4, fig = fig)
    ax5 = plt.subplot2grid((100,100), (20, 0), colspan = 40, rowspan = 4, fig = fig)
    ax6 = plt.subplot2grid((100,100), (5, 0), colspan = 40, rowspan = 4, fig = fig)
    ax7 = plt.subplot2grid((100,100), (25, 0), colspan = 40, rowspan = 70, fig = fig)

    ax7.imshow(np.log(real_vals/rand_vals)) 


    for ax, name in zip([ax2, ax3, ax4, ax5, ax6], ['ctcf', 'trl', 'pser', 'cp190', "zelda",]):
        bw_val = bw_vals[name]
        ax.plot(np.arange(len(bw_val)), bw_val, linewidth = .8)
        ax.set_xticks([])
        ax.set_yticks([int(np.max(bw_val))])
        # ax.yaxis.set_label_position("right")
        h = ax.set_ylabel(name, fontsize=4, rotation=0, labelpad=15)
        ax.yaxis.tick_right()
        ax.tick_params(axis='both', which='major', labelsize=6)    




def plot_loop_pileup(L_bw_vals, R_bw_vals, real_vals, rand_vals):
    fig, ax = plt.subplots(figsize = (5, 3), dpi=200)

    ax1 = plt.subplot2grid((100,100), (35, 0), colspan = 40, rowspan = 70, fig = fig)

    ax2 = plt.subplot2grid((100,100), (0, 0), colspan = 40, rowspan = 4, fig = fig)
    ax3 = plt.subplot2grid((100,100), (10, 0), colspan = 40, rowspan = 4, fig = fig)
    ax4 = plt.subplot2grid((100,100), (15, 0), colspan = 40, rowspan = 4, fig = fig)
    ax5 = plt.subplot2grid((100,100), (20, 0), colspan = 40, rowspan = 4, fig = fig)
    ax6 = plt.subplot2grid((100,100), (5, 0), colspan = 40, rowspan = 4, fig = fig)
    ax7 = plt.subplot2grid((100,100), (25, 0), colspan = 40, rowspan = 70, fig = fig)

    ax1.imshow(np.log(real_vals/rand_vals))

    for ax, name in zip([ax2, ax3, ax4, ax5, ax6], ['ctcf', 'trl', 'pser', 'cp190', "zelda",]):
        L_bw_val = L_bw_vals[name]
        R_bw_val = R_bw_vals[name]
        print(L_bw_val.shape)
        ax.plot(np.arange(len(L_bw_val)), L_bw_val, linewidth = .8)
        ax.plot(np.arange(len(R_bw_val)), R_bw_val, linewidth = .8)
        ax.set_xticks([])
        ax.set_yticks([int(max(np.max(L_bw_val), np.max(R_bw_val)))])
        
        # ax.yaxis.set_label_position("right")
        h = ax.set_ylabel(name, fontsize=4, rotation=0, labelpad=15)
        ax.yaxis.tick_right()
        ax.tick_params(axis='both', which='major', labelsize=6)    

def scatter_loops(df, prot1, prot2, indices=None, cs='blue', colorname=""):
    if indices is None:
        indices = np.asarray(df.index)
    fig, axs = plt.subplots(1, 2, figsize= (10, 5))
    ax = axs[0]
    v1 = df[prot1][indices]
    v2 = df[prot2][indices]
    colors = cs[indices]

    rand_indices = np.random.choice(range(len(df)), len(indices))
    rand_colors = cs[rand_indices]

    vmin, vmax = np.min(np.min(colors)), np.max(np.max(colors))

    pos = ax.scatter(v1, v2, c = colors, alpha = .7, vmin=vmin, vmax = vmax)
    ax.set_title(f"{prot1} co-expression over selected loops. \nColor: {colorname}")
    ax.set_xlabel(f"{prot1}, anchor 1 (Normalized bigwig)")
    ax.set_ylabel(f"{prot2}, anchor 2 (Normalized bigwig)")
    ax.set_ylim([-.1, 1.1])
    ax.set_xlim([-.1, 1.1])
    # plt.colorbar(pos)


    ax = axs[1]
    v1 = df[prot1].values[rand_indices]
    v2 = df[prot2].values[rand_indices]
    pos = ax.scatter(v1, v2, c = rand_colors, alpha = .7, vmin = vmin, vmax = vmax)
    ax.set_title(f"{prot1} co-expression over randomly selected loops. \nColor: {colorname}")
    ax.set_xlabel(f"{prot1}, anchor 1 (Normalized bigwig)")
    ax.set_ylabel(f"{prot2}, anchor 2 (Normalized bigwig)")
    ax.set_ylim([-.1, 1.1])
    ax.set_xlim([-.1, 1.1])

    fig.subplots_adjust(right=0.95)
    cbar_ax = fig.add_axes([0.97, 0.15, 0.025, 0.7])
    fig.colorbar(pos, cax=cbar_ax)


def scatter_anchors(df, prot1, prot2, indices=None, cs='blue', colorname=""):
    if indices is None:
        indices = np.asarray(df.index)
        rand_indices = np.asarray(df.index)
    else:
        rand_indices = np.random.choice(range(len(df)), len(indices))


    if cs != 'blue':
        colors = cs[indices]
        vmin, vmax = np.min(np.min(colors)), np.max(np.max(colors))
        rand_colors = cs[rand_indices]
    else: 
        colors = 'blue'
        rand_colors = 'blue'
        vmin = 0
        vmax = 1

    fig, axs = plt.subplots(1, 2, figsize= (10, 5))        
    ax = axs[0]
    v1 = df[prot1].values[indices]
    v2 = df[prot2].values[indices]

    pos = ax.scatter(v1, v2, c = colors, alpha = .7, vmin=vmin, vmax = vmax)
    ax.set_title(f"{prot1} co-expression over selected loops. \nColor: {colorname}")
    ax.set_xlabel(f"{prot1}, anchor (Normalized bigwig)")
    ax.set_ylabel(f"{prot2}, anchor (Normalized bigwig)")
    ax.scatter(np.mean(v1), np.mean(v2), c='red')
    # ax.set_ylim([-.1, 1.1])
    # ax.set_xlim([-.1, 1.1])
    # plt.colorbar(pos)


    ax = axs[1]
    v1 = df[prot1].values[rand_indices]
    v2 = df[prot2].values[rand_indices]
    pos = ax.scatter(v1, v2, c = rand_colors, alpha = .7, vmin = vmin, vmax = vmax)
    ax.set_title(f"{prot1} co-expression over randomly selected loops. \nColor: {colorname}")
    ax.set_xlabel(f"{prot1}, anchor (Normalized bigwig)")
    ax.set_ylabel(f"{prot2}, anchor (Normalized bigwig)")
    # ax.set_ylim([-.1, 1.1])
    # ax.set_xlim([-.1, 1.1])
    ax.scatter(np.mean(v1), np.mean(v2), c='red')

    fig.subplots_adjust(right=0.95)
    cbar_ax = fig.add_axes([0.97, 0.15, 0.025, 0.7])
    fig.colorbar(pos, cax=cbar_ax)


import seaborn as sns
def side_heatmap(dfs, name, order, vmax_keys, mask, vmax):
    for key in dfs:
        fig, axs = plt.subplots(1, 2, figsize = (15, 5))
        sns.heatmap(dfs[name].iloc[order, :].loc[mask, :], vmax=vmax, ax=axs[0], yticklabels=False)
        axs[0].set_title(name)
        axs[1].set_title(key)
        sns.heatmap(dfs[key].iloc[order, :].loc[mask, :], vmax=vmax_keys.get(key, 25), ax=axs[1], yticklabels=False)
        fig.suptitle(key)




def plot(loops, chrom, start, end, c, tads = {}, plotLoops = True, balance=True, show = False, save = False, save_dir = '../standard_plots', show_tads = True, tad_colors = None, cmap = plt.cm.bwr, vmin = 0, vmax = 17500, axs = None):
    if start < 0 or end > len(chromos[chrom]):
        print("Those coordinates are not valid!")
        return False
    res = c.binsize
    mat = np.rot90(c.matrix(balance=balance).fetch((chrom, start, end)))

    matmax = np.max(mat[~np.isnan(mat)])
    if axs == None:
        fig, ax = plt.subplots(figsize = (10, 10), dpi=200)

        ax2 = plt.subplot2grid((100,100), (0, 0), colspan = 90, rowspan = 2, fig = fig)
        ax3 = plt.subplot2grid((100,100), (3, 0), colspan = 90, rowspan = 2, fig = fig)
        ax4 = plt.subplot2grid((100,100), (6, 0), colspan = 90, rowspan = 2, fig = fig)
        ax5 = plt.subplot2grid((100,100), (9, 0), colspan = 90, rowspan = 2, fig = fig)
        ax6 = plt.subplot2grid((100,100), (12, 0), colspan = 90, rowspan = 2, fig = fig)
        ax7 = plt.subplot2grid((100,100), (15, 0), colspan = 90, rowspan = 2, fig = fig)
        ax8 = plt.subplot2grid((100,100), (18, 0), colspan = 90, rowspan = 2, fig = fig)
        ax9 = plt.subplot2grid((100,100), (21, 0), colspan = 90, rowspan = 2, fig = fig)
        ax10 = plt.subplot2grid((100,100), (25, 0), colspan = 90, rowspan = 2, fig = fig)
        ax11 = plt.subplot2grid((100,100), (28, 0), colspan = 90, rowspan = 2, fig = fig)
        ax12 = plt.subplot2grid((100,100), (31, 0), colspan = 90, rowspan = 2, fig = fig)
        ax13 = plt.subplot2grid((100,100), (34, 0), colspan = 90, rowspan = 2, fig = fig)
        ax14 = plt.subplot2grid((100,100), (37, 0), colspan = 90, rowspan = 2, fig = fig)
        ax15 = plt.subplot2grid((100,100), (40, 0), colspan = 90, rowspan = 2, fig = fig)
        ax16 = plt.subplot2grid((100,100), (43, 0), colspan = 90, rowspan = 2, fig = fig)
        ax17 = plt.subplot2grid((100,100), (46, 0), colspan = 90, rowspan = 2, fig = fig)

        ax1 = plt.subplot2grid((100,100), (50, 0), colspan = 90, rowspan = 50, fig = fig)

    else:
        fig = 1
        save = False
        ax1, ax2, ax3, ax4, ax5, ax6, ax7 = axs

    ax2.set_xticks([])
    ax2.set_yticks([])


    paths = ["../tfs/bws/GSE152770_GAF_GFP_COMBINED.bw", "../tfs/bws/CP190.bw", "../tfs/bws/GSE62925_pol2_COMBINED.bw", 
            '../tfs/bws/ctcf.bw', '../tfs/bws/zelda.bw', '../tfs/bws/BEAF32_SRR1636749.bw', '../tfs/bws/DREF_SRR1636771.bw', 
            '../tfs/bws/GSE60428_PC_COMBINED.bw', '../tfs/bws/GSE60428_PH_COMBINED.bw', '../tfs/bws/rad21_COMBINED.bw', 
            '../histones/bws/H3K27me3.bw', '../histones/bws/H3K27ac.bw', '../tfs/bws/GSE57876_dCP_STARR_COMBINED.bw',
            '../tfs/bws/GSE57876_hkCP_STARR_COMBINED.bw'
            ]
    names = ['TRL2', 'CP190', 'PSER5_2', 'CTCF', 'ZELDA', "BEAF32", "DREF", "PC_2", "PH", "RAD21", "27me3", "H3K27ac", 'dcp', 'hkp']
    for ax, path, name in zip([ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, 
                               ax11, ax12, ax13, ax14, ax15, ax16], paths, names):
        bw = pyBigWig.open(path)
        xs = [start, end]
        ys = np.asarray(bw.values(chrom, xs[0], xs[1]))
        ymax = np.max(ys)
        ax.plot(range(*xs), ys, label='TRL')
        ax.set_xlim(xs[0], xs[-1])
        ax.set_xticks([])
        ax.set_yticks([int(ymax)])
        ax.set_ylabel(name, fontsize=6, rotation=0, labelpad=15)
        ax.tick_params(axis='both', which='major', labelsize=6)
    for ax in [ax17]:
        ax.axis('off')

    try:
        cmd = f"pyGenomeTracks --tracks ../loops/test.ini --region {chrom}:{start}-{end}  --width 38 --trackLabelFraction=0 --dpi 130 -o gtf_{chrom}_{start}_{end}.png"
        p = subprocess.call(cmd, shell=True)        
        img = mpimg.imread(f'gtf_{chrom}_{start}_{end}.png')
        ax2.imshow(img, aspect='auto')
    except:
        pass
    ax3.set_xticks([])
    
    ax1.pcolormesh(mat, cmap=cm.gist_heat_r, norm = colors.LogNorm())


    if show_tads == True:
        if tads.get(chrom):
            ymaxes = []
            ymins = []
            xs = []
            vs = []
            bob = 0
            for tad in tads.get(chrom):
                ss = tad[0]
                es = tad[1]
                if int(ss) > start and int(ss) < end:
                    bob += 1
                    if bob == 1:
                        print(f"tad: {tad}, ss: {ss}, es: {es}")
                    vs.append(tad[2])                    
                    x = (int(ss)-start)//res
                    y = mat.shape[1] - x
                    xs.append(x)
                    ymins.append(y-10)
                    ymaxes.append(y+10)
            if tad_colors is None:
                cs = []
                if len(vs) > 1:
                    if vs[0] == '+' or vs[0] == '-':
                        for v in vs:
                            if v=='+':
                                cs.append('teal')
                            else:
                                cs.append('black')
                        ax1.vlines(xs, ymins, ymaxes, colors=cs)
                else:
                    ax1.vlines(xs, ymins, ymaxes, colors='teal')
            else:
                cs = cmap(np.asarray(tad_colors)/vmax)
                ax1.vlines(xs, ymins, ymaxes, colors=cs) 

    for key in loops.keys():
        if key[0] != chrom:
            continue
        ss = key[1]
        if int(ss) > start:
            for partner in loops[key]:
                partner = partner
                if int(partner[1]) < end:
                    x = (int(ss)-start)//res
                    y = mat.shape[1]-(partner[1]-start)//res
                    ax1.scatter(x, y, marker='s', alpha = .4, c='gray', 
                                edgecolor='black', s=200)
    
    if show == False:
        plt.close()
    print(save)
    if save == True:
        print('saving to')
        fig.savefig(f"{save_dir}{chrom}_{start}_{end}.png")
    return fig, ax1, ax2, ax3



# Plotting a looplist
def plot_looplist(looplist, chrom, start, end, c, tads = {}, plotLoops = True, show = False, save = False, save_dir = '../standard_plots', show_tads = True, tad_colors = None, cmap = plt.cm.bwr, vmin = 0, vmax = 17500, axs = None):
    if start < 0 or end > len(chromos[chrom]):
        print("Those coordinates are not valid!")
        return False
    res = c.binsize
    mat = np.rot90(c.matrix(balance=True).fetch((chrom, start, end)))

    matmax = np.max(mat[~np.isnan(mat)])
    if axs == None:
        fig, ax = plt.subplots(figsize = (10, 10), dpi=200)

        ax2 = plt.subplot2grid((100,100), (0, 0), colspan = 90, rowspan = 2, fig = fig)
        ax3 = plt.subplot2grid((100,100), (3, 0), colspan = 90, rowspan = 2, fig = fig)
        ax4 = plt.subplot2grid((100,100), (6, 0), colspan = 90, rowspan = 2, fig = fig)
        ax5 = plt.subplot2grid((100,100), (9, 0), colspan = 90, rowspan = 2, fig = fig)
        ax6 = plt.subplot2grid((100,100), (12, 0), colspan = 90, rowspan = 2, fig = fig)
        ax7 = plt.subplot2grid((100,100), (15, 0), colspan = 90, rowspan = 2, fig = fig)
        ax8 = plt.subplot2grid((100,100), (18, 0), colspan = 90, rowspan = 2, fig = fig)
        ax9 = plt.subplot2grid((100,100), (21, 0), colspan = 90, rowspan = 2, fig = fig)
        ax10 = plt.subplot2grid((100,100), (25, 0), colspan = 90, rowspan = 2, fig = fig)
        ax11 = plt.subplot2grid((100,100), (28, 0), colspan = 90, rowspan = 2, fig = fig)
        ax12 = plt.subplot2grid((100,100), (31, 0), colspan = 90, rowspan = 2, fig = fig)
        ax13 = plt.subplot2grid((100,100), (34, 0), colspan = 90, rowspan = 2, fig = fig)
        ax14 = plt.subplot2grid((100,100), (37, 0), colspan = 90, rowspan = 2, fig = fig)
        ax15 = plt.subplot2grid((100,100), (40, 0), colspan = 90, rowspan = 2, fig = fig)
        ax16 = plt.subplot2grid((100,100), (43, 0), colspan = 90, rowspan = 2, fig = fig)
        ax17 = plt.subplot2grid((100,100), (46, 0), colspan = 90, rowspan = 2, fig = fig)

        ax1 = plt.subplot2grid((100,100), (50, 0), colspan = 90, rowspan = 50, fig = fig)

    else:
        fig = 1
        save = False
        [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12, ax13, ax14, ax15, ax16, ax17] = axs

    ax2.set_xticks([])
    ax2.set_yticks([])

    paths = ["../tfs/bws/GSE152770_GAF_GFP_COMBINED.bw", "../tfs/bws/CP190.bw", "../tfs/bws/GSE62925_pol2_COMBINED.bw", 
            '../tfs/bws/ctcf.bw', '../tfs/bws/zelda.bw', '../tfs/bws/BEAF32_SRR1636749.bw', '../tfs/bws/DREF_SRR1636771.bw', 
            '../tfs/bws/GSE60428_PC_COMBINED.bw', '../tfs/bws/GSE60428_PH_COMBINED.bw', '../tfs/bws/rad21_COMBINED.bw', 
            '../histones/bws/H3K27me3.bw', '../histones/bws/H3K27ac.bw', '../tfs/bws/GSE57876_dCP_STARR_COMBINED.bw',
            '../tfs/bws/GSE57876_hkCP_STARR_COMBINED.bw'
            ]
    names = ['TRL2', 'CP190', 'PSER5_2', 'CTCF', 'ZELDA', "BEAF32", "DREF", "PC_2", "PH", "RAD21", "27me3", "H3K27ac", 'dcp', 'hkp']
    for ax, path, name in zip([ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, 
                               ax11, ax12, ax13, ax14, ax15, ax16], paths, names):
        bw = pyBigWig.open(path)
        xs = [start, end]
        ys = np.asarray(bw.values(chrom, xs[0], xs[1]))
        ymax = np.max(ys)
        ax.plot(range(*xs), ys, label='TRL')
        ax.set_xlim(xs[0], xs[-1])
        ax.set_xticks([])
        ax.set_yticks([int(ymax)])
        ax.set_ylabel(name, fontsize=6, rotation=0, labelpad=15)
        ax.tick_params(axis='both', which='major', labelsize=6)
    for ax in [ax17]:
        ax.axis('off')

    try:
        cmd = f"pyGenomeTracks --tracks ../loops/test.ini --region {chrom}:{start}-{end}  --width 38 --trackLabelFraction=0 --dpi 130 -o gtf_{chrom}_{start}_{end}.png"
        p = subprocess.call(cmd, shell=True)        
        img = mpimg.imread(f'gtf_{chrom}_{start}_{end}.png')
        ax2.imshow(img, aspect='auto')
    except:
        pass
    ax3.set_xticks([])
    
    ax1.pcolormesh(mat, cmap=cm.gist_heat_r, norm = colors.LogNorm(vmin = 0.00005, vmax = matmax, ))


    if show_tads == True:
        if tads.get(chrom):
            ymaxes = []
            ymins = []
            xs = []
            for tad in tads.get(chrom):
                ss = tad[1]
                es = tad[2]
                if int(ss) > start and int(ss) < end:
                    x = (int(ss)-start)//res
                    y = mat.shape[1] - x
                    xs.append(x)
                    ymins.append(y-10)
                    ymaxes.append(y+10)
            if tad_colors is None:
                ax1.vlines(xs, ymins, ymaxes, colors='teal')
            else:
                cs = cmap(np.asarray(tad_colors)/vmax)
                ax1.vlines(xs, ymins, ymaxes, colors=cs) 


    markers = ['s','D','o','<']
    for i, loops in enumerate(looplist):
        xs, ys = [], []
        for key in loops.keys():
            if key[0] != chrom:
                continue
            ss = key[1]
            if int(ss) > start:
                for partner in loops[key]:
                    partner = partner
                    if int(partner[1]) < end:
                        x = (int(ss)-start)//res
                        y = mat.shape[1]-(partner[1]-start)//res
                        xs.append(x)
                        ys.append(y)
        ax1.scatter(xs, ys, marker=markers[i], alpha = .4,
                    edgecolor='black', s=200)

        
    if show == False:
        plt.close()
    print(save)
    if save == True:
        print('saving to')
        fig.savefig(f"{save_dir}{chrom}_{start}_{end}.png")
    return fig, ax1, ax2, ax3


