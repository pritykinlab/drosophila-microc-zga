import pandas as pd
import pybedtools as pbt
import numpy as np
from copy import deepcopy
# PROCESS DIFF_MUSTACHE OUTPUT (NOT DESEQ2)

def get_anchors(loops):
    anchors = {}
    dups = 0
    for key in loops.keys():
        if int(key[0][1]) > int(key[1][1]):
            continue
        val = loops[key]
        for i in [0, 1]:
            if anchors.get(key[i]):
                anchors[key[i]] += val
                dups += 1
            else:
                anchors[key[i]] = val
    return anchors

def dump_ancs(anchors, name, delim='\t'):
    with open(name, 'w') as f:
        for anc in anchors:
            chrom, s, e = anc
            val = anchors[anc]
            chrom, s, e, = map(str, [chrom, s, e])
            line = f"{chrom}{delim}{s}{delim}{e}\n"
            f.write(line)

def make_int(loop):
    return loop[0], int(loop[1]), int(loop[2])

import cooler
import scipy
import scipy.ndimage
# Processes the results of the mustache loops, where all loops are piled into 'out.tsv'
import networkx as nx
# Processes the results of the mustache loops, where all loops are piled into 'out.tsv'
import networkx as nx
def remove_adjacent_loops(pbt_loops, wsz):
    loopset = set()
    for i in pbt_loops:
        loop1 = tuple(list(make_int(tuple(i[:3]))) + list(make_int(tuple(i[3:6]))))
        loopset.add(loop1)

    loop_to_ind = {}
    ind_to_loop = {}
    for c, loop1 in enumerate(loopset):
        loop_to_ind[loop1] = c
        ind_to_loop[c] = loop1

    n = len(loopset)
    loop_adjacency = np.zeros((n, n))

    for loop1 in loopset:
        delta = loop1[2] - loop1[1]
        for i in range(-wsz, wsz):
            for j in range(-wsz, wsz):
                newloop = list(loop1)
                newloop[1] += delta*i
                newloop[2] += delta*i

                newloop[4] += delta*j
                newloop[5] += delta*j

                if tuple(newloop) in loopset:
                    c1 = loop_to_ind[loop1]
                    c2 = loop_to_ind[tuple(newloop)]
                    loop_adjacency[c1, c2] = 1
    
    G = nx.from_numpy_array(loop_adjacency)

    newloops = []
    for comp in nx.connected_components(G):
        starts = []
        ends = []
        for i in comp:
            l = ind_to_loop[i]
            chrom, s, e = l[0], l[1], l[4]
            starts.append(s)
            ends.append(e)

        start = int(np.mean(starts)//delta*delta)
        end = int(np.mean(ends)//delta*delta)
        newloop = (chrom, start, start+delta, chrom, end, end+delta)
        newloops.append(newloop)

    newloops = sorted(list(newloops))

    newancs = set()
    for l in sorted(list(newloops)):
        l1, l2 = l[:3], l[3:6]
        newancs.add(l1)
        newancs.add(l2)
    newancs = sorted(list(newancs))
    return newloops, newancs



def remove_adjacent_ancs(pbt_loops, wsz):

    ancset = set()
    for i in pbt_loops:
        anc1, anc2 = tuple(make_int(tuple(i[:3]))), tuple(make_int(tuple(i[3:6])))
        ancset.add(anc1)
        ancset.add(anc2)

    anc_to_ind = {}
    ind_to_anc = {}
    for c, anc1 in enumerate(ancset):
        anc_to_ind[anc1] = c
        ind_to_anc[c] = anc1

    n = len(ancset)
    anc_adjacency = np.zeros((n, n))

    for anc1 in ancset:
        delta = anc1[2] - anc1[1]
        for i in range(-wsz, wsz):
            newanc = list(anc1)
            newanc[1] += delta*i
            newanc[2] += delta*i

            if tuple(newanc) in ancset:
                c1 = anc_to_ind[anc1]
                c2 = anc_to_ind[tuple(newanc)]
                anc_adjacency[c1, c2] = 1
    
    G = nx.from_numpy_array(anc_adjacency)

    anc_mapping = {}
    for comp in nx.connected_components(G):
        starts = []
        for i in comp:
            anc = ind_to_anc[i]
            chrom, s = anc[0], anc[1]
            starts.append(s)

        start = int(np.mean(starts)//delta*delta)
        newanc = (chrom, start, start+delta)
        for i in comp:
            anc = ind_to_anc[i]
            anc_mapping[anc] = newanc

    newloops = []
    for i in pbt_loops:
        anc1, anc2 = tuple(make_int(tuple(i[:3]))),  tuple(make_int(tuple(i[3:6])))
        newanc1, newanc2 = anc_mapping[anc1], anc_mapping[anc2]
        newloops.append(list(newanc1) + list(newanc2)) 
    newloops = sorted(list(newloops))

    newancs = set()
    for l in newloops:
        l1, l2 = tuple(l[:3]), tuple(l[3:6])
        newancs.add(l1)
        newancs.add(l2)
    newancs = sorted(list(newancs))
    return newloops, newancs

def make_int2(loop):
    return (loop[0], int(loop[1]), int(loop[2]), loop[3], int(loop[4]), int(loop[5]))

import networkx as nx
# Processes the results of the mustache loops, where all loops are piled into 'out.tsv'
def process_loops(directory, file, cool, expected, res, merge_wsz, center_wsz):
    df = pd.read_csv(directory + file, sep='\t').iloc[:, :6].astype(object)

    df = df[~(df.iloc[:, 0] == 'BIN1_CHR')]
    df['BIN1_START'] = df['BIN1_START'].astype(int)
    df['BIN1_END'] = df['BIN1_END'].astype(int)
    df['BIN2_START'] = df['BIN2_START'].astype(int)
    df['BIN2_END'] = df['BIN2_END'].astype(int)
    df.to_csv(directory + 'all_loops_tmp.loops', sep='\t', header=False, index=None)

    pbt_loops = pbt.BedTool(directory + 'all_loops_tmp.loops')

    newloops, newancs = remove_adjacent_ancs(pbt_loops, wsz=merge_wsz)
    centered_loops = center_loops(newloops, cool, expected, res=res, wsz=center_wsz)

    newancs = set()
    for l in sorted(list(centered_loops)):
        l1, l2 = l[:3], l[3:6]
        newancs.add(l1)
        newancs.add(l2)
    newancs = sorted(list(newancs))

    not_nan_loops = []
    for loop in centered_loops:
        chrom, s, e = loop[0], loop[1], loop[4]
        try:
            v = cool.matrix(balance=True).fetch((chrom, s-5*res, s+5*res), (chrom, e-5*res, e+5*res))
        except Exception as e:
            if 'Genomic region out of bounds' in e.args[0]:
                pass
            else:
                raise Exception
        if np.isnan(v).mean() > .2:
            pass
        else:
            not_nan_loops.append(loop)
    

    not_nan_loops = pbt.BedTool([make_int2(x) for x in not_nan_loops])
    final_loops, final_ancs = remove_adjacent_loops(not_nan_loops, center_wsz)

    return final_ancs, final_loops


# Processes the results of the mustache loops, where all loops are piled into 'out.tsv'
def process_loops_without_moving(directory, extend_bins, res):
    df = pd.read_csv(directory + '/out.tsv', sep='\t').iloc[:, :6].astype(object)

    df = df[~(df.iloc[:, 0] == 'BIN1_CHR')]
    df['BIN1_START'] = df['BIN1_START'].astype(int)
    df['BIN1_END'] = df['BIN1_END'].astype(int)
    df['BIN2_START'] = df['BIN2_START'].astype(int)
    df['BIN2_END'] = df['BIN2_END'].astype(int)

    df['BIN1_START'] = df['BIN1_START'] - extend_bins*res
    df['BIN1_END'] = df['BIN1_END'] + extend_bins*res
    df['BIN2_START'] = df['BIN2_START'] - extend_bins*res
    df['BIN2_END'] = df['BIN2_END'] + extend_bins*res

    df.to_csv(directory + 'full_loops.loops', sep='\t', header=False, index=None)
 
    full_loops = pbt.BedTool(directory + 'full_loops.loops')
    loops = {}
    for loop in full_loops:
        l1 = tuple(loop[:3])
        l2 = tuple(loop[3:6])
        loops[(l1, l2)] = 1
    dump_ancs(get_anchors(loops), directory + 'full_ancs.csv')
    anchors = pbt.BedTool(directory + 'full_ancs.csv')
    return anchors, full_loops



def make_loopset(loops1, res):
    loopset = set()
    for l in loops1:
        l1, l2 = l[:3], l[3:6]
        chrom, s, e = l1[0], l1[1], l2[1]
        s, e = map(int, [s, e])
        s = s//res*res
        e = e//res*res
        newl1 = (chrom, s, s+res)
        newl2 = (chrom, e, e+res)
        loopset.add((newl1, newl2))
    return loopset
        
def check_if_in_loops(loop, loopset, res, wsz=5):
    l1, l2 = loop
    chrom, s, e = l1[0], l1[1], l2[1]
    s, e = map(int, [s, e])
    for i in range(-wsz, wsz):
        for j in range(-wsz, wsz+1):
            newl1 = (chrom, s-i*res, s-i*res+res)
            newl2 = (chrom, e-j*res, e-j*res+res)
            if (newl1, newl2) in loopset:
                return 1
    return 0

def compare_loops(loops1, loops2, res):
    loops1_set = make_loopset(loops1, res)
    loops2_set = make_loopset(loops2, res)
    c = 0 
    not_shared_loops = []
    shared_loops = []
    for loop in loops1_set:
        v = check_if_in_loops(loop, loops2_set, res)
        if v == 0:
            not_shared_loops.append(loop)
        elif v == 1:
            shared_loops.append(loop)
        c = c + v        
    return c, (not_shared_loops, shared_loops)


import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors


import scipy
def center_loops(loops, cool, expected_df, res, wsz=20):
    roll_sz = 6
    mat_dict = {}
    expecteds = {}
    pc_rat = .5
    n=2*wsz+1

    newloops = []    
    for c, i in enumerate(loops):
        l1, l2 = i[:3], i[3:6]
        chrom = l1[0]
        s = ((int(l1[1]) + int(l1[2]))//2)//res*res
        e = ((int(l2[1]) + int(l2[2]))//2)//res*res

        if (e-s)/res < wsz:
            continue
        expected = make_expected(expected_df, chrom, s, e, wsz, res=res, balanced=True)
        new_l1 = (l1[0], int(s)-res*wsz, int(s)+res*(wsz+1))

        new_l2 = (l2[0], int(e)-res*wsz, int(e)+res*(wsz+1))

        try:
            val = np.asarray(cool.matrix(balance=True).fetch(new_l1, new_l2))
        except Exception as e:
            if 'Genomic region out of bounds' in e.args[0]:
                continue
            else:
                raise Exception

        if val.shape != (2*wsz+1, 2*wsz+1):
            continue
        obs_div_exp = np.log2((val + expected*pc_rat)/(expected + expected*pc_rat))
        midpoint = len(val)//2
        smooth_mat_input = deepcopy(obs_div_exp)
        smooth_mat_input[np.isnan(smooth_mat_input)] = 0
        smooth_mat = scipy.ndimage.gaussian_filter(smooth_mat_input, sigma=1)
        smooth_mat[np.isnan(obs_div_exp)] = -10
        v = np.asarray(np.unravel_index(np.nanargmax(smooth_mat[midpoint-roll_sz:midpoint+roll_sz+1, midpoint-roll_sz:midpoint+roll_sz+1]), (roll_sz*2+1, roll_sz*2+1)))
        left, right = roll_sz-v[0], roll_sz-v[1]


        newmat = np.roll(val, left, axis=0)
        newmat = np.roll(newmat, right, axis=1)


        fulloop = tuple(i[:6])
        newstart = int(i[1]) - res*left
        newend = int(i[4]) - res*right


        new_l1 = (l1[0], newstart, newstart+res)
        new_l2 = (l2[0], newend, newend+res)
        new_fulloop = tuple(list(new_l1) + list(new_l2))
        newloops.append(new_fulloop)


        # val2 = np.asarray(cool.matrix(balance=True).fetch(new_l1, new_l2))
        # fig, axs = plt.subplots(1, 4, figsize=(15, 5))
        # axs[1].matshow(val, cmap=cm.gist_heat_r)
        # axs[0].matshow(smooth_mat, cmap=cm.bwr, vmin=-4, vmax=4)
        # axs[2].matshow(newmat, cmap=cm.gist_heat_r)
        # axs[3].matshow(val2, cmap=cm.gist_heat_r)


    return newloops


def make_expected(df, chrom, s, e, wsz, res=400, balanced=True):
    chrom = str(chrom)
    diag = (e-s)//res
    expected = np.zeros((wsz*2+1, wsz*2+1))
    n = expected.shape[0]

    test = np.zeros((wsz*2+1, wsz*2+1))
    n = test.shape[0]

    averages = df[(df.region == chrom)]

    bottom = diag
    
    diags = np.flip(np.indices((n, n))[0]) + (np.indices((n, n))[1]) - (n-1)
    assert (np.diag(diags) == 0).all()
    diags = diags + diag
    for val in np.unique(diags):
        if balanced==True:
            row = averages[averages['diag'] == val]
            if len(row) < 1:
                expected[diags==val] = np.nan
            else:
                expected[diags==val] = float(row['balanced.avg'])
        elif balanced==False:
            row = averages[averages['diag'] == val]
            if len(row) < 1:
                expected[diags==val] = np.nan
            else:
                expected[diags==val] = float(row['count.avg'])
    return expected


import os
import statsmodels
import statsmodels.stats
import statsmodels.stats.multitest
def process_deseq_output(deseq_output_file):
    deseq_results = pd.read_csv(deseq_output_file, sep=' ')

    shifted = [int(x) for x in list(zip(*pd.Series(deseq_results.index).apply(lambda x: x.split("_")).values))[0]]
    new_frame = pd.DataFrame()
    for i, column in enumerate(['BIN1_CHR', 'BIN1_START', 'BIN1_END', 'BIN2_CHR', 'BIN2_START', 'BIN2_END']):
        new_frame[column] = list(zip(*pd.Series(deseq_results.index).apply(lambda x: x.split("_")).values))[i+1]
        if 'START' in column or 'END' in column:
            new_frame[column] = np.asarray(list(zip(*pd.Series(deseq_results.index).apply(lambda x: x.split("_")).values))[i+1]).astype(int)
    new_frame['TYPE'] = (deseq_results['log2FoldChange'] > 0).values.astype(float)
    new_frame['padj'] = deseq_results['padj'].values.astype(float)
    new_frame['pval'] = deseq_results['pvalue'].values.astype(float)
    new_frame['l2fc'] = deseq_results['log2FoldChange'].values.astype(float)
    new_frame['baseMean'] = deseq_results['baseMean'].values.astype(float)
    new_frame['shifted'] = shifted
    new_frame = new_frame[new_frame['shifted'] == 0]
    new_frame['padj'] = statsmodels.stats.multitest.fdrcorrection(new_frame['pval'])[1]
    return new_frame

    # new_frame.to_csv(directory + 'all_loops_tmp.loops', sep='\t', header=None, index=None)
    # new_frame[(new_frame['padj'] < .05)].to_csv(directory + 'deseq_loops_tmp.loops', sep='\t', header=None, index=None)
    # new_frame[(new_frame['padj'] > .05)].to_csv(directory + 'non_significant_tmp.loops', sep='\t', header=None, index=None)
    # new_frame = new_frame[new_frame['padj'] < .05]
        
        
        
        
    # loops_diff = pd.read_csv(directory + 'all_loops_tmp.loops', delimiter='\t', 
    #                      names=['BIN1_CHR', 'BIN1_START', 'BIN1_END', 'BIN2_CHR', 'BIN2_START', 'BIN2_END', 'TYPE' ,'padj', 'l2fc', 'shifted'])
    # loops_diff['BIN1_CHR'] = loops_diff['BIN1_CHR'].apply(str) 
    # loops_diff['BIN2_CHR'] = loops_diff['BIN2_CHR'].apply(str) 
    # loops_diff['FILLER'] = 1
    # loops_diff['STRAND'] = '*'
    
        

    # newloops['l2fc'] = -newloops['l2fc']
    # significant_loops = newloops[newloops['padj'] <= .05]
    # nonsignificant_loops = newloops[~(newloops['padj'] <= .05)]

    # treg_loops = significant_loops[significant_loops['l2fc'] < 0]
    # tcon_loops = significant_loops[significant_loops['l2fc'] > 0]

    # newloops.to_csv(directory + 'full_loops.loops', sep='\t', header=None, index=None)

    # significant_loops.to_csv(directory + 'significant_loops.loops', sep='\t', header=None, index=None)
    # nonsignificant_loops.to_csv(directory + 'nonsignificant_loops.loops', sep='\t', header=None, index=None)
    # treg_loops.to_csv(directory + 'significant_treg_loops.loops', sep='\t', header=None, index=None)
    # tcon_loops.to_csv(directory + 'significant_tcon_loops.loops', sep='\t', header=None, index=None)
    
    # full_loops = pbt.BedTool(directory + 'full_loops.loops')
    # peaks_loops = pbt.BedTool(directory + 'significant_loops.loops')
    # nonsig_loops = pbt.BedTool(directory + 'nonsignificant_loops.loops')
    # treg_loops = pbt.BedTool(directory + 'significant_treg_loops.loops')
    # tcon_loops = pbt.BedTool(directory + 'significant_tcon_loops.loops')

    # loops = {}
    # for loop in full_loops:
    #     l1 = tuple(loop[:3])
    #     l2 = tuple(loop[3:6])
    #     loops[(l1, l2)] = 1
    # dump_ancs(get_anchors(loops), directory + 'full_ancs.csv')
    
    
    # loops = {}
    # for loop in nonsig_loops:
    #     l1 = tuple(loop[:3])
    #     l2 = tuple(loop[3:6])
    #     loops[(l1, l2)] = 1
    # dump_ancs(get_anchors(loops), directory + 'nonsignificant_ancs.csv')
    
    
    # loops = {}
    # for loop in peaks_loops:
    #     l1 = tuple(loop[:3])
    #     l2 = tuple(loop[3:6])
    #     loops[(l1, l2)] = 1
    # dump_ancs(get_anchors(loops), directory + 'significant_ancs.csv')

    # loops = {}
    # for loop in treg_loops:
    #     l1 = tuple(loop[:3])
    #     l2 = tuple(loop[3:6])
    #     loops[(l1, l2)] = 1
    # dump_ancs(get_anchors(loops), directory + 'treg_ancs.csv')

    # loops = {}
    # for loop in tcon_loops:
    #     l1 = tuple(loop[:3])
    #     l2 = tuple(loop[3:6])
    #     loops[(l1, l2)] = 1
    # dump_ancs(get_anchors(loops), directory + 'tcon_ancs.csv')

    
    # full_anchors = pbt.BedTool(directory + 'full_ancs.csv')
    # peaks_anchors = pbt.BedTool(directory + 'significant_ancs.csv')
    # treg_anchors = pbt.BedTool(directory + 'treg_ancs.csv').sort()
    # tcon_anchors = pbt.BedTool(directory + 'tcon_ancs.csv').sort()
    # nonsig_anchors = pbt.BedTool(directory + 'nonsignificant_ancs.csv')
    
    # merged_peak_anchors = peaks_anchors.sort()

    # loops_diff = pd.read_csv(directory + 'significant_loops.loops', delimiter='\t', 
    #                      names=['BIN1_CHR', 'BIN1_START', 'BIN1_END', 'BIN2_CHR', 'BIN2_START', 'BIN2_END', 'TYPE' ,'padj', 'shifted', 'FILLER', 'STRAND', 'l2fc'])
    # loops_diff['BIN1_CHR'] = loops_diff['BIN1_CHR'].apply(str) 
    # loops_diff['BIN2_CHR'] = loops_diff['BIN2_CHR'].apply(str)   

    # return (full_loops, full_anchors, peaks_loops, treg_loops, tcon_loops, nonsig_loops,
    #         merged_peak_anchors, treg_anchors, tcon_anchors, nonsig_anchors, loops_diff)



