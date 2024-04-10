from copy import deepcopy
from itertools import product
import numpy as np
import pandas as pd 
import pyBigWig
import scipy.stats
import cooler
import scipy
import random
import glob

chromo_keys = ['2L', '2R', '3L', '3R', '4', 'M', 'X', 'Y']
chromos = {}
for chromo in chromo_keys:
    with open(f'./dm6/chr{chromo}.txt', 'r') as f:
        chromos[chromo] = f.readline()

import time
import networkx as nx


def dump_ancs(loopset, name, delim='\t'):
    with open(name, 'w') as f:
        for anc in loopset.anchors:
            chrom, s, e = anc
            s, e, = map(str, [s, e])
            line = f"{chrom}{delim}{s}{delim}{e}\n"
            f.write(line)


def deposit(otherloops, name, delim='\t'):
    with open(f'{name}.txt', 'w') as f:
        for loop in otherloops.keys():
            loop1, loop2 = loop
            loop1 = [str(l) for l in loop1]
            loop2 = [str(l) for l in loop2]
            line = f"{delim}".join(loop1) + f"{delim}" + f"{delim}".join(loop2) + "\n"
            f.write(line)

def get_intersection(dic1, dic2):
    res = {}
    for key1 in dic1.keys():
        v1 = dic1.get(key1)
        if v1:
            v2 = dic2.get(key1)
            if v2:
                res[key1] = (v1, v2)
    return res

def get_union(dic1, dic2):
    res = {}
    for key1 in dic1.keys():
        v1 = dic1.get(key1)
        if v1:
            res[key1] = (v1)

    for key2 in dic2.keys():
        v2 = dic2.get(key2)
        if v2:
            res[key2] = (v2)
    return res    

def get_ss(loop):
    if type(loop) != list:
        start = loop.split("|")[1]
    else:
        start = loop[1]
    return start


def write_ancs(anchors, name, buffer, res):
    print(f"Writing to ../meme/for_meme/{name}.txt")
    print(f"Writing with buffer size {buffer}")
    with open(f"../meme/for_meme/{name}.txt", 'w') as anc_file:
        for anc in anchors:
            chrom, start, end = anc[:3]
            seq = get_sequence(chrom, start, end, res=res, buffer = buffer)
            write_file(anc, seq, anc_file)


def write_file(coords, seq, loop_file):
    coords = [str(c) for c in coords]
    loop_file.write(f">{'|'.join(coords)}\n")
    loop_file.write(f'{seq}\n')

def get_sequence(chrm, x1, x2, res=400, buffer = 0 ):
    x1 = int(x1)-res*buffer
    x2 = int(x2)+res*buffer
    return chromos[chrm][x1:x2]

def write_loops(loops, name, buffer):
    print(f"Writing to ../tfs/{name}.txt")
    print(f"Writing with buffer size {buffer}")
    written = {}
    with open(f"../tfs/{name}.txt", 'w') as loop_file:
        for key in loops.keys():
            for partner in loops[key]:
                if written.get(partner):
                    pass
                else:
                    seq = get_sequence(*partner, buffer = buffer)
                    write_file(partner, seq, loop_file)
                written[partner] = 1
            if written.get(key):
                pass
            else:
                seq = get_sequence(*key, buffer = buffer)
                write_file(key, seq, loop_file)
            written[key] = 1


# Function for making other loops
def make_otherloops(loops, inp = '../loops/v2_400kbLoops.txt', out = '../tfs/trl_db.txt', skip_header = False, val_index = -1):
    written = {}
    params = {}
    print(f'skip_header = {skip_header}')
    with open(out, 'w+') as loop_file:
        with open(inp, 'r') as f:
            for i, line in enumerate(f.readlines()):
                if skip_header == True and i==0:
                    continue
                vals = line.strip().split('\t')
                coordsx, coordsy = vals[:3], vals[3:6]                
                coordsx[1] = int(coordsx[1])
                coordsx[2] = int(coordsx[2])

                coordsy[1] = int(coordsy[1])
                coordsy[2] = int(coordsy[2])

                x_ss = get_ss(coordsx)
                y_ss = get_ss(coordsy)
                
                
                to_sort = sorted([(x_ss, coordsx), (y_ss, coordsy)])
                coordsx, coordsy = to_sort[0][1], to_sort[1][1]
                if coordsx == coordsy:
                    print(1)
                    continue
                assert(coordsx[1] < coordsy[1])

                seqx = get_sequence(*coordsx)
                seqy = get_sequence(*coordsy)
                if not written.get(tuple(coordsx)): # remove duplicates
                    write_file(coordsx, seqx, loop_file)
                if not written.get(tuple(coordsy)): # remove duplicates
                    write_file(coordsy, seqy, loop_file)
                written[tuple(coordsx)] = 1
                written[tuple(coordsy)] = 1

                key1 = tuple(coordsx)
                key2 = tuple(coordsy)

                val = float(vals[val_index])+.1

                loops[(key1, key2)] = val
    return params

# Getting tads
def get_tads(filename = '../tads/tads_domains.bed'):
    tads = {}
    with open(filename, 'r') as f:
        for line in f.readlines():
            vals = line.split('\t')
            chrom = vals[0]
            start, end = int(vals[1]), int(vals[2])
            score = float(vals[4])
            tads.setdefault(chrom, [])
            tads[chrom].append((start, end, score))
    return tads

# Function for making loops
def make_loops(loops, inp = '../loops/v2_400kbLoops.txt', out = '../tfs/trl_db.txt' ):
    written = {}
    params = {}
    with open(out, 'w+') as loop_file:
        with open(inp, 'r') as f:
            for i, line in enumerate(f.readlines()):
                if i == 0: # Header
                    continue
                vals = line.strip().split('\t')
                coordsx, coordsy = vals[:3], vals[3:6]                
                coordsx[1] = int(coordsx[1])
                coordsx[2] = int(coordsx[2])

                coordsy[1] = int(coordsy[1])
                coordsy[2] = int(coordsy[2])

                x_ss = get_ss(coordsx)
                y_ss = get_ss(coordsy)
                
                
                to_sort = sorted([(x_ss, coordsx), (y_ss, coordsy)])
                coordsx, coordsy = to_sort[0][1], to_sort[1][1]
                coordsx[0] = coordsx[0].replace('arm_', 'chr')
                coordsy[0] = coordsy[0].replace('arm_', 'chr')                
                if coordsx == coordsy:
                    continue
                assert(coordsx[1] < coordsy[1])

                seqx = get_sequence(*coordsx)
                seqy = get_sequence(*coordsy)
                if not written.get(tuple(coordsx)): # remove duplicates
                    write_file(coordsx, seqx, loop_file)
                if not written.get(tuple(coordsy)): # remove duplicates
                    write_file(coordsy, seqy, loop_file)
                written[tuple(coordsx)] = 1
                written[tuple(coordsy)] = 1


                key1 = tuple(coordsx)
                key2 = tuple(coordsy)

                if key1 == ('chrX', 1026400, 1027200):
                    print(coordsx)
                    print(coordsy)

                loops.setdefault(key1, [])
                loops[key1].append(key2)
                params[(key1, key2)] = vals[6:]
    return params


# def bw_pileup_anchors(anchors, wsz, loopset, binsize = 800,):
#     res = anchors[0][2] - anchors[0][1]
#     shape = (2*wsz+1, 2*wsz+1)
#     length = (2*wsz+1)*res//binsize
#     print(res, shape, length)
#     bws_mean = {'trl':np.zeros(length), 'ctcf':np.zeros(length), 'pser' : np.zeros(length), 'cp190': np.zeros(length), 'zelda': np.zeros(length)}
#     bws = {'trl':[], 'ctcf':[], 'pser' : [], 'cp190':[], 'zelda': []}

#     c = 0
#     for anchor in anchors:
#         l1 = anchor
#         chrom, loop1_start, loop1_end = l1
#         res = loop1_end-loop1_start
#         start1 = loop1_start-res*wsz
#         end1 = loop1_end+res*wsz
        
#         for name, bw in zip(('trl', 'ctcf', 'pser', 'cp190', 'zelda'), (loopset.pybw_trl, loopset.pybw_ctcf, loopset.pybw_pser5,
#                                                                     loopset.pybw_cp190, loopset.pybw_zelda)):
#             v1 = np.asarray(bw.values(chrom, start1, end1))
#             v1 = np.reshape(v1, (length, binsize))
#             v1 = np.sum(v1, axis=1)
#             bws_mean[name] += v1
#             bws[name].append(v1)
#         c+=1
#     for name in ('trl', 'ctcf', 'pser', 'cp190', 'zelda'):
#         bws_mean[name] /= c
#     return bws_mean, bws


# def bw_pileup(otherloops, wsz, loopset, binsize = 800,):
#     res = list(otherloops.keys())[0][0][2] - list(otherloops.keys())[0][0][1]
#     shape = (2*wsz+1, 2*wsz+1)
#     length = (2*wsz+1)*res//binsize
#     print(res, shape, length)
#     L_bws_mean = {'trl':np.zeros(length), 'ctcf':np.zeros(length), 'pser' : np.zeros(length), 'cp190': np.zeros(length), 'zelda': np.zeros(length)}
#     R_bws_mean = {'trl':np.zeros(length), 'ctcf':np.zeros(length), 'pser' : np.zeros(length), 'cp190': np.zeros(length), 'zelda': np.zeros(length)}
#     L_bws = {'trl':[], 'ctcf':[], 'pser' : [], 'cp190':[], 'zelda': []}
#     R_bws = {'trl':[], 'ctcf':[], 'pser' : [], 'cp190':[], 'zelda': []}


#     c = 0
#     for loop in otherloops.keys():
#         l1, l2 = loop
#         chrom, loop1_start, loop1_end = l1
#         res = loop1_end-loop1_start
#         start1 = loop1_start-res*wsz
#         end1 = loop1_end+res*wsz

#         chrom, loop2_start, loop2_end = l2
#         res = loop2_end-loop2_start
#         start2 = loop2_start-res*wsz
#         end2 = loop2_end+res*wsz
        
#         for name, bw in zip(('trl', 'ctcf', 'pser', 'cp190', 'zelda'), (loopset.pybw_trl, loopset.pybw_ctcf, loopset.pybw_pser5,
#                                                                     loopset.pybw_cp190, loopset.pybw_zelda)):
#             v1 = np.asarray(bw.values(chrom, start1, end1))
#             v2 = np.asarray(bw.values(chrom, start2, end2))
#             v1 = np.reshape(v1, (length, binsize))
#             v1 = np.sum(v1, axis=1)
#             v2 = np.reshape(v2, (length, binsize))
#             v2 = np.sum(v2, axis=1)
#             L_bws_mean[name] += v1
#             R_bws_mean[name] += v2

#             L_bws[name].append(v1)
#             R_bws[name].append(v2)

#         c+=1
#     for name in ('trl', 'ctcf', 'pser', 'cp190', 'zelda'):
#         L_bws_mean[name] /= c
#         R_bws_mean[name] /= c
#     return L_bws_mean, R_bws_mean, L_bws, R_bws





def quant_anchor_df(anchors, bigwigs, res, buffer):
    anchor_df = pd.DataFrame()
    window=res*buffer
    print("window is", window)

    anclist = []
    extended_chrom, extended_start, extended_end, = [], [], []  
    for i in anchors:
        chrom, s, e = i[:3]
        s, e = map(int, [s, e])
        anclist.append((chrom, s, e))
        extended_chrom.append(chrom), extended_start.append(s-window), extended_end.append(e+window)

    for i, bigwig_name in enumerate(list(bigwigs)):
        vals = []
        bbi = bigwigs[bigwig_name] 
        #vs = bbi.stackup(extended_chrom, extended_start, extended_end)
        #vs[vs < 0] = 0
        #vs = np.log2((bbi.stackup(extended_chrom, extended_start, extended_end)+.1)/(bbi.info['summary']['mean']+.1))
        vs = np.log2(np.nanmean(bbi.stackup(extended_chrom, extended_start, extended_end)+.1, axis=1))
        #vs[vs < .11] = np.nan
        vals = vs
        anchor_df[bigwig_name] = vals
    anchor_df['anchor'] = anclist
    return anchor_df

def make_anchor_df(loopset, graph, ind2anc):
    anc_list = list(loopset.anchors.keys())
    anchor_df = pd.DataFrame()
    
    # columns = ['trl', 'ctcf', 'pser', 'cp190', 'zelda', 'housekeeping', ]
    # peak_dicts = [loopset.peak_trl, loopset.peak_ctcf, loopset.peak_pser, loopset.peak_cp190, loopset.peak_zelda, loopset.peak_housekeeping]

    columns = []
    peak_dicts = []
    
    for prot in loopset.peak_names:
        columns.append(prot)
        pybw = getattr(loopset, f'{prot}')
        peak_dicts.append(pybw)
    
    for col, peak_dict in zip(columns, peak_dicts):
        vals = []
        for i, anc in enumerate(anc_list):  
            if col =='housekeeping':
                v = peak_dict.get(anc, 0)
            else:
                v = peak_dict.get(anc, 0)  > 0
            vals.append(v)
        anchor_df[col] = vals

    anchor_df['anchor'] = anc_list
    degs = np.zeros(len(anchor_df))

    for node in graph.nodes:
        anc = ind2anc[node]
        deg = graph.degree[node]
        ind = np.where(anchor_df.anchor == anc) 
        degs[ind] = deg
    anchor_df['degs'] = degs
    # enhancer_vals =  np.asarray([loopset.peak_enhancer.get(anc) != None for anc in anc_list])
    # anchor_df['enhancer'] = enhancer_vals  
    return anchor_df


def annotation_to_dict(annotation, rowname='is_housekeeping'):
    dic = {}
    for i, row in annotation.iterrows():
        name = row['gene_name']
        val = row[rowname]
        if val==1:
            dic[name] = 1
    return dic

# Used to get housekeeping genes
def get_bed(filename= '../tfs/all_crms.bed'):
    chrom_macs_trl = {}
    with open(filename, 'r') as f:
        for i, line in enumerate(f.readlines()):
            values = line.strip().split('\t')
            chrom, start, end = values[2], int(values[4]), int(values[5])
            name = values[12]
            chrom_macs_trl.setdefault(chrom, [])
            chrom_macs_trl[chrom].append((start, end, name))
    for chrom in chrom_macs_trl.keys():
        chrom_macs_trl[chrom] = sorted(chrom_macs_trl[chrom])
    return chrom_macs_trl    

def get_TSS(filename= '../tfs/all_crms.bed'):
    chrom_macs_trl = {}
    with open(filename, 'r') as f:
        for i, line in enumerate(f.readlines()):
            values = line.strip().split('\t')
            chrom, start, end = values[0], int(values[1]), int(values[2])
            name = values[3]
            reads = values[4]
            orientation = values[5]
            chrom_macs_trl.setdefault(chrom, [])
            chrom_macs_trl[chrom].append((start, end, orientation))
    for chrom in chrom_macs_trl.keys():
        chrom_macs_trl[chrom] = sorted(chrom_macs_trl[chrom])
    return chrom_macs_trl    


def make_values_array(loopset, wsz,):
    print(f"Making values with wsz = {wsz}")
    anc_list = list(loopset.anchors)
    og_values = np.zeros((len(anc_list), 5))

    for i, anc in enumerate(anc_list):
        trl = get_anc_values([anc], loopset.pybw_trl, wsz)[0]
        ctcf = get_anc_values([anc], loopset.pybw_ctcf, wsz)[0]
        pser = get_anc_values([anc], loopset.pybw_pser5, wsz)[0]
        cp190 = get_anc_values([anc], loopset.pybw_cp190, wsz)[0]
        zelda = get_anc_values([anc], loopset.pybw_zelda, wsz)[0]

        og_values[i, 0] = trl
        og_values[i, 1] = ctcf
        og_values[i, 2] = pser
        og_values[i, 3] = cp190
        og_values[i, 4] = zelda

    # values = np.log(og_values)
    values = (og_values)
    where = np.where(np.isinf(values))
    for (i, j) in zip(where[0], where[1]):
        newvals = values[:, j]
        values[i, j] = np.min(newvals[~np.isinf(newvals)])

    return values

def make_loop_df(loopset, anchor_df, anc2ind):
    loop_df = pd.DataFrame()
    otherloops = loopset.other_loops
    loop_values = np.empty(shape=(len(otherloops), len(anchor_df.columns)*2), dtype=object)
    for i, loop in enumerate(otherloops.keys()):
        anc1, anc2 = loop
        ind1, ind2 = anc2ind[anc1], anc2ind[anc2]
        loop_values[i, :] = pd.concat((anchor_df.iloc[ind1, :], anchor_df.iloc[ind2, :])).to_numpy()

    n_cols = len(anchor_df.columns)
    for i, column in enumerate(anchor_df.columns):
        loop_df[column + "_1"] = (loop_values[:, i])
        loop_df[column + "_2"] = (loop_values[:, i+n_cols])

        if column != 'anchor':    
            loop_df[column + "_1"] = loop_df[column + "_1"].astype(float)
            loop_df[column + "_2"] = loop_df[column + "_2"].astype(float)
            
    dists = []
    for i, row in loop_df.iterrows():
        s1, s2 = row['anchor_1'][1], row['anchor_2'][1]
        assert s2 > s1
        dists.append(s2-s1)
    loop_df['dists'] = dists

    for name in anchor_df.columns:
        if name == 'anchor':
            continue
        else:
            loop_df[name] = loop_df[name + "_1"].values + loop_df[name + '_2'].values
    return loop_df

def df2loops(df, indices):
    otherloops = {}
    for ind in indices:
        anc1, anc2 = df.loc[ind, 'anchor_1'], df.loc[ind, 'anchor_2'],
        anc1, anc2 = sorted([anc1, anc2])
        assert(anc1[1] < anc2[1])
        otherloops[(anc1, anc2)] = 1
    return otherloops


def skip_tad_pileup(anchors, wsz, loopset, binsize=800,):
    shape = (2*wsz+1, 2*wsz+1)
    length = (2*wsz+1)
    print(shape, length)
    print("NOT: Dividing each pileup by its sum – interpreting loop as probability distribution")
    bws_mean = {'trl':np.zeros(length), 'ctcf':np.zeros(length), 'pser' : np.zeros(length), 'cp190': np.zeros(length), 'zelda': np.zeros(length)}

    c = 0
    mats = []
    meanval = np.zeros(shape)

    #### CAN WE PARALLELIZE THIS? ####
    for anchor in anchors:
        try:
            val = fetch_anchor_window(anchor, wsz, binsize)
        except:
            print("Error fetching anchor window")
            continue
        val[np.isnan(val)] = 0
        val /= np.sum(val)

        meanval += val
        # mats.append(val)
        chrom, anc_start, _ = anchor
        start = anc_start-binsize*wsz
        end = anc_start+binsize*(wsz+1)
        
        for name, bw in zip(('trl', 'ctcf', 'pser', 'cp190', 'zelda'), (loopset.pybw_trl, loopset.pybw_ctcf, loopset.pybw_pser5,
                                                                    loopset.pybw_cp190, loopset.pybw_zelda)):
            v1 = bw.values(chrom, start, end)
            v1 = np.reshape(v1, (length, binsize))
            v1 = np.sum(v1, axis=1)

            bws_mean[name] += v1

        c+=1
    for name in ('trl', 'ctcf', 'pser', 'cp190', 'zelda'):
        bws_mean[name] /= c
        bws_mean[name] /= c
    meanval /= c
    return meanval, bws_mean

def anchor_pileup(anchors, wsz, loopset, binsize=800,):
    shape = (2*wsz+1, 2*wsz+1)
    length = (2*wsz+1)
    print(shape, length)
    print("Dividing each pileup by its sum – interpreting loop as probability distribution")
    bws_mean = {'trl':np.zeros(length), 'ctcf':np.zeros(length), 'pser' : np.zeros(length), 'cp190': np.zeros(length), 'zelda': np.zeros(length)}

    c = 0
    mats = []
    meanvals = []
    #### CAN WE PARALLELIZE THIS? ####
    for anchor in anchors:
        try:
            val = fetch_anchor_window(anchor, wsz, binsize)
        except:
            print("Error fetching anchor window")
            continue
        val[np.isnan(val)] = 0

        meanvals.append(val)
        # mats.append(val)
        chrom, anc_start, _ = anchor
        start = anc_start-binsize*wsz
        end = anc_start+binsize*(wsz+1)
        
        for name, bw in zip(('trl', 'ctcf', 'pser', 'cp190', 'zelda'), (loopset.pybw_trl, loopset.pybw_ctcf, loopset.pybw_pser5,
                                                                    loopset.pybw_cp190, loopset.pybw_zelda)):
            v1 = bw.values(chrom, start, end)
            v1 = np.reshape(v1, (length, binsize))
            v1 = np.sum(v1, axis=1)

            bws_mean[name] += v1

        c+=1
    for name in ('trl', 'ctcf', 'pser', 'cp190', 'zelda'):
        bws_mean[name] /= c
        bws_mean[name] /= c
    return meanvals, bws_mean

def df_subset(prot_name, dfs, loopset):
    peak_df = pd.DataFrame([*dfs[prot_name].index], columns=['chrom', 'start', 'end']).astype(object)
    peak_df['signal'] = list(zip(peak_df.chrom, peak_df.start, peak_df.end))
    peaks = df_to_dict(peak_df)
    c = 0
    ancs = get_peak_dict(loopset.anchors, peaks, buffer = loopset.buffer, method='append')
    return ancs

def fetch_loop_window(loop, wsz, res):
    co = cooler.Cooler(f'../All_nc14.mcool::resolutions/{res}')
    coolmats = co.matrix(balance=True)
    chrom, start, end = loop[0], int(loop[1]), int(loop[4])
    s = res
    v = coolmats.fetch(f'{chrom}: {start-wsz*s}-{start+(wsz+1)*s}',f'{chrom}: {end-wsz*s}-{end+(wsz+1)*s}')
    v[np.isnan(v)] = 0
    return v

def loops_to_pileup(loops, wsz, res=400, binsize=200):
    vals = []
    close_to_diag = 0
    for i, loop in enumerate(loops):
        loop = loop[:6]
        if np.abs(int(loop[1]) - int(loop[4]))//res < 30:
            close_to_diag += 1
            continue
        val = fetch_loop_window(loop, wsz, binsize)
        vals.append(val)
    vals  = np.asarray(vals)
    return vals


import scipy

def normalize_loop_pileup(mat_list, zero_cutoff = 1):   
    print("Aligning matrices based on argmax, taking the mean, then normalizing by center 3 value")
    meanmats = []
    wsz = 8
    newmeanvals = []
    midpoint = len(mat_list[0])//2
    length = len(mat_list[0])
    for mat in mat_list:
        smooth_mat = scipy.ndimage.gaussian_filter(mat, sigma=1)
        v = np.asarray(np.unravel_index(np.argmax(smooth_mat[midpoint-wsz:midpoint+wsz+1, midpoint-wsz:midpoint+wsz+1]), (wsz*2+1, wsz*2+1)))
        left, right = wsz-v[0], wsz-v[1]
        newmat = np.roll(mat, left, axis=0)
        newmat = np.roll(newmat, right, axis=1)

        zeromat = np.zeros(mat_list[0].shape)
        zeromat[wsz:length-wsz, wsz:length-wsz] = newmat[wsz:length-wsz, wsz:length-wsz]
        newmat = zeromat
        if newmat[midpoint, midpoint] == 0:
            continue
        newmeanvals.append(newmat)
    newmeanvals = np.asarray(newmeanvals)

    maxes = (newmeanvals[:, midpoint, midpoint])

    newmeanvals = (newmeanvals.T/maxes).T

    # m2 = np.mean(newmeanvals, axis=0)
    # m2 /= np.sum(m2[20-2:20+2+1, 20-2:20+2+1])
    return newmeanvals







def df_to_ancvals(prot_name, labels, label_keys, dfs, loopset, method='mean'):
    print("Aligning anchors with pre-called tad boundaries when possible")
    vs = []
    tads = get_tads()
    tad_boundaries = {}
    for chrom in tads:
        tad_boundaries[chrom] = []
        for item in tads[chrom]:
            tad_boundaries[chrom].append((item[0]-800, item[0]+800, item[2]))        

    for keys in label_keys:
        mask = pd.Series(labels).isin(keys).values

        signal = np.ones(len(dfs[prot_name].loc[mask, :]))
        peak_df = pd.DataFrame([*dfs[prot_name].loc[mask, :].index], columns=['chrom', 'start', 'end']).astype(object)
        peak_df['signal'] = signal
        peaks = df_to_dict(peak_df)
            
        peak_ancs = get_peak_dict(loopset.anchors, peaks, res=400, buffer=3)
        
        newancs = {}
        print(len(peak_ancs))
        for anc in peak_ancs:
            chrom = anc[0]
            peak = peak_recursion(anc, tad_boundaries[chrom], 400, 0)
            if peak:
                peak = peak[0]
                newanc = chrom, peak[0]+800, peak[1]-400
                assert newanc[2]-newanc[1] == 400
                newancs[anc]=1
            else:
                newancs[anc]=1    

        print("Piling up")
        meanval, bws_mean = anchor_pileup(newancs, 30, loopset, binsize=400,)
        vs.append(np.asarray(meanval))
    return vs

import time
def cube(x):
    loop = (('chr3R', 21061600, 21062000), ('chr3R', 21131200, 21131600))
    return fetch_loop_window(loop, 5, 400)





def loop_pileup(otherloops, wsz, loopset, method, binsize=800, zero_cutoff = 1):
    shape = (2*wsz+1, 2*wsz+1)
    length = (2*wsz+1)
    print(shape, length)
    print("Using method: ", method)
    print("Dividing each pileup by its sum – interpreting loop as probability distribution")
    L_bws_mean = {'trl':np.zeros(length), 'ctcf':np.zeros(length), 'pser' : np.zeros(length), 'cp190': np.zeros(length), 'zelda': np.zeros(length)}
    R_bws_mean = {'trl':np.zeros(length), 'ctcf':np.zeros(length), 'pser' : np.zeros(length), 'cp190': np.zeros(length), 'zelda': np.zeros(length)}
    L_bws = {'trl':[], 'ctcf':[], 'pser' : [], 'cp190':[], 'zelda': []}
    R_bws = {'trl':[], 'ctcf':[], 'pser' : [], 'cp190':[], 'zelda': []}
    c = 0
    mats = []
    meanval = np.zeros(shape)

    #### CAN WE PARALLELIZE THIS? ####
    for loop in otherloops.keys():
        try:
            val = fetch_loop_window(loop, wsz, binsize)
        except:
            continue
        val[np.isnan(val)] = 0
        if (val == 0).sum()/np.prod(val.shape) > zero_cutoff:
            continue
        # val /= np.sum(val)
        meanval += val
        mats.append(val)
        l1, l2 = loop
        chrom, loop1_start, _ = l1
        start1 = loop1_start-binsize*wsz
        end1 = loop1_start+binsize*(wsz+1)

        chrom, loop2_start, _ = l2
        start2 = loop2_start-binsize*wsz
        end2 = loop2_start+binsize*(wsz+1)
        
        for name, bw in zip(('trl', 'ctcf', 'pser', 'cp190', 'zelda'), (loopset.pybw_trl, loopset.pybw_ctcf, loopset.pybw_pser5,
                                                                    loopset.pybw_cp190, loopset.pybw_zelda)):
            v1 = bw.values(chrom, start1, end1)
            v1 = np.reshape(v1, (length, binsize))
            v1 = np.sum(v1, axis=1)

            L_bws_mean[name] += v1

            v2 = bw.values(chrom, start2, end2)
            v2 = np.reshape(v2, (length, binsize))
            v2 = np.sum(v2, axis=1)
            R_bws_mean[name] += v2

            # L_bws[name].append(bw.values(chrom, start1, end1))
            # R_bws[name].append(bw.values(chrom, start2, end2))

        c+=1
    for name in ('trl', 'ctcf', 'pser', 'cp190', 'zelda'):
        L_bws_mean[name] /= c
        R_bws_mean[name] /= c
    meanval /= c

    return mats, L_bws_mean, R_bws_mean, 


def loop_to_otherloop(loops):
    res = {}
    for key in loops.keys():
        for partner in loops[key]:
            res[(key, partner)] = 1
    return res

def otherloop_to_loop(otherloops):
    loops = {}
    for loop in otherloops.keys():
        l1, l2 = loop
        loops.setdefault(l1, [])
        loops[l1].append(l2)
    return loops

def loop_to_symloop(otherloops):
    loops = {}
    for loop in otherloops.keys():
        l1, l2 = loop
        loops.setdefault(l1, [])
        loops.setdefault(l2, [])
        loops[l1].append((l1, l2))
        loops[l2].append((l2, l1))
    return loops

def get_a_not_b(dic1, dic2):
    res = {}
    for key1 in dic1.keys():
        v1 = dic1.get(key1)
        if v1:
            if dic2.get(key1):
                pass
            else:
                res[key1] = v1
    return res

def get_anc_values(ancsoi, bw, wsz):
    vals = []
    for anc in ancsoi:
        chrom, s, e = anc
        res = e-s
        s2 = s - wsz*res
        e2 = e + wsz*res
        if s2 < 0:
            s2 = 0
            print(f"s2 < 0")
        if e2 > len(chromos[chrom]):
            e2 = len(chromos[chrom])
            print(f"e2 > chromsize")
        val = np.asarray(bw.values(chrom, s2, e2))
        val[np.isnan(val)] = 0
        vals.append(np.sum(val))
    return vals

def get_anchors(loops):
    anchors = {}
    for key in loops.keys():
        anchors[key[0]] = 1
        anchors[key[1]] = 1
    return anchors


def make_notnet_df(loopset, subgraphs, ind2anc, wsz = 5):
    network_df = pd.DataFrame()

    agrees = []
    rand_agrees = []
    anclist = list(loopset.anchors.keys())
    ancsois = []
    inds = []

    all_vals = {}
    mean_vals = {}

    rand_all_vals = {}
    rand_mean_vals = {}
    graph_inds = []
    for i, g in enumerate(subgraphs):
        ancsoi = [ind2anc[node] for node in g.nodes]
        ancsois.append(ancsoi)
        inds.append([node for node in g.nodes])
        randanc = [random.choice(anclist) for a in range(len(ancsoi))]
        graph_inds.append(i)
        assert len(randanc)==len(ancsoi)
        for name in loopset.peak_names.keys():
            if name != "peak_housekeeping":
                peak_dic = getattr(loopset, name)
                vals = np.asarray([peak_dic.get(anc, 0) > 0 for anc in ancsoi])
                randvals = np.asarray([peak_dic.get(anc, 0) >0 for anc in randanc])

                mean_vals.setdefault(name, [])
                mean_vals[name].append(np.mean(vals))

                rand_mean_vals.setdefault(name, [])
                rand_mean_vals[name].append(np.mean(randvals))

                all_vals.setdefault(name, [])
                all_vals[name].append(vals)

                rand_all_vals.setdefault(name, [])
                rand_all_vals[name].append(randvals)
            else:
                peak_dic = getattr(loopset, name)
                vals = np.asarray([peak_dic.get(anc, 0)  for anc in ancsoi])
                randvals = np.asarray([peak_dic.get(anc, 0) for anc in randanc])

                mean_vals.setdefault(name, [])
                mean_vals[name].append(np.mean(vals))

                rand_mean_vals.setdefault(name, [])
                rand_mean_vals[name].append(np.mean(randvals))

                all_vals.setdefault(name, [])
                all_vals[name].append(vals)

                rand_all_vals.setdefault(name, [])
                rand_all_vals[name].append(randvals)                

    for name, val in mean_vals.items():
        df_col = f'mean_{name}'
        network_df[df_col] = val

    for name, val in rand_mean_vals.items():
        df_col = f'rand_mean_{name}'
        network_df[df_col] = val

    for name, val in all_vals.items():
        df_col = f'all_{name}'
        network_df[df_col] = val                

    for name, val in rand_all_vals.items():
        df_col = f'rand_all_{name}'
        network_df[df_col] = val                

    network_df['ancs'] = ancsois
    network_df['inds'] = inds
    network_df['graph_inds'] = graph_inds

    return network_df


def make_net_df(loopset, subgraphs, ind2anc, wsz = 5):
    network_df = pd.DataFrame()

    agrees = []
    rand_agrees = []
    anclist = list(loopset.anchors.keys())
    ancsois = []
    inds = []

    all_vals = {}
    mean_vals = {}

    rand_all_vals = {}
    rand_mean_vals = {}
    graph_inds = []
    for i, g in enumerate(subgraphs):
        ancsoi = [ind2anc[node] for node in g.nodes]
        ancsois.append(ancsoi)
        inds.append([node for node in g.nodes])
        randanc = [random.choice(anclist) for a in range(len(ancsoi))]
        graph_inds.append(i)
        assert len(randanc)==len(ancsoi)
        for name in loopset.peak_names.keys():
            if name != "peak_housekeeping":
                peak_dic = getattr(loopset, name)
                vals = np.asarray([peak_dic.get(anc, 0) > 0 for anc in ancsoi])
                randvals = np.asarray([peak_dic.get(anc, 0) >0 for anc in randanc])

                mean_vals.setdefault(name, [])
                mean_vals[name].append(np.mean(vals))

                rand_mean_vals.setdefault(name, [])
                rand_mean_vals[name].append(np.mean(randvals))

                all_vals.setdefault(name, [])
                all_vals[name].append(vals)

                rand_all_vals.setdefault(name, [])
                rand_all_vals[name].append(randvals)
            else:
                peak_dic = getattr(loopset, name)
                vals = np.asarray([peak_dic.get(anc, 0)  for anc in ancsoi])
                randvals = np.asarray([peak_dic.get(anc, 0) for anc in randanc])

                mean_vals.setdefault(name, [])
                mean_vals[name].append(np.mean(vals))

                rand_mean_vals.setdefault(name, [])
                rand_mean_vals[name].append(np.mean(randvals))

                all_vals.setdefault(name, [])
                all_vals[name].append(vals)

                rand_all_vals.setdefault(name, [])
                rand_all_vals[name].append(randvals)                

    for name, val in mean_vals.items():
        df_col = f'mean_{name}'
        network_df[df_col] = val

    for name, val in rand_mean_vals.items():
        df_col = f'rand_mean_{name}'
        network_df[df_col] = val

    for name, val in all_vals.items():
        df_col = f'all_{name}'
        network_df[df_col] = val                

    for name, val in rand_all_vals.items():
        df_col = f'rand_all_{name}'
        network_df[df_col] = val                

    network_df['ancs'] = ancsois
    network_df['inds'] = inds
    network_df['graph_inds'] = graph_inds

    return network_df
        

def make_subgraphs(graph, S, most_edges, l =5):
    # print("YOU NEED TO GET CONNECTED COMPONENTS RIGHT AFTER BIG_SUBGRAPHS!!! OTHERWISE YOU ARE MESSING IT UP!!!!! ")
    # print("HUGE ERROR HUGE ERROR, YOU ARE ONLY SUPPOSED TO HAVE GRAPHS OF SIZE >= 4! ")
    big_subgraphs = []
    small_subgraphs = []
    print("Making HSC subgraphs")
    for ind in most_edges:
        too = nx.Graph(S[ind])
        t = labelled_HCS(too, l=l)
        new_subgraphs, pruned_subgraphs = get_big_subgraphs(t, n = 4)
        for subg in new_subgraphs:
            big_subgraphs.append(graph.subgraph(subg))
        for subg in pruned_subgraphs:
            small_subgraphs.append(graph.subgraph(subg))
    print("Removing single-node connectors")
    bsubs2 = []
    smallsubs2 = []
    for i, g in enumerate(big_subgraphs):
        comps = list(nx.biconnected_components(g))
        doo = 0
        for comp in comps:
            bsubs2.append(graph.subgraph(comp))    
    return big_subgraphs, bsubs2, small_subgraphs

def categorize_graphs(S):
    sizes = []
    for s in S:
        sizes.append(s.size()) # Get the ones with lots of nodes

    edge_sizes = []
    for s in S:
        adj_mat = nx.adjacency_matrix(s).toarray()
        adj_mat = np.ravel(adj_mat[adj_mat > 0])
        if (np.sum(adj_mat)) == 0:
            edge_sizes.append(0)
            continue
        else:
            metr = scipy.stats.hmean(np.ravel(adj_mat))    
            edge_sizes.append(metr) # Get ones with lots of edges
        
    largest = np.argsort(sizes)[::-1]
    most_edges = np.argsort(edge_sizes)[::-1]
    return largest, most_edges

def make_OG_graph(conn_matrix, method, ):
    if method=='genomicdistance':
        matmax = np.max(conn_matrix)
        adjacency = 1-conn_matrix/matmax
        adjacency[adjacency == 1] = 0
    elif method=='loopvalue':
        matmax = np.median(conn_matrix[conn_matrix > 0])
        adjacency = conn_matrix/matmax
    graph = nx.from_numpy_matrix(adjacency)
    return graph, adjacency


def make_graph(conn_matrix, method, thresh = .2):
    if method=='genomicdistance':
        matmax = np.max(conn_matrix)
        adjacency = 1-conn_matrix/matmax
        adjacency[adjacency == 1] = 0
        adjacency[adjacency < thresh] = 0
    elif method=='loopvalue':
        matmax = np.median(conn_matrix[conn_matrix > 0])
        adjacency = conn_matrix/matmax
        adjacency[adjacency < thresh] = 0
    graph = nx.from_numpy_matrix(adjacency)
    return graph, adjacency

# def peak_recursion(loop, peaks, res, buffer):
#     loop_start, loop_end = loop[1]-buffer*res, loop[2] + buffer*res
#     assert(loop_start < loop_end)
#     midpt = len(peaks)//2
#     if type(peaks[0][2]) == int or type(peaks[0][2]) == float:
#         method = 'ints'
#         l = 0
#     else:
#         method = 'strings'
#         l = []
#     if len(peaks) < 50:
#         for peak in peaks:
#             trl_start, trl_end = peak[:2]
#             if (trl_start > loop_start and trl_start < loop_end) or (trl_end < loop_end and trl_end > loop_start) or \
#                 (trl_start < loop_start and trl_end > loop_end):
#                 if method == 'ints':
#                     l += peak[2]
#                 elif method=='strings':
#                     l.append(peak[2])
#         return l
#     else:
#         peak = peaks[midpt]
#         start = peak[0]
#         end = peak[1]
#         if start > loop_end:
#             return peak_recursion(loop, peaks[:midpt+2], res, buffer)
#         elif end < loop_start:
#             return peak_recursion(loop, peaks[midpt-2:], res, buffer)
#         else:
#             for peak in peaks:
#                 trl_start, trl_end = peak[:2]
#                 if trl_start > loop_end:
#                     break
#                 if (trl_start > loop_start and trl_start < loop_end) or (trl_end < loop_end and trl_end > loop_start) or \
#                     (trl_start < loop_start and trl_end > loop_end):
#                     if method == 'ints':
#                         l += peak[2]
#                     elif method=='strings':
#                         l.append(peak[2])
#             return l


def peak_recursion(anchor, peaks, res, buffer):
    loop_start, loop_end = anchor[1]-buffer*res, anchor[2] + buffer*res
    assert(loop_start <= loop_end)
    if type(peaks[0][2]) == int or type(peaks[0][2]) == float:
        method = 'ints'
        l = []
    else:
        method = 'strings'
        l = []

    midpt = len(peaks)//2
    while (peaks[midpt][0] > loop_end) and (len(peaks) > 10):
        peaks = peaks[:midpt]
        midpt = len(peaks)//2


    peaks = sorted(peaks, key=lambda x: x[1])
    while (peaks[midpt][1] < loop_start) and (len(peaks) > 10):
        peaks = peaks[midpt:]
        midpt = len(peaks)//2

    for peak in peaks:
        trl_start, trl_end = peak[:2]
        if (trl_start > loop_start and trl_start < loop_end) or (trl_end < loop_end and trl_end > loop_start) or \
            (trl_start < loop_start and trl_end > loop_end):
            if method == 'ints':
                l.append(peak)
            elif method=='strings':
                l.append(peak)
    return l


def binary_mac_intersection(loop, chrom_macs, res = 400, buffer = 3, method="add"):
    chrom = loop[0]
    peaks = chrom_macs[chrom]
    ls = peak_recursion(loop, peaks, res, buffer)
    if method=='add':
        l=0
    elif method == 'append':
        l = []
    if len(ls) > 0:
        if method == 'add':
            for peak in ls:
                l += peak[2]
        elif method == 'append':
            for peak in ls:
                l.append(peak)
    return l

def get_peak_dict(anchors, chrom_macs,  res = 400, buffer = 3, method='add'):
    for chrom in chrom_macs:
        assert sorted(chrom_macs[chrom]) == chrom_macs[chrom]
    peak_dict = {}
    print(f'buffer size of {buffer} bins is {buffer*res}')
    
    if type(anchors) == dict:
        for anc in anchors.keys():
            chrom = anc[0]
            if not chrom_macs.get(chrom):
                continue
            is_peak = binary_mac_intersection(anc, chrom_macs, res, buffer, method=method)
            if (is_peak != 0) and (is_peak != []):
                val = is_peak
                if val == []:
                    print(is_peak)
                peak_dict[anc] = val
            else:
                val = 0
    elif type(anchors) == int:
        for anc in anchors:
            chrom = anc[0]
            if not chrom_macs.get(chrom):
                continue
            is_peak = binary_mac_intersection(anc, chrom_macs, res, buffer, method=method)
            if is_peak != 0 and is_peak != []:
                val = is_peak
                peak_dict[anc] = val
            else:
                val = 0
    elif type(anchors) == list:
        for anc in anchors:
            chrom = anc[0]
            if not chrom_macs.get(chrom):
                continue
            is_peak = binary_mac_intersection(anc, chrom_macs, res, buffer, method=method)
            if is_peak != 0 and is_peak != []:
                val = is_peak
                peak_dict[anc] = val
            else:
                val = 0

    return peak_dict

def split_peaks(df, bw,):
    dups = df.duplicated(subset=['chrom', 'start', 'end'], keep=False)
    for v, group in df[dups].groupby(['chrom', 'start', 'end']):
        chrom = group['chrom'].iloc[0]
        ground_start = group['start'].iloc[0]
        ground_end = group['end'].iloc[0]

        group.iloc[0, :]
        l = len(group)
        if l > 1:
            for i in range(l):
                if i != l-1:
                    index = group.iloc[i, :].name
                    next_index = group.iloc[i+1, :].name
                    current_summit = group['summit'].iloc[i] + ground_start
                    next_summit = group['summit'].iloc[i+1] + ground_start
                    between_values = np.asarray(bw.values(chrom, current_summit, next_summit))
                    valley = np.argmin(between_values) + current_summit
                    current_end = valley
                    next_start = valley+1
                    df.at[index, 'end']  = current_end
                    df.at[next_index, 'start']  = next_start
    return df

# Used on original bedgraph file
def df_to_dict(df):
    peaks = {}
    for chrom in chromos.keys():
        peaks[chrom] = sorted(df[df['chrom'] == chrom][['start', 'end', 'signal']].values.tolist())
    return peaks


# peaks = {}
# for i, peak in enumerate(sorted_peaks[:5000-1]):
#     if clusters.labels_[i] != 0:
#         continue
#     else:
#         chrom, start, end = peak[1:4]
#         peak_vals = fetch_bw(loopset_400_small.pybw_pc_rj, chrom, start, end)
#         max_point = np.argmax(peak_vals)
#         new_start = start + max_point
#         window = 10
#         new_end = new_start + window
#         new_start = new_start - window

#         peaks.setdefault(chrom, [])
#         peaks[chrom].append((new_start, new_end, _))

# for c in peaks:
#     peaks[c] = sorted(peaks[c])
    
# places = []
# place_to_peak = {}
# for peak in sorted_peaks_vp:    
#     chrom, start, end = peak[1:4]
#     peak_vals = fetch_bw(loopset_400_small.pybw_pc_vp, chrom, start, end)
#     max_point = np.argmax(peak_vals)
#     new_start = start + max_point
#     window = 10
#     new_end = new_start + window
#     new_start = new_start - window
#     places.append((chrom, new_start, new_end))
#     place_to_peak[(chrom, new_start, new_end)] = peak[1:4]


# def anc_peaks_and_dfs(prot_name, chrom_macs, loopset, buffer, res, align_start = True):
#     print(buffer, res)
#     peak_ancs = get_peak_dict(loopset.anchors, chrom_macs[f"chrom_macs_{prot_name}"], buffer, res)
#     sorted_peaks = []
#     for peak_anc in peak_ancs:
#         chrom, s, e = peak_anc
#         s = s - buffer*res
#         e = e + buffer*res
#         v = peak_ancs[peak_anc]
#         dist = e - s
#         sorted_peaks.append((dist, chrom, s, e, v))
#     sorted_peaks = sorted(sorted_peaks)
#     dfs = make_bw_dfs(loopset, sorted_peaks, prot_name, align_start = align_start)    

#     ag = AgglomerativeClustering(n_clusters=6)
#     clusters = ag.fit(dfs[prot_name])
#     order = np.argsort(clusters.labels_)
#     return sorted_peaks, dfs, clusters, order

from sklearn.cluster import AgglomerativeClustering
import pybedtools as pbt


def peaks_and_dfs(prot_name, bigwigs, summits, anchors, align_start = False):
    # peaks = loopset.pbts['prot_name']
    # anchor_peaks = pbt_anchors.slop(b=1000, g=pbt.genome_registry.dm6).intersect(peaks, u=True)
    intersecting = anchors.intersect(summits[prot_name], u=True)
    dfs = make_bw_dfs(prot_name, bigwigs, summits, intersecting, align_start = align_start)
    ag = AgglomerativeClustering(n_clusters=6)
    clusters = ag.fit(dfs[prot_name])
    order = np.argsort(clusters.labels_)
    return intersecting, dfs, clusters, order



def make_bw_dfs(prot_name, bigwigs, summits, intersecting, align_start = True, bw_list=None, window=200, res=400):
    dfs = {}
    lists = {}
    peaks = []

    veclength = window*2+res
    for i, bigwig_name in enumerate(list(bigwigs)):        
        bigwig = bigwigs[bigwig_name]
        lists[bigwig_name] = np.zeros((len(intersecting), veclength))
        for j, peak in enumerate(intersecting):    
            chrom, start, end = peak[:3]
            start, end = map(int, [start, end])         
            if align_start == True:
                peak_vals = fetch_bw(bigwig, chrom, start, end)  
                max_point = np.argmax(peak_vals)
                new_peak = pbt.BedTool(f'{chrom} {start+max_point} {start+max_point}', from_string=True).slop(b=window+res//2, g='../dm6/dm6.chrom.sizes')[0]

            else:
                new_peak = pbt.BedTool(f'{chrom} {start} {end}', from_string=True).slop(window, g='../dm6/dm6.chrom.sizes')[0]

            chrom, new_start, new_end = new_peak[:3]
            new_start, new_end = map(int, [new_start, new_end])         


            vals = fetch_bw(bigwig, chrom, new_start, new_end)
            assert len(vals) <= veclength
            start_ind = (veclength-len(vals))//2

            lists[bigwig_name][j, start_ind:len(vals)+start_ind] = vals
            if i == 0:
                peaks.append((chrom, new_start, new_end))
            # if i==0:
            # except Exception as e:
            #     # print(e)
            #     continue

            #     for i, attr in enumerate(bw_list):
            #         attr = 'pybw_' + attr
            #         try:
            #             new_bw = getattr(loopset, attr)
            #             vals = fetch_bw(new_bw, chrom, new_start, new_end)
            #             name = "_".join(attr.split("_")[1:])
            #             names.add(name)
            #             lists.setdefault(name, [])
            #             lists[name].append(vals)
            #             if i==0:
            #                 peaks.append((chrom, start, end))                
            #         except Exception as e:
            #             print(e)
            #             continue  
    for name in bigwigs:
        dfs[name] = pd.DataFrame(lists[name], index=peaks)
    return dfs  


import seaborn as sns
def graph_peaks(dfs, prot_name, clusters, order):
    labels = pd.Series(clusters.labels_[order])
    node_pal = sns.cubehelix_palette(labels.unique().size)
    node_lut = dict(zip(map(int, labels.unique()), node_pal))
    node_colors = labels.map(node_lut)

    g = sns.clustermap(dfs[prot_name].iloc[order, :], vmax=10, yticklabels=False, 
                       row_colors=node_colors.values,
                      # Turn off the clustering
                        row_cluster=False, col_cluster=False,
                      # Make the plot look better when many rows/cols
                      linewidths=0, xticklabels=False,)


# # Used
# def peaks_to_dict(df):
#     peaks = {}
#     for chrom in chromos.keys():
#         peaks[chrom] = df[df['chrom'] == chrom][['start', 'end', ]].values.tolist()
#         peaks[chrom] = sorted(peaks[chrom])
#     return peaks



import subprocess

def motifs(name, p = .0005, tf_name = 'long_trl'):
    a = f"fimo --thresh {p} --oc ../tfs/{name}_{tf_name}_out ../tfs/{tf_name}.meme ../tfs/{name}.txt"
    subprocess.call(a, shell=True) 
    motifs = {}
    
    with open(f'../tfs/{name}_{tf_name}_out/fimo.txt', 'r') as f:
        for i, line in enumerate(f.readlines()):
            if i == 0 or 'MA' not in line[:2]: # Ignore header
                continue 
            if i == 1:
#                 print(f'starting with {line}')
                pass
            line = line.strip()
            values = line.split('\t')
            [loop, start, stop, strand] = values[1:5]
            loop = loop.split("|")
            loop[1] = int(loop[1])
            loop[2] = int(loop[2])
            loop = tuple(loop)
            chrom = loop[0]
            motif_start = loop[1] + int(start)
            motif_end = loop[1] + int(stop)
            motifs.setdefault(chrom, [])
            motifs[chrom].append((motif_start, motif_end, strand))
    return motifs


from itertools import product
def loop_is_near(loop, loops):
    key1 = list(loops.keys())[0][0]
    res1 = key1[2] - key1[1]
    
    res2 = loop[0][2]-loop[0][1]
    assert res2==res1, loop
    
    res = res2
    l1 = loop[0]
    l2 = loop[1]

    chrom1 = l1[0]
    anch1 = l1[1:]
    chrom2 = l2[0]
    assert(chrom2==chrom1)
    anch2 = l2[1:]
    for d1, d2 in product(range(-2, 3), range(-2, 3)):
        new_anch1 = anch1[0]+d1*res, anch1[1] + d1*res
        new_anch2 = anch2[0]+d2*res, anch2[1] + d2*res
        new_loop = ((chrom1, *new_anch1), (chrom2, *new_anch2))
        v = loops.get(new_loop)
        if v:
            return v
    return False

# Fetches BW, sets NaN = 0
def fetch_bw(bw, chrom, L, R):
    values = np.asarray(bw.values(chrom, L, R))
    values[np.isnan(values)] = 0
    return values


class LoopSet:
    def __init__(self, other_loops, name, buffer, p, res, chrom_macs, annotations, wsz=4):
        self.res = res
        self.p = p
        self.buffer = buffer
        self.name = name
        self.anchors = get_anchors(other_loops)
        self.collapse_anchors(wsz=wsz, res=self.res)
        self.other_loops = self.collapse_otherloops(other_loops)
        self.annotations = annotations

        loops = otherloop_to_loop(self.other_loops)
        write_loops(loops, name, buffer)        
        self.loops = loops

        symloops = loop_to_symloop(self.other_loops)
        self.symloops = symloops
        
        self.chrom_macs = chrom_macs
        for key in chrom_macs.keys():
            setattr(self, key, chrom_macs[key])

        for key in annotations.keys():
            setattr(self, key, annotations[key])

    def set_all(self):
        self.run()
        self.get_motifs()
        self.set_bws()
        self.set_loop_peaks()
        
    def run(self):
        p = self.p
        self.peak_names = {}
        buffer = self.buffer
        name = self.name
        print("Peak Dicts")
        t1 = time.time()
        for key in self.chrom_macs.keys():
            name = key.split("_")[-1]
            value = get_peak_dict(self.anchors, self.chrom_macs[key], buffer = buffer,)
            peak_name = f'peak_{name}'
            setattr(self, peak_name, value)
            self.peak_names[peak_name] = 1


        print("Setting housekeeping genes")
        gene_dict = get_peak_dict(self.anchors, self.annotations['genes'], buffer = buffer, method='append')
        for gene in gene_dict:
            gene_dict[gene] = np.unique(gene_dict[gene])
        

        self.peak_housekeeping = {}
        self.peak_names["peak_housekeeping"] = 1
        housekeeping = pd.read_csv("../for_annotations/fly_features.csv", usecols = [0, 25])
        housekeeping_dict = annotation_to_dict(housekeeping, 'is_housekeeping')
        for anc in self.anchors:
            genes = gene_dict.get(anc, [])
            houses = 0
            counter = 0
            for gene in genes:
                houses += housekeeping_dict.get(gene, 0)
            self.peak_housekeeping[anc] = houses

        self.gene_dict = gene_dict
                
    def get_motifs(self):
        print("Motifs")
        t1 = time.time()        
        # self.trl_motif = motifs(self.loops, self.name, tf_name = 'long_trl', p = self.p)
        # self.ctcf_motif = motifs(self.loops, self.name, tf_name = 'ctcf', p = self.p)
        t1f = time.time()
        print(t1f - t1)
        
    def set_loop_peaks(self):
        self.loop_peak_trl = peak_with_loops(self.other_loops, self.peak_trl)
        self.loop_peak_pser = peak_with_loops(self.other_loops, self.peak_pser)
        self.loop_peak_ctcf = peak_with_loops(self.other_loops, self.peak_ctcf)
        self.loop_peak_zelda = peak_with_loops(self.other_loops, self.peak_zelda)

        loop_peak_trl = self.loop_peak_trl
        loop_peak_pser = self.loop_peak_pser
        
        loop_peak_ctcf = self.loop_peak_ctcf      
        
        one_trl_peak = {}
        both_trl_peak = {}
        one_ctcf_peak = {}
        both_ctcf_peak = {}
        one_pser_peak = {}
        both_pser_peak = {}
        for loop in self.other_loops.keys():
            if loop_peak_trl.get(loop) == 1:
                one_trl_peak[loop] = 1
            elif loop_peak_trl.get(loop) == 2:
                both_trl_peak[loop] = 1
            
            if loop_peak_pser.get(loop) == 1:
                one_pser_peak[loop] = 1
            elif loop_peak_pser.get(loop) == 2:
                both_pser_peak[loop] = 1
            
            if loop_peak_ctcf.get(loop) == 1:
                one_ctcf_peak[loop] = 1
            elif loop_peak_ctcf.get(loop) == 2:
                both_ctcf_peak[loop] = 1
                
        self.one_trl_peak = one_trl_peak
        self.both_trl_peak = both_trl_peak
        self.one_ctcf_peak = one_ctcf_peak
        self.both_ctcf_peak = both_ctcf_peak
        self.one_pser_peak = one_pser_peak
        self.both_pser_peak = both_pser_peak        

        
    def set_bws(self):
        print(f"Setting bigwigs using padding: {self.buffer}*{self.res}")
        pybw_trl = pyBigWig.open("../tfs/bws/trl.bw")
        pybw_trl2 = pyBigWig.open("../tfs/bws/GSE152770_GAF_GFP_COMBINED.bw")

        pybw_ctcf = pyBigWig.open("../tfs/bws/ctcf.bw")
        pybw_pser5 = pyBigWig.open("../tfs/bws/pser.bw")
        pybw_pol22 = pyBigWig.open("../tfs/bws/GSE62925_pol2_COMBINED.bw")

        pybw_cp190 = pyBigWig.open("../tfs/bws/cp190.bw")
        pybw_zelda = pyBigWig.open("../tfs/bws/zelda.bw")

        pybw_pc_rj = pyBigWig.open('../tfs/bws/PC_RJ_SRR1636792.bw')
        pybw_pc_vp = pyBigWig.open('../tfs/bws/PC_VP_SRR1636793.bw')
        pybws_dref = pyBigWig.open('../tfs/bws/DREF_SRR1636771.bw')
        pybws_beaf = pyBigWig.open('../tfs/bws/BEAF32_SRR1636749.bw')
        pybw_rad21 = pyBigWig.open('../tfs/bws/rad21_COMBINED.bw')

        pybw_dcp = pyBigWig.open('../tfs/bws/GSE57876_dCP_STARR_COMBINED.bw')
        pybw_hkp = pyBigWig.open('../tfs/bws/GSE57876_hkCP_STARR_COMBINED.bw')
        
        pybw_pc2 = pyBigWig.open('../tfs/bws/GSE60428_PC_COMBINED.bw')
        pybw_ph = pyBigWig.open('../tfs/bws/GSE60428_PH_COMBINED.bw')


        self.pybw_trl = pybw_trl
        self.pybw_trl2 = pybw_trl2

        self.pybw_ctcf = pybw_ctcf
        self.pybw_pser5 = pybw_pser5
        self.pybw_cp190 = pybw_cp190
        self.pybw_zelda = pybw_zelda     
        
        self.pybw_pc_rj = pybw_pc_rj
        self.pybw_pc_vp = pybw_pc_vp
        self.pybw_dref = pybws_dref
        self.pybw_beaf = pybws_beaf
        self.pybw_rad21 = pybw_rad21

        self.pybw_dcp = pybw_dcp
        self.pybw_hkp = pybw_hkp

        self.pybw_pol22 = pybw_pol22
        self.pybw_pc2 = pybw_pc2
        self.pybw_ph = pybw_ph

        bw_names = []
        for file in glob.glob("../histones/bws/*"):
            print(file)
            name = file.split("/")[-1].split(".")[0]
            bw = pyBigWig.open(file)
            attr=f'pybw_{name}'
            setattr(self, attr, bw)
            bw_names.append(attr)
        bw_names += ['pybw_trl', 'pybw_ctcf', 'pybw_pser5', 'pybw_cp190', 'pybw_zelda', 'pybw_pc_rj', 'pybw_pc_vp', 'pybw_dref', 'pybw_beaf']
        self.bw_names = bw_names

        # self.bw_trl = {}
        # self.bw_ctcf = {}
        # self.bw_pser5 = {}
        # self.bw_cp190 = {}
        # self.bw_zelda = {}

        # for anchor in self.anchors.keys():
        #     start, end = anchor[1]-self.buffer*self.res, anchor[2]+self.buffer*self.res
        #     chrom = anchor[0]
        #     maxsize = len(chromos[chrom])
        #     if end > maxsize:
        #         end = maxsize-10
        #         print("Boundary error")
        #     if start < 0:
        #         start = 2
        #         print("Boundary error")
        #     v_trl = np.asarray(pybw_trl.values(chrom, start, end))
        #     v_trl[np.isnan(v_trl)] = 0
        #     self.bw_trl[anchor] = v_trl
            
        #     v_ctcf = np.asarray(pybw_ctcf.values(chrom, start, end))
        #     v_ctcf[np.isnan(v_ctcf)] = 0
        #     self.bw_ctcf[anchor] = v_ctcf

        #     v_pser5 = np.asarray(pybw_pser5.values(chrom, start, end))
        #     v_pser5[np.isnan(v_pser5)] = 0
        #     self.bw_pser5[anchor] = v_pser5

        #     v_cp190 = np.asarray(pybw_cp190.values(chrom, start, end))
        #     v_cp190[np.isnan(v_cp190)] = 0
        #     self.bw_cp190[anchor] = v_cp190

        #     zelda = np.asarray(pybw_zelda.values(chrom, start, end))
        #     zelda[np.isnan(zelda)] = 0
        #     self.bw_zelda[anchor] = zelda

    def collapse_anchors(self, wsz, res):
        print(f"Collapsing with wsz = {wsz}, res = {res}")
        collapsed_ancs = {}
        checked = {}
        anc2canc = {}
        canc2anc = {}
        sorted_anchors = sorted(list(self.anchors.keys()))
        running_anchor_list = []
        for i, anchor in enumerate(sorted_anchors):
            current_chrom, current_start = anchor[0], anchor[1]
            if i != len(sorted_anchors) - 1:
                next_anchor = sorted_anchors[i+1]
                next_chrom, next_start = next_anchor[0], next_anchor[1]
                if current_chrom == next_chrom:
                    if (next_start - current_start)/res < wsz:
                        running_anchor_list.append(anchor)
                    else:
                        running_anchor_list.append(anchor)
                        starts = list(zip(*running_anchor_list))[1]
                        canc_start = starts[len(starts)//2]
                        canc = (current_chrom, canc_start, canc_start + res)
                        canc2anc[canc] = []
                        canc2anc[canc] = running_anchor_list
                        for anc in running_anchor_list:
                            anc2canc[anc] = canc
                        running_anchor_list = []    
                else:
                    running_anchor_list.append(anchor)
                    starts = list(zip(*running_anchor_list))[1]
                    canc_start = starts[len(starts)//2]
                    canc = (current_chrom, canc_start, canc_start + res)
                    canc2anc[canc] = []
                    for anc in running_anchor_list:
                        anc2canc[anc] = canc
                        canc2anc[canc].append(anc)
                    running_anchor_list = []
            else:
                running_anchor_list.append(anchor)
                starts = list(zip(*running_anchor_list))[1]
                canc_start = starts[len(starts)//2]
                canc = (current_chrom, canc_start, canc_start + res)
                canc2anc[canc] = []
                for anc in running_anchor_list:
                    anc2canc[anc] = canc
                    canc2anc[canc].append(anc)
                running_anchor_list = []            
        self.anc2canc = anc2canc
        self.canc2anc = canc2anc
        new_anchors = {}
        for new_anc in self.canc2anc:
            new_anchors[new_anc] = 1
        self.anchors = new_anchors


    def collapse_otherloops(self, otherloops):
        print(f"Collapsing loops")
        collapsed_otherloops = {}
        for loop, val in otherloops.items():
            l1, l2 = loop
            newl1, newl2 = self.anc2canc[l1], self.anc2canc[l2]
            collapsed_otherloops[(newl1, newl2)] = val
        return collapsed_otherloops
    

class RandomLoopSet(LoopSet):
    def __init__(self, other_loops, name, buffer, p, res, chrom_macs, annotations, shift, ):
        print("STILL NEED TO EDGE-CASE: WHAT IF SHIFT IS TOO BIG?")
        print(f"Shifting {name} by {shift}")
        other_loops = uniform_shift_loops(other_loops, shift=shift)
        loops = otherloop_to_loop(other_loops)
        super().__init__(other_loops, name, buffer, p, res, chrom_macs, annotations)


def uniform_shift_loops(otherloops, shift = 10000):
    print(f'Shifting by {shift}')
    newloops = {}
    for loop in otherloops:
        l1, l2 = loop
        chrom1 = l1[0]
        chrom2 = l2[0]
        assert(chrom1==chrom2)
        maxsize = len(chromos[chrom1])
        chrom = chrom1
        newl1 = l1[1]+shift, l1[2]+shift
        newl2 = l2[1]+shift, l2[2]+shift
        if newl2[1] > maxsize:
            newl1 = l1[1]-shift, l1[2]-shift
            newl2 = l2[1]-shift, l2[2]-shift
        if newl1[0] < 0:
            continue
        oldv = otherloops[loop]
        newloop = ((chrom, newl1[0], newl1[1]), (chrom, newl2[0], newl2[1]))
        if loop_is_near(newloop, otherloops):
#             print("overlap!")
            pass
        else:
            newloops[newloop] = oldv
    return newloops


def make_peaks(prot_name, chrom_macs):
    peaks = []
    chrom_macs_pc_rj = chrom_macs[f"chrom_macs_{prot_name}"]
    for chrom in chrom_macs_pc_rj:
        for peak in chrom_macs_pc_rj[chrom]:
            dist = peak[1]-peak[0]
            signal = peak[2]
            peaks.append((dist, chrom, peak[0], peak[1], signal))
    sorted_peaks = sorted(peaks)[::-1]    
    return sorted_peaks



def peak_with_loops(other_loops, peaks):
    peak_loops = {}
    c=0
    for loop in other_loops.keys():
        anch1, anch2 = loop
        if peaks.get(anch1) and peaks.get(anch2):
            peak_loops[loop] = 2
        elif (peaks.get(anch1) and (not peaks.get(anch2))) or (peaks.get(anch2) and (not peaks.get(anch1))):
            peak_loops[loop] = 1
        elif (not peaks.get(anch1)) and (not peaks.get(anch2)):
            peak_loops[loop] = 0
        else:
            print("SCRUBB")
        c+=1
    return peak_loops
    

def fetch_anchor_window(anchor, wsz, res):
    co = cooler.Cooler(f'../All_nc14.mcool::resolutions/{res}')
    coolmats = co.matrix(balance=True)
    chrom, start, end = anchor[0], anchor[1], anchor[2]

    res = co.binsize
    s = res
    v = coolmats.fetch((chrom, start-wsz*s, start+(wsz+1)*s,))
    return v


def rinse_loops_from_nan(loops, nan_cutoff):
    looplist = list(loops.keys())
    c=0
    badvs = []
    badloops = []
    co = cooler.Cooler('../All_nc14.mcool::resolutions/800')
    coolmats = co.matrix(balance=True)
    print(len(looplist))
    for k, loop in enumerate(looplist):
        if k % (len(looplist)//10) == 0:
            print(k, len(looplist))
        
        chrom, start, end = loop[0][0], loop[0][1], loop[1][1]

        wsz = 12
        s = 800
        res = co.binsize
        try:
            v = coolmats.fetch(f'{chrom}: {start-wsz*s}-{start+(wsz+1)*s}',f'{chrom}: {end-wsz*s}-{end+(wsz+1)*s}')
        except Exception as e:
            print(e, chrom, start, end)
            continue
        if np.sum(np.isnan(v)) > nan_cutoff:
            c+=1
            badvs.append(v)
            badloops.append(loop)
        
    for loop in badloops:
        del loops[loop]
    return badloops, badvs





def make_conn_matrix(loopset, method='distance'):
    print(f"Using method {method}")
    anc2ind = {}
    ind2anc = {}
    for i, anchor in enumerate(loopset.anchors):
        anc2ind[anchor] = i
        ind2anc[i] = anchor

    n_anchors = len(loopset.anchors)
    conn_matrix = np.zeros((n_anchors, n_anchors))

    for anc in loopset.anchors:
        loops_with_anc = loopset.symloops[anc]
        for loop in loops_with_anc:
            anc1, anc2 = loop
            res = anchor[2]-anchor[1]
            if method=='genomicdistance':
                dist = np.abs(anc1[1]-anc2[1])
                i1, i2 = anc2ind[anc1], anc2ind[anc2]
                conn_matrix[i1, i2] = dist                    
            elif method=='loopvalue':
                loop = (anc1, anc2)
                i1, i2 = anc2ind[anc1], anc2ind[anc2]
                val = fetch_loop_window(loop, 5, res=res)
                conn_matrix[i1, i2] = np.sum(val)
    return conn_matrix, anc2ind, ind2anc    




# Align Loops that are nearby
def align_loops(loop_list, wsz = 4):
    res = 0
    for loops in loop_list:
        key1 = list(loops.keys())[0][0]

        res1 = key1[2] - key1[1]
        if res == 0:
            res = res1
        elif res != res1:
            print(res)
            print(res1)
            print("Different loop resolutions.")
            raise Exception()
            
    checked = {}
    for i in range(len(loop_list)):
        loops = deepcopy(loop_list[i])
        for loop in loops.keys():
            if checked.get(loop) == 1:
                continue
            else:
                ps = []
                ps.append((loop[0][1], loop[1][1]))
                checked[loop] = 1
                l1 = loop[0]
                l2 = loop[1]
                chrom1 = l1[0]
                anch1 = l1[1:]
                chrom2 = l2[0]
                assert(chrom2==chrom1)
                anch2 = l2[1:]

                for other_loops in loop_list[i:]:
                    for d1, d2 in product(range(-wsz, wsz+1), range(-wsz, wsz+1)):
                        if d1 == 0 and d2 == 0:
                            continue
                        new_anch1 = anch1[0]+d1*res, anch1[1] + d1*res
                        new_anch2 = anch2[0]+d2*res, anch2[1] + d2*res
                        new_loop = ((chrom1, *new_anch1), (chrom2, *new_anch2))
                        v = other_loops.get(new_loop)
                        if v:
                            ps.append(deepcopy((new_loop[0][1], new_loop[1][1])))
                            del other_loops[new_loop]
                            other_loops[loop] = v
                p1s, p2s = list(zip(*ps))
                position = int(np.round(np.mean(p1s)/res)*res), int(np.round(np.mean(p2s)/res)*res)
                
                new_loop = ((chrom1, position[0], position[0] + res), (chrom1, position[1], position[1] + res))
                
                for other_loops in loop_list[i:]:
                    v = other_loops.get(loop)
                    if v:
                        del other_loops[loop]
                        other_loops[new_loop] = v




import networkx as nx
import numpy as np

"""Python implementation of basic HCS
Implementation of Highly Connected Subgraphs (HCS) clustering which is introduced by "Hartuv, E., & Shamir, R. (2000).
 A clustering algorithm based on graph connectivity. Information processing letters, 76(4-6), 175-18"
 
Based on NetworkX and Numpy
Notation:
    G = Graph
    E = Edge
    V = Vertex
    
    |V| = Number of Vertices in G
    |E| = Number of Edges in G
"""


def create_example_graph():
    """Create example graph used in the paper
    :return: NetworkX Graph
    """
    v = {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7, 'I': 8, 'X': 9, 'Y': 10, 'Z': 11}

    adjacency = np.zeros(shape=(12, 12), dtype=np.uint8)
    adjacency[v['A'], [v['B'], v['C'], v['D']]] = 1
    adjacency[v['B'], [v['A'], v['D'], v['E'], v['Y']]] = 1
    adjacency[v['C'], [v['A'], v['D'], v['E']]] = 1
    adjacency[v['D'], [v['A'], v['B'], v['C'], v['E']]] = 1
    adjacency[v['E'], [v['B'], v['C'], v['D'], v['F']]] = 1
    adjacency[v['F'], [v['E'], v['Y'], v['G'], v['I'], v['H']]] = 1
    adjacency[v['G'], [v['Z'], v['F'], v['I'], v['H']]] = 1
    adjacency[v['H'], [v['F'], v['G'], v['I']]] = 1
    adjacency[v['I'], [v['H'], v['F'], v['G']]] = 1
    adjacency[v['X'], [v['Y'], v['Z']]] = 1
    adjacency[v['Y'], [v['B'], v['X'], v['Z'], v['F']]] = 1
    adjacency[v['Z'], [v['X'], v['Y'], v['G']]] = 1

    return nx.from_numpy_matrix(adjacency)


def highly_connected(G, E, l=5):
    """Checks if the graph G is highly connected
    Highly connected means, that splitting the graph G into subgraphs needs more than 0.5*|V| edge deletions
    This definition can be found in Section 2 of the publication.
    :param G: Graph G
    :param E: Edges needed for splitting G
    :return: True if G is highly connected, otherwise False
    """
    # np.median(nx.adjacency_matrix(G))
    # print(nx.adjacency_matrix(G).shape[0])
    # vals = 0
    # for edge in list(E):
    #     i = (np.where(G.nodes)==edge[0])[0]
    #     j = (np.where(G.nodes)==edge[1])[0]
    #     val =

    return len(E) > len(G.nodes) / l


def remove_edges(G, E):
    """Removes all edges E from G
    Iterates over all edges in E and removes them from G
    :param G: Graph to remove edges from
    :param E: One or multiple Edges
    :return: Graph with edges removed
    """

    for edge in E:
        G.remove_edge(*edge)
    return G


def HCS(G, l):
    """Basic HCS Algorithm
    cluster labels, removed edges are stored in global variables
    :param G: Input graph
    :return: Either the input Graph if it is highly connected, otherwise a Graph composed of
    Subgraphs that build clusters
    """

    E = nx.algorithms.connectivity.cuts.minimum_edge_cut(G)
    if not highly_connected(G, E, l):
        G = remove_edges(G, E)
        sub_graphs = [G.subgraph(c).copy() for c in nx.connected_components(G)]
        if len(sub_graphs) == 2:
            H = HCS(sub_graphs[0], l)
            _H = HCS(sub_graphs[1], l)

            G = nx.compose(H, _H)

    return G


def improved_HCS(G):
    """
    Implements improvements mentioned in the paper
    1. Iterated HCS
    2. Singleton adoption
    3. Removing Low Degree Vertices
    """
    pass


def labelled_HCS(G, l):
    """
    Runs basic HCS and returns Cluster Labels
    :param G: Input graph
    :return: List of cluster assignments for the single vertices
    """
    G = deepcopy(G)

    _G = HCS(G, l)

    sub_graphs = (G.subgraph(c).copy() for c in nx.connected_components(_G))


    labels = {}

    for _class, _cluster in enumerate(sub_graphs, 1):
        c = list(_cluster.nodes)   
        labels[_class] = c
    return labels

def get_big_subgraphs(labels, n = 5):
    keys = []
    notkeys = []
    for key, v in labels.items():
        if len(v) > n:
            keys.append(v)
        else:
            notkeys.append(v)
    return keys, notkeys
