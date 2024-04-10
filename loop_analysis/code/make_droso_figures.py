import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.cm as cm
import matplotlib.colors as colors
import pandas as pd

import seaborn as sns
import scipy
def make_order(matrix):
    linkage = scipy.cluster.hierarchy.linkage(matrix, method='average', metric='cosine')
    dendro = scipy.cluster.hierarchy.dendrogram(linkage, no_plot=True,
                                        color_threshold=-np.inf)
    order = dendro['leaves']
    return order   

def make_order2_cosine(matrix):
    linkage = scipy.cluster.hierarchy.linkage(matrix, method='average', metric='cosine')
    dendro = scipy.cluster.hierarchy.dendrogram(linkage, no_plot=True,
                                        color_threshold=-np.inf)

    ordering = scipy.cluster.hierarchy.optimal_leaf_ordering(linkage, matrix)
    order = scipy.cluster.hierarchy.leaves_list(ordering)
    # order = dendro['leaves']
    return order


 
def make_order2(matrix):
    linkage = scipy.cluster.hierarchy.linkage(matrix, method='average', metric='euclidean')
    dendro = scipy.cluster.hierarchy.dendrogram(linkage, no_plot=True,
                                        color_threshold=-np.inf)

    ordering = scipy.cluster.hierarchy.optimal_leaf_ordering(linkage, matrix)
    order = scipy.cluster.hierarchy.leaves_list(ordering)
    # order = dendro['leaves']
    return order        

def make_order2_and_cluster_complete(matrix, n_clusters):
    linkage = scipy.cluster.hierarchy.linkage(matrix, method='complete', metric='euclidean')
    dendro = scipy.cluster.hierarchy.dendrogram(linkage, no_plot=True,
                                        color_threshold=-np.inf)

    ordering = scipy.cluster.hierarchy.optimal_leaf_ordering(linkage, matrix)
    order = scipy.cluster.hierarchy.leaves_list(ordering)
    po_ = fcluster(linkage, t=n_clusters, criterion='maxclust')-1

    return order, po_


def make_order2_and_cluster_ward(matrix, n_clusters):
    linkage = scipy.cluster.hierarchy.linkage(matrix, method='ward', metric='euclidean')
    dendro = scipy.cluster.hierarchy.dendrogram(linkage, no_plot=True,
                                        color_threshold=-np.inf)

    ordering = scipy.cluster.hierarchy.optimal_leaf_ordering(linkage, matrix)
    order = scipy.cluster.hierarchy.leaves_list(ordering)
    po_ = fcluster(linkage, t=n_clusters, criterion='maxclust')-1

    # order = dendro['leaves']
    return order, po_


def make_order2_and_cluster(matrix, n_clusters):
    linkage = scipy.cluster.hierarchy.linkage(matrix, method='average', metric='euclidean')
    dendro = scipy.cluster.hierarchy.dendrogram(linkage, no_plot=True,
                                        color_threshold=-np.inf)

    ordering = scipy.cluster.hierarchy.optimal_leaf_ordering(linkage, matrix)
    order = scipy.cluster.hierarchy.leaves_list(ordering)
    po_ = fcluster(linkage, t=n_clusters, criterion='maxclust')-1

    # order = dendro['leaves']
    return order, po_



from scipy.cluster.hierarchy import fcluster
def make_order_and_cluster(mat, n_clusters):
    linkage = scipy.cluster.hierarchy.linkage(mat, method='ward', metric='euclidean')
    dendro = scipy.cluster.hierarchy.dendrogram(linkage, no_plot=True,
                                        color_threshold=-np.inf)
    # ordering = scipy.cluster.hierarchy.optimal_leaf_ordering(linkage, mat)
    # order = scipy.cluster.hierarchy.leaves_list(ordering)

    order = dendro['leaves']
    po_ = fcluster(linkage, t=n_clusters, criterion='maxclust')-1

    return order, po_

def make_order_and_cluster_perf(mat, n_clusters):
    linkage = scipy.cluster.hierarchy.linkage(mat, method='ward', metric='euclidean')
    dendro = scipy.cluster.hierarchy.dendrogram(linkage, no_plot=True,
                                        color_threshold=-np.inf)
    ordering = scipy.cluster.hierarchy.optimal_leaf_ordering(linkage, mat)
    order = scipy.cluster.hierarchy.leaves_list(ordering)

    po_ = fcluster(linkage, t=n_clusters, criterion='maxclust')-1

    return order, po_


def make_df_from_dict(dic):
    df = pd.DataFrame()
    values = []
    labels = []
    for key in dic:
        values = values + list(dic[key])
        labels = labels + [key]*len(dic[key])
    df['values']=values
    df['labels']=labels
    return df    


def get_cs_from_inds(nc14_vs_nc18, inds, nc18_bw, nc14_bw, window_size = 800, use_log = True):
    Ls, Rs = [], []
    for j, row in nc14_vs_nc18[inds].iterrows():
        chrom, s1, s2 = row['BIN1_CHR'], row['BIN1_START'], row['BIN2_START']
        L = (chrom, s1-window_size, s1+window_size)
        R = (chrom, s2-window_size, s2+window_size)
        Ls.append(L)
        Rs.append(R)

    chroms1, ss1, es1 = list(zip(*Ls))
    chroms2, ss2, es2 = list(zip(*Rs))    
    if use_log ==True:
        vs_18_L = np.log2(nc18_bw.stackup(chroms1, ss1, es1)+1e-1)
        vs_14_L = np.log2(nc14_bw.stackup(chroms1, ss1, es1)+1e-1)
    else:
        vs_18_L = nc18_bw.stackup(chroms1, ss1, es1)
        vs_14_L = nc14_bw.stackup(chroms1, ss1, es1)

    cs_L = np.nanmean(vs_14_L - vs_18_L, axis=1)

    chroms1, ss1, es1 = list(zip(*Ls))
    chroms2, ss2, es2 = list(zip(*Rs))    
    if use_log ==True:
        vs_18_R = np.log2(nc18_bw.stackup(chroms2, ss2, es2)+1e-1)
        vs_14_R = np.log2(nc14_bw.stackup(chroms2, ss2, es2)+1e-1)
    else:
        vs_18_R = nc18_bw.stackup(chroms2, ss2, es2)
        vs_14_R = nc14_bw.stackup(chroms2, ss2, es2)

    cs_R = np.nanmean(vs_14_R - vs_18_R, axis=1)
    cs = cs_L + cs_R
    return cs


def get_peaks_from_inds(nc14_vs_nc18, inds, ancset,):
    Ls, Rs = [], []
    peakvals = []
    for j, row in nc14_vs_nc18[inds].iterrows():
        chrom, s1, e1, _, s2, e2 = row['BIN1_CHR'], row['BIN1_START'], row['BIN1_END'], row['BIN2_CHR'], row['BIN2_START'], row['BIN2_END']
        L = (chrom, s1, e1)
        R = (chrom, s2, e2)

        peakvals.append((L in ancset) + (R in ancset))
    return peakvals

arr = np.asarray
def get_raw_from_inds(nc14_vs_nc18, inds, bw, slop=0, use_log = True):
    Ls, Rs = [], []
    for j, row in nc14_vs_nc18[inds].iterrows():
        chrom, s1, e1, _, s2, e2 = row['BIN1_CHR'], row['BIN1_START'], row['BIN1_END'], row['BIN2_CHR'], row['BIN2_START'], row['BIN2_END']
        L = (chrom, s1-slop, e1+slop)
        R = (chrom, s2-slop, e2+slop)
        Ls.append(L)
        Rs.append(R)

    chroms1, ss1, es1 = list(zip(*Ls))
    chroms2, ss2, es2 = list(zip(*Rs))
    if use_log ==True:
        vs_18_L = np.log2(np.nanmean(bw.stackup(chroms1, ss1, es1)+.1, axis=1))
    else:
        vs_18_L = np.nanmean(bw.stackup(chroms1, ss1, es1), axis=1)

    cs_L = vs_18_L

    chroms1, ss1, es1 = list(zip(*Ls))
    chroms2, ss2, es2 = list(zip(*Rs))
    if use_log ==True:
        vs_18_R = np.log2(np.nanmean(bw.stackup(chroms2, ss2, es2)+.1, axis=1))
    else:
        vs_18_R = np.nanmean(bw.stackup(chroms2, ss2, es2), axis=1)
    cs_R = vs_18_R

    cs = cs_L + cs_R
    return cs

def make_base_plot(loops, df, base_valdict, cond1, cond2):
	c = make_c_for_loops(loops, df)
	z = np.linspace(1e-4, 1e-1)
	fig, ax = plt.subplots(1, 1, figsize=(3, 3))
	a, b = arr(base_valdict[cond1]), arr(base_valdict[cond2])
	x = a+b
	y = np.log(b/a)
	ax.scatter(a, b, c=c, cmap=cm.bwr,  vmin=-.9, vmax=1.5)
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlabel(f"Base,  {cond1}")
	ax.set_ylabel(f"Base, {cond2}")
	ax.set_title(f'Raw read counts: {cond2} vs. {cond1}')
	ax.plot(z, z)


def make_oe_plot(loops, df, oe_valdict, cond1, cond2, c=None):
        if c is None:
            c = make_c_for_loops(loops, df)
        z = np.linspace(0, 6)
        a, b = arr(oe_valdict[cond1]), arr(oe_valdict[cond2])
        fig, ax = plt.subplots(1, figsize=(3, 3))
        ax.scatter(a, b, c=c, cmap=cm.bwr,  vmin=-.9, vmax=1.5)
        ax.plot(z, z, color='red')
        ax.set_xlabel(f"O/E,  {cond1}")
        ax.set_ylabel(f"O/E, {cond2}")
        ax.set_title(f'O/E: {cond2} vs. {cond1}')

def make_c_from_df(df, pco=.05):
    cs = []
    for _, row in df.iterrows():
        if (row['padj'] < pco) and (row['l2fc'] > 0):
            cs.append(1)
        elif (row['padj'] < pco) and (row['l2fc'] < 0):
            cs.append(-1)
        elif (row['padj'] > pco):
            cs.append(0)
    return cs

def make_c_for_loops(loops, df ):
    cs = []
    for loop in loops:
        chrom, a1, a2, chrom, b1, b2 = loop
        row = df[(df['BIN1_CHR'].values == chrom) & 
                     (df['BIN1_START'].values == int(a1)) & 
                     (df['BIN2_START'].values == int(b1))]
        if len(row) == 0:
            cs.append(0)
            continue
        if (row['padj'] < .05).values and (row['l2fc'] > 0).values:
            cs.append(1)
        elif (row['padj'] < .05).values and (row['l2fc'] < 0).values:
            cs.append(-1)
        elif (row['padj'] > .05).values:
            cs.append(0)
        else:
            cs.append(0)
            continue
    return cs

import itertools
def cdf(fcs, names, xlabel="", folder="", filename="tmp", xmin=None, xmax=None):
    
    lens ={}
    if xmin is None:
        xmin = np.min(fcs[0])
    if xmax is None:
        xmax = np.max(fcs[0])
    xs = np.linspace(xmin, xmax, 1000)
    fc_dict = {}
    for i, fc in enumerate(fcs):
        name = names[i]
        fc_dict[name] = []
        for j, x in enumerate(xs):
            fc_dict[name].append((np.asarray(fc) < x).mean())
        lens[name] = len(fc)
    fig, axs = plt.subplots(figsize=(8, 6))
    for name in fc_dict:
        axs.plot(xs, fc_dict[name], label=f'{name}, {lens[name]}')
    axs.legend()

    axs.set_xlabel(xlabel)
    axs.set_xlim([xmin, xmax])
    plt.tight_layout()
    fig.savefig(f'plots/{folder}/{filename}')


    dict_for_ks = dict(zip(names, fcs))
    for pair in itertools.combinations(list(dict_for_ks.keys()), 2):
        name1, name2 = pair
        ks_test = scipy.stats.ks_2samp(dict_for_ks[name1], dict_for_ks[name2])
        print(name1, name2, ks_test[1])
    # pval = ks_test[1]
    # axs.set_title(pval.round(3))
    # axs[1].plot(xs, ys_treatment-ys_control, label=f'{name1} - {name2}')

    
def make_decay_correlation_plot(chromo_keys_no_chr, expected_dict, condition_dict):
    n = len(chromo_keys_no_chr)
    n_conditions = len(np.unique(list(condition_dict.values())))
    palette = sns.color_palette('gist_ncar', n_colors = n_conditions)
    
    chrom_to_ind = dict(zip(chromo_keys_no_chr, range(n)))
    condition_to_ind = dict(zip(np.unique(list(condition_dict.values())), range(n_conditions)))


    fig, axs = plt.subplots(3, 2, figsize=(30, 20))
    axs = np.ravel(axs)
    conditions = set()
    for _, replicate1 in enumerate(list(expected_dict.keys())):
        df1 = expected_dict[replicate1]
        for _, replicate2 in enumerate(list(expected_dict.keys())):
            df2 = expected_dict[replicate]
            df1 = df1[df1['region']=='3R'].sort_values('diag')
            df2 = df2[df2['region']=='3R'].sort_values('diag')


            for c, (chrom, x) in enumerate(df.groupby('region')):
                if chrom == 'Y':
                    continue
                ax = axs[chrom_to_ind[chrom]]
                condition = condition_dict[replicate]
                if condition not in conditions:
                    ax.plot(x.diag, x['balanced.avg'], label=condition, color=palette[condition_to_ind[condition]])
                else:
                    ax.plot(x.diag, x['balanced.avg'], color=palette[condition_to_ind[condition]])
                conditions.add(condition)
                
                ax.set_title(chrom)
                ax.set_yscale('log')
                ax.set_xscale('log')
                ax.set_ylim([1e-6, 5e-2])
    for ax in axs:
        ax.legend()
    fig.savefig(f'./plots/qc/decay/decay_balanced.png')


def make_decay_plot(chromo_keys_no_chr, expected_dict, condition_dict):
    n = len(chromo_keys_no_chr)
    n_conditions = len(np.unique(list(condition_dict.values())))
    palette = sns.color_palette('gist_ncar', n_colors = n_conditions)
    
    chrom_to_ind = dict(zip(chromo_keys_no_chr, range(n)))
    condition_to_ind = dict(zip(np.unique(list(condition_dict.values())), range(n_conditions)))


    fig, axs = plt.subplots(3, 2, figsize=(30, 20))
    axs = np.ravel(axs)
    conditions = set()
    for _, replicate in enumerate(list(expected_dict.keys())):
        df = expected_dict[replicate]
        for c, (chrom, x) in enumerate(df.groupby('region')):
            if chrom == 'Y':
                continue
            ax = axs[chrom_to_ind[chrom]]
            condition = condition_dict[replicate]
            if condition not in conditions:
                ax.plot(x.diag, x['balanced.avg'], label=condition, color=palette[condition_to_ind[condition]])
            else:
                ax.plot(x.diag, x['balanced.avg'], color=palette[condition_to_ind[condition]])
            conditions.add(condition)
            
            ax.set_title(chrom)
            ax.set_yscale('log')
            ax.set_xscale('log')
            ax.set_ylim([1e-6, 5e-2])
    for ax in axs:
        ax.legend()
    fig.savefig(f'./plots/qc/decay/decay_balanced.png')

    for _, replicate in enumerate(list(expected_dict.keys())):
        df = expected_dict[replicate]
        for c, (chrom, x) in enumerate(df.groupby('region')):
            if chrom == 'Y':
                continue
            ax = axs[chrom_to_ind[chrom]]
            condition = condition_dict[replicate]
            if condition not in conditions:
                ax.plot(x.diag, x['count.avg'], label=condition, color=palette[condition_to_ind[condition]])
            else:
                ax.plot(x.diag, x['count.avg'], color=palette[condition_to_ind[condition]])
            conditions.add(condition)
            
            ax.set_title(chrom)
            ax.set_yscale('log')
            ax.set_xscale('log')
            ax.set_ylim([1e-4, 1e3])
    for ax in axs:
        ax.legend()
    fig.savefig(f'./plots/qc/decay/decay_count.png')    


def binding_site_anchor_size(full_anchors, tfs, res=400, shift=50):
    frac_dict = {}
    shifted_frac_dict = {}
    for tf_name in tfs:
        if "H3" in tf_name:
            continue
        fracs = []
        shifted_fracs = []
        tf = tfs[tf_name]
        xs = []
        for i in range(0, 20*res, res):
            frac = len(tf.intersect(full_anchors.slop(b=i, g=chromsizes), u=True))/len(tf)
            shifted_frac = len(tf.intersect(full_anchors.shift(s=shift*res, g=chromsizes).slop(b=i, g=chromsizes), u=True))/len(tf)    
            xs.append(i)
            fracs.append(frac)
            shifted_fracs.append(shifted_frac)
        frac_dict[tf_name]= fracs
        shifted_frac_dict[tf_name] = shifted_fracs
    xs = np.asarray(xs)

    fig, axs = plt.subplots(1, 2, figsize=(10, 4))   
    for tf_name in frac_dict:
        axs[0].plot(xs/res, frac_dict[tf_name], label=tf_name, marker='o')
    axs[0].legend()   
    for tf_name in frac_dict:
        axs[1].plot(xs/res, shifted_frac_dict[tf_name], label=tf_name, marker='o')
    axs[1].legend()
    fig.savefig('./plots/qc/fraction_TF_covered.png')


    frac_dict = {}
    shifted_frac_dict = {}
    for tf_name in tfs:
        if "H3" in tf_name:
            continue        
        fracs = []
        shifted_fracs = []
        tf = tfs[tf_name]
        xs = []
        for i in range(0, 20*res, res):
            frac = len(full_anchors.slop(b=i, g=chromsizes).intersect(tf, u=True))/len(full_anchors)
            shifted_frac = len(full_anchors.shift(s=shift*res, g=chromsizes).slop(b=i, g=chromsizes).intersect(tf, u=True))/len(full_anchors)    
            xs.append(i)
            fracs.append(frac)
            shifted_fracs.append(shifted_frac)
        frac_dict[tf_name]= fracs
        shifted_frac_dict[tf_name] = shifted_fracs
    xs = np.asarray(xs)

    fig, axs = plt.subplots(1, 2, figsize=(10, 4))   
    for tf_name in frac_dict:
        axs[0].plot(xs/res, frac_dict[tf_name], label=tf_name, marker='o')
    axs[0].legend()   
    for tf_name in frac_dict:
        axs[1].plot(xs/res, shifted_frac_dict[tf_name], label=tf_name, marker='o')
    axs[1].legend()
    fig.savefig('./plots/qc/fraction_LOOPS_covered.png') 


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

from copy import deepcopy
from mpl_toolkits.axes_grid1 import ImageGrid
#from statannotations.Annotator import Annotator
def plot_boxplots(keys, bbi_bws, df, cond1, cond2, pco = .1, use_log=True, savepath='./', save=False):
    figs = []
    for key in keys:
        bw = bbi_bws[key]
        all_cs = {}

        inds = ((df['padj'] < pco) & (df['l2fc'] < 0))
        all_cs[f'{cond1} < {cond2}'] = get_raw_from_inds(df, inds,  bw, slop=0, use_log=use_log)

        inds = ((df['padj'] > pco))
        all_cs[f'{cond1} = {cond2}'] = get_raw_from_inds(df, inds,  bw, slop=0, use_log=use_log)

        inds = ((df['padj'] < pco) & (df['l2fc'] > 0))
        all_cs[f'{cond1} > {cond2}'] = get_raw_from_inds(df, inds,  bw, slop=0, use_log=use_log)

        tmp = make_df_from_dict(all_cs)

        pvalues = []
        comparisons = []
        for key1 in all_cs:
            for key2 in all_cs:
                if key1 <= key2:
                    continue
                else:
                    pval = scipy.stats.ranksums(all_cs[key1], all_cs[key2])[1]
                    pvalues.append(pval)
                    comparisons.append((key1, key2))


        fig, ax = plt.subplots(figsize=(5, 4))
        plotting_parameters = {
            'data':    tmp,
            'x':       'labels',
            'y':       'values',
            'width' : .4,
        }

        sns.boxplot(ax=ax, **plotting_parameters)

        annotator = Annotator(ax, comparisons, **plotting_parameters)
        annotator.set_pvalues(pvalues)
        annotator.annotate()
        ax.set_ylabel(f"Raw {key} signal")
        ax.set_xlabel(f"")

        ax.set_title(f"Raw {key} signal")
        figs.append(fig)
    return figs


def plot_apa(og_apa, keys=None, pc_ratio=.1, cond2name={}, vmin=-2, vmax=2, delbad=True, badfunc=None):
    vs = og_apa[0]
    if keys is None:
        keys = vs.keys(); print(keys, len(keys))
    n = len(keys);
    fig = plt.figure(figsize=(3*n, 3))
    ax = ImageGrid(fig, 111, nrows_ncols=(1, n), axes_pad=1, share_all=False, cbar_location="right", cbar_mode="single", cbar_size="7%", cbar_pad=0.15,)
    
    counts = og_apa[0]
    expecteds = og_apa[1]

    for c, (cond) in enumerate(keys):
        count, exp = og_apa[0][cond][:, 3:-3, 3:-3], og_apa[1][cond][:, 3:-3, 3:-3]
        pc = exp*.1
        oe = ((count+pc)/(exp+pc))
        if badfunc:
            bad = badfunc(oe)
        else:
            bad = ((oe>50).sum(axis=(1, 2)) > 8)

        if c == 0:
            allbad = bad
        else:
            allbad[bad==1]=1
    if delbad==False:
        allbad = np.zeros_like(allbad)
    for c, (cond) in enumerate(keys):
        a = ax[c]
        count = counts[cond]
        exp = expecteds[cond]
        pc = exp*pc_ratio
        oe = ((count+pc)/(exp+pc))[~allbad]
        oe = np.clip(oe, -20, 20)

        full_v = np.log2(np.nanmean(oe, axis=0)[3:-3, 3:-3])
        if c < n-1:
            a.matshow(full_v, cmap=cm.bwr, vmin=vmin, vmax=vmax)
        elif c == n-1:
            pos = a.matshow(full_v, cmap=cm.bwr, vmin=vmin, vmax=vmax)
            a.cax.colorbar(pos)
            a.cax.toggle_label(True)
            a.cax.set_yticks([vmin, 0, vmax])

        a.set_title(cond2name.get(cond, cond))
        a.text(0, 33, f"n={len(oe)}", fontsize=12)

    for c, a in enumerate(ax):
        a.set_xticks([-.5, 17, 34.5])
        a.set_xticklabels(['-14kb', 'Anchor', '+14kb'])        
        a.set_yticks([-.5, 17, 34.5])
        a.set_yticklabels(['-14kb', 'Anchor', '+14kb'])
        a.tick_params(axis="x", bottom=True, top=False, labelbottom=True, labeltop=False)
        a.tick_params(axis="y", left=True, labelleft=True)
    return fig, allbad


def scatter_raw_apa(full_apa, bad, c1, c2, cond2name, d=2):
    avoid = bad
    m1, m2 = full_apa[0][c1][~avoid], full_apa[0][c2][~avoid]
    exp1, exp2 = full_apa[1][c1][~avoid], full_apa[1][c2][~avoid]

    x, y = (m1)[:, 20-d:20+d+1, 20-d:20+d+1], (m2)[:, 20-d:20+d+1, 20-d:20+d+1]
    x = np.nanmean(x, axis=(1,2))
    y = np.nanmean(y, axis=(1,2))
    z = np.linspace(0, .006)
    fig, ax = plt.subplots(1, 2, figsize=(6, 3))
    ax[0].scatter(x, y, s=3)
    ax[0].plot(z, z, color='red')
    ax[0].set_xlim([-.0005, .006])
    ax[0].set_ylim([-.0005, .006])
    ax[0].set_xlabel(f"{cond2name[c1]}")
    ax[0].set_ylabel(f"{cond2name[c2]}")
    ax[0].set_title(f"Interaction frequency at loops: \n {cond2name[c2]} vs. {cond2name[c1]}")

    z = np.linspace(0, max(max(x), max(y)))
    ax[1].scatter(x, y, s=3)
    ax[1].plot(z, z, color='red')
    ax[1].set_xlabel(f"{cond2name[c1]}")
    ax[1].set_ylabel(f"{cond2name[c2]}")
    ax[1].set_title(f"Interaction frequency at loops: \n {cond2name[c2]} vs. {cond2name[c1]}")
    plt.tight_layout()
    return fig, x, y

def plot_decay_curves(expected_dict_all, cond2name, chromlist, cdict={}):
    n = len(chromlist)
    m = n//2 + n%2
    fig, axs = plt.subplots(2, m, figsize=(4*n, 8))
    axs = np.ravel(axs)
    for c, chrom in enumerate(chromlist):
        ax = axs[c]
        for key, df in expected_dict_all.items():
            x = deepcopy((df[df['region'] == chrom]['diag'].values*800)/1_000)
            y = deepcopy(df[df['region'] == chrom]['balanced.avg'].values)
            bad = np.isnan(x) | np.isnan(y)
            x, y = x[~bad], y[~bad]
            y = scipy.ndimage.gaussian_filter1d(y, sigma=4)
            ax.plot(x, y, c=cdict.get(key, None), label=cond2name[key])
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylim(1e-8, 5e-2)
        ax.legend()
        ax.set_xlabel('Distance')
        ax.set_ylabel('Interaction Frequency')
        if 'hr' in chrom:
           chrom = chrom[3:]
        ax.set_title(f'Chr{chrom}')
    plt.tight_layout()
    return fig


import scipy
def apa(loops, cool_dict, expected_dict, res, wsz=20):
    roll_sz = 3
    print(roll_sz)
    mat_dict = {}
    expecteds = {}
    # pc_rat = .5
    pc_rat = .1
    n=2*wsz+1

    for condition in cool_dict:
        cool = cool_dict[condition]
        expected_df = expected_dict[condition]
        mat_dict[condition] = []
        places = []
        for c, i in enumerate(loops):
            l1, l2 = i[:3], i[3:6]
            chrom = l1[0]
            s = ((int(l1[1]) + int(l1[2]))//2)//res*res
            e = ((int(l2[1]) + int(l2[2]))//2)//res*res
            if (e-s)/res < wsz:
                continue
            places.append([i[:6]])
            expected = make_expected(expected_df, chrom, s, e, wsz, res=res, balanced=True)

            new_l1 = (l1[0], int(s)-res*wsz, int(s)+res*(wsz+1))

            new_l2 = (l2[0], int(e)-res*wsz, int(e)+res*(wsz+1))
            val = np.asarray(cool.matrix(balance=True).fetch(new_l1, new_l2))
            if val.shape != (2*wsz+1, 2*wsz+1):
                print(l1, l2)
                print(val.shape)
                continue

            obs_div_exp = np.log2((val + expected*pc_rat)/(expected + expected*pc_rat))
            # obs_div_exp = val

            midpoint = len(val)//2
            smooth_mat_input = deepcopy(obs_div_exp)
            smooth_mat_input[np.isnan(smooth_mat_input)] = 0
            smooth_mat = scipy.ndimage.gaussian_filter(smooth_mat_input, sigma=1)
            smooth_mat[np.isnan(obs_div_exp)] = -10
            v = np.asarray(np.unravel_index(np.nanargmax(smooth_mat[midpoint-roll_sz:midpoint+roll_sz+1, midpoint-roll_sz:midpoint+roll_sz+1]), (roll_sz*2+1, roll_sz*2+1)))
            left, right = roll_sz-v[0], roll_sz-v[1]

            newmat = np.roll(obs_div_exp, left, axis=0)
            newmat = np.roll(newmat, right, axis=1)
            zeromat = np.zeros(newmat.shape)
            zeromat[roll_sz:n-roll_sz, roll_sz:n-roll_sz] = newmat[roll_sz:n-roll_sz, roll_sz:n-roll_sz]

            obs_div_exp = zeromat
            obs_div_exp[np.isinf(obs_div_exp)] = np.nan
            mat_dict[condition].append(obs_div_exp)
            expecteds[condition] = expected
            c += 1
    for key in mat_dict:
        mat_dict[key] = np.asarray(mat_dict[key])
    return mat_dict, expecteds, places

def apa_noroll(loops, cool_dict, expected_dict, wsz=20, skip_small=True):
    roll_sz = 3
    print(roll_sz)
    mat_dict = {}
    expecteds = {}
    # pc_rat = .5
    pc_rat = .1
    n=2*wsz+1

    for condition in cool_dict:
        cool = cool_dict[condition]
        res = cool.info['bin-size']
        expected_df = expected_dict[condition]
        mat_dict[condition] = []
        expecteds[condition] = []

        places = []
        for c, i in enumerate(loops):

            l1, l2 = i[:3], i[3:6]
            chrom = l1[0]
            s = ((int(l1[1]) + int(l1[2]))//2)//res*res
            e = ((int(l2[1]) + int(l2[2]))//2)//res*res
            if ((e-s)/res < wsz) and skip_small==True:
                continue
            places.append([i[:6]])
            expected = make_expected(expected_df, chrom, s, e, wsz, res=res, balanced=True)

            new_l1 = (l1[0], int(s)-res*wsz, int(s)+res*(wsz+1))

            new_l2 = (l2[0], int(e)-res*wsz, int(e)+res*(wsz+1))
            try:
                val = np.asarray(cool.matrix(balance=True).fetch(new_l1, new_l2))
            except Exception as e:
                print("Could not fetch in ", condition, ' error: ', e)
                bad = np.zeros((2*wsz+1, 2*wsz+1))
                mat_dict[condition].append(bad)
                expecteds[condition].append(bad)
                continue
            if val.shape != (2*wsz+1, 2*wsz+1):
                print(l1, l2)
                print(val.shape)
                bad = np.zeros((2*wsz+1, 2*wsz+1))
                mat_dict[condition].append(bad)
                expecteds[condition].append(bad)
                continue
            diag = (e-s)//res
            
            inds = -np.arange(0, 2*wsz+1)
            good_inds = np.subtract.outer(inds, inds)+diag >= 2
            
            val[~good_inds] = np.nan
            expected[~good_inds] = np.nan

            mat_dict[condition].append(val)
            expecteds[condition].append(expected)
            c += 1
            
    for key in mat_dict:
        mat_dict[key] = np.asarray(mat_dict[key])
    for key in expecteds:
        expecteds[key] = np.asarray(expecteds[key])

    return mat_dict, expecteds, places


def boundary_pileup_analysis(boundaries, cool_dict, 
                            expected_dict, wsz=20, 
                            skip_small=True, pc_rat = .1, add_chr = False):
    mat_dict = {}
    expecteds = {}
    n=2*wsz+1
    for condition in cool_dict:
        cool = cool_dict[condition]
        expected_df = expected_dict[condition]
        mat_dict[condition] = []
        expecteds[condition] = []
        res = cool.info['bin-size']

        places = []
        for c, i in enumerate(boundaries):
            l1 = i[:3]
            chrom = l1[0]
            if add_chr == True:
                chrom = 'chr' + chrom
            s = int(l1[1])
            e = int(l1[2])
            expected = make_expected(expected_df, chrom, s, e, wsz, res=res, balanced=True)
            new_l1 = (chrom, int(s)-res*wsz, int(s)+res*wsz)
            try:
                val = np.asarray(cool.matrix(balance=True).fetch(new_l1))
            except Exception as e:
                continue
            if val.shape != (2*wsz+1, 2*wsz+1):
                # print(l1)
                # print(val.shape)
                continue
            diag = 0
            inds = -np.arange(0, 2*wsz+1)
            good_inds = np.abs(np.subtract.outer(inds, inds)+diag) >= 0
            val[np.isinf(val)] = np.nan
            val[~good_inds] = np.nan
            mat_dict[condition].append(val)
            expecteds[condition].append(expected)
            c += 1
            places.append([i[:3]])
 
    for key in mat_dict:
        mat_dict[key] = np.asarray(mat_dict[key])
    for key in expecteds:
        expecteds[key] = np.asarray(expecteds[key])

    return mat_dict, expecteds, places
 

# def shifted_tf_binding(partners, mustache_tcon_anchors, mustache_treg_anchors, mustache_nonsig_anchors):
#     for name in partners:
#         po_ = partners[name]

#         xs, ys = [], []
#         xs_2, ys_2 = [], []
#         xs_3, ys_3 = [], []
#         for i in range(-20, 20, 2):
#             shift = i
#             places = mustache_tcon_anchors.shift(s=i*5000, g='./annotations/chromsizes').intersect(po_, u=True)
#             x = shift
#             y = len(places)
#             xs.append(x)
#             ys.append(y)

#             places2 = mustache_treg_anchors.shift(s=i*5000, g='./annotations/chromsizes').intersect(po_, u=True)
#             xs_2.append(shift)
#             ys_2.append(len(places2))

#             places3 = mustache_nonsig_anchors.shift(s=i*5000, g='./annotations/chromsizes').intersect(po_, u=True)
#             xs_3.append(shift)
#             ys_3.append(len(places3))
            
#         [xs, ys, xs_2, ys_2, xs_3,  ys_3] = map(np.asarray, [xs, ys, xs_2, ys_2, xs_3, ys_3])
#         fig, axs = plt.subplots( figsize=(6, 4))
#         axs.set_xlabel('Scale = 5kb')
#         axs.plot(xs, ys/len(mustache_tcon_anchors), marker='o', label='Tcon anchors')
#         axs.plot(xs_2, ys_2/len(mustache_treg_anchors), marker='o', label='Treg anchors')
#         axs.plot(xs_3, ys_3/len(mustache_nonsig_anchors), marker='o', label='Non-sig anchors')
#         axs.set_title(name)
#         # ax.set_title('Anchors')
#         axs.legend()
#         fig.savefig(f'plots/shifted_tf_binding/{name}.png')


def examine_deseq_output_unnormalized(deseq_dicts, w, cool_dict, prefix, res):
    for key in deseq_dicts:
        c1, c2 = key.split("_")[0], key.split("_")[2]
        print(c1, c2)
        df = deseq_dicts[key]
        balanced = {}
        nonbalanced = {}
        for i in range(3):
            balanced[i] = []
            nonbalanced[i] = []

        chroms = df['BIN1_CHR']
        ss = df['BIN1_START']
        es = df['BIN2_START']
        lfcs = df.l2fc.values
        padjs = df.padj.values
        diags = (es-ss)//res    
        for i in range(len(chroms)):
            chrom, diag = chroms[i], diags[i]
            s = ss[i]
            e = es[i]
            lfc = np.round(lfcs[i], 2)
            padj = np.round(padjs[i], 2)
            if padj < .05 and lfc  > 0:
                ax_ind = 0
            if padj < .05 and lfc  < 0:
                ax_ind = 1
            if padj > .05:
                ax_ind = 2
            diag = (e-s)//res

            v1 = (cool_dict[c1].matrix(balance=True).fetch((chrom, s-w*res, s+w*res), (chrom, e-w*res, e+w*res)
                                      ).mean())
            v2 = (cool_dict[c2].matrix(balance=True).fetch((chrom, s-w*res, s+w*res), (chrom, e-w*res, e+w*res)
                                      ).mean())

            balanced[ax_ind].append((np.log2(v2), np.log2(v1)))

            v1 = (cool_dict[c1].matrix(balance=False).fetch((chrom, s-w*res, s+w*res), (chrom, e-w*res, e+w*res)
                                      ).mean())

            v2 = (cool_dict[c2].matrix(balance=False).fetch((chrom, s-w*res, s+w*res), (chrom, e-w*res, e+w*res)
                                      ).mean())

            nonbalanced[ax_ind].append((np.log2(v2), np.log2(v1)))

        titles = {
            0 : 'LFC > 0',
            1 : 'LFC < 0',
            2 : 'NS',
        }
        z = np.linspace(-20, 20, 1000)
        fig, ax = plt.subplots(1, 3, figsize=(9, 3))
        for i in nonbalanced:
            x, y = list(zip(*nonbalanced[i]))
            ax[i].scatter(x, y)
            ax[i].plot(z, z)
            ax[i].set_xlim([1, 9])
            ax[i].set_ylim([1, 9])
            ax[i].set_xlabel(f"{c2} (raw)")
            ax[i].set_ylabel(f"{c1} (raw)")
            ax[i].set_title(f'{titles[i]}: {len(x)}')
        fig.suptitle(f"nonbalanced interaction frequency {key}")
        plt.tight_layout()
        fig.savefig(f'./plots/qc/deseq/{prefix}_{c1}_vs_{c2}_nonbalanced.png')
        fig, ax = plt.subplots(1, 3, figsize=(18, 6))    
        for i in balanced:
            x, y = list(zip(*balanced[i]))
            ax[i].scatter(x, y)
            ax[i].plot(z, z)
            ax[i].set_xlim([-13, -4])
            ax[i].set_ylim([-13, -4])
            ax[i].set_xlabel(f"{c2} (balanced)")
            ax[i].set_ylabel(f"{c1} (balanced)")
            ax[i].set_title(f'{titles[i]}: {len(x)}')
        fig.suptitle(f"balanced interaction frequency {key}")
        plt.tight_layout()
        fig.savefig(f'./plots/qc/deseq/{prefix}_{c1}_vs_{c2}_balanced.png')

def examine_deseq_output_div_by_exp(deseq_dicts, w, cool_dict, prefix,
            balanced_expected_value_dict, unbalanced_expected_value_dict, res):
    balanced_vals = {}
    nonbalanced_vals = {}

    for c1 in balanced_expected_value_dict:
        df = deseq_dicts['nc14_vs_nc18']
        balanced = []
        non_balanced = []

        chroms = df['BIN1_CHR']
        ss = df['BIN1_START']
        es = df['BIN2_START']
        diags = (es-ss)//res    
        for i in range(len(chroms)):
            chrom, diag = chroms[i], diags[i]
            s = ss[i]
            e = es[i]
            diag = (e-s)//res
            v1 = (cool_dict[c1].matrix(balance=True).fetch((chrom, s-w*res, s+w*res), (chrom, e-w*res, e+w*res)
                                      ).mean()/balanced_expected_value_dict[c1][(chrom, diag)])

            balanced.append(float(np.log2(v1)))
            v1 = (cool_dict[c1].matrix(balance=False).fetch((chrom, s-w*res, s+w*res), (chrom, e-w*res, e+w*res)
                                      ).mean()/(unbalanced_expected_value_dict[c1][(chrom, diag)]))
            non_balanced.append(float(np.log2(v1)))

        z = np.linspace(0, 10, 1000)
        fig, ax = plt.subplots()
        ax.scatter(balanced, non_balanced)
        ax.set_title(f'{c1}, OBS/EXP')
        ax.set_xlabel(f'Balanced O/E')
        ax.set_ylabel(f'Unbalanced O/E')
        ax.set_xlim([-1.5, 6])
        ax.set_ylim([-1.5, 6])
        fig.savefig(f'./plots/qc/deseq/{prefix}_{c1}_obs_div_exp_balanced_vs_unbalanced.png')
        balanced_vals[c1] = np.asarray(balanced)
        nonbalanced_vals[c1] = np.asarray(non_balanced)

    for key in deseq_dicts:
        c1, c2 = key.split("_")[0], key.split("_")[2]
        print(c1, c2)
        df = deseq_dicts[key]

        v_balanced = balanced_vals[c1]-balanced_vals[c2]
        v_deseq = df.l2fc
        fig, ax = plt.subplots()
        ax.scatter(v_balanced, v_deseq)
        ax.set_ylim(-6, 6)
        ax.set_xlim(-6, 6)
        ax.set_xlabel(f"{c1} O/E - {c2} O/E, balanced")
        ax.set_ylabel("DESeq LFC")
        ax.set_title(f"Balanced: {c1} vs. {c2}")
        z = np.linspace(-4, 4, 100)
        ax.plot(z, z)
        fig.savefig(f'./plots/qc/deseq/{prefix}_{c1}_vs_{c2}_deseq_vs_OE_balanced.png')
   
    for key in deseq_dicts:
        c1, c2 = key.split("_")[0], key.split("_")[2]
        print(c1, c2)
        df = deseq_dicts[key]

        v_nonbalanced = nonbalanced_vals[c1]-nonbalanced_vals[c2]
        v_deseq = df.l2fc
        fig, ax = plt.subplots()
        ax.scatter(v_nonbalanced, v_deseq)
        ax.set_ylim(-6, 6)
        ax.set_xlim(-6, 6)
        ax.set_xlabel(f"{c1} O/E - {c2} O/E, unbalanced")
        ax.set_ylabel("DESeq LFC")
        ax.set_title(f"Non-balanced: {c1} vs. {c2}")
        z = np.linspace(-4, 4, 100)
        ax.plot(z, z)
        fig.savefig(f'./plots/qc/deseq/{prefix}_{c1}_vs_{c2}_deseq_vs_OE_unbalanced.png')



import random
def both_looping_heatmap(loopdict, ancdict, prot_dict):    
    for i, (loop_name, loops_of_interest) in enumerate(loopdict.items()):
        ancsets = {}
        for name in prot_dict:
            ancsets[name] = set()
            for k in ancdict[loop_name].intersect(prot_dict[name], u = True):
                anc = tuple(k[:3])
                ancsets[name].add(anc)
            print(len(ancsets[name]), name)
        listified = list(loops_of_interest)
        anclist_test = []
        left_anchors = []
        right_anchors = []
        for x in listified:
            l1, l2 = x[:3], x[3:6]
            anclist_test.append(tuple(l1))
            anclist_test.append(tuple(l2))    
            left_anchors.append(tuple(l1))
            right_anchors.append(tuple(l2))
        left_anchors = np.asarray(left_anchors)
        right_anchors = np.asarray(right_anchors)

        ancpo_ = list(anclist_test)

        enrichment_df = pd.DataFrame()
        pvalue_df = pd.DataFrame()
        pvalue2_df = pd.DataFrame()

        for partner1_name in list(ancsets.keys()):
            enrichments = []
            pvals = []
            pvals2 = []
            for partner2_name in list(ancsets.keys()):
                if partner1_name < partner2_name:
                    enrichments.append(0)
                    pvals.append(0)
                    pvals2.append(0)
                    continue

                cross_binders = 0
                for loop in listified:
                    l1 = tuple(loop[:3])
                    l2 = tuple(loop[3:6])
                    if ((l1 in ancsets[partner1_name]) and (l2 in ancsets[partner2_name])) or \
                        ((l2 in ancsets[partner1_name]) and (l1 in ancsets[partner2_name])):
                            cross_binders += 1
                print(cross_binders, partner1_name, partner2_name)
                n = 1000
                binding_prop = []
                for k in range(n):
                    permutation = np.random.permutation(len(left_anchors))
                    left_anchors = left_anchors[permutation]
                    rand_cross_binders = 0
                    for _, l1 in enumerate(right_anchors):
                        l2 = left_anchors[_]
                        if ((tuple(l1) in ancsets[partner1_name]) and (tuple(l2) in ancsets[partner2_name])) or \
                            ((tuple(l2) in ancsets[partner1_name]) and (tuple(l1) in ancsets[partner2_name])):
                                rand_cross_binders += 1
                    binding_prop.append(rand_cross_binders/len(listified))
                loop_p = cross_binders/len(listified)
                rand_ps = np.asarray(binding_prop)

                print(np.min(rand_ps), np.max(rand_ps).round(2), loop_p, partner1_name, partner2_name)
                pval2 = min((rand_ps < loop_p).mean(), (rand_ps > loop_p).mean())
                pvals2.append(pval2)

                n = 20000
                rand_cross_binders = 0
                for _ in range(n):
                    randloop = random.sample(ancpo_, 2)
                    l1, l2 = randloop
                    l1, l2 = tuple(l1[:3]), tuple(l2[:3])
                    if ((l1 in ancsets[partner1_name]) and (l2 in ancsets[partner2_name])) or \
                        ((l2 in ancsets[partner1_name]) and (l1 in ancsets[partner2_name])):
                            rand_cross_binders += 1
                pval = one_sample_proportion_significance_test(cross_binders, len(listified), rand_cross_binders/n)


                rand_bothp = rand_cross_binders/n
                real_bothp = cross_binders/len(listified)
                enrichment = real_bothp/ rand_bothp
                print("Enrichment:", enrichment)
                enrichments.append(enrichment)      
                pvals.append(pval)

            enrichment_df[partner1_name] = enrichments
            pvalue_df[partner1_name] = pvals
            pvalue2_df[partner1_name] = pvals2
        enrichment_df.index = list(ancsets.keys())
        pvalue_df.index = list(ancsets.keys())
        pvalue2_df.index = list(ancsets.keys())

        pvalue_df += pvalue_df.T - np.diag(np.diag(pvalue_df))
        pvalue2_df += pvalue2_df.T - np.diag(np.diag(pvalue2_df))

        enrichment_df += enrichment_df.T - np.diag(np.diag(enrichment_df))
        sns.set(font_scale=2)
        enrichment_df[enrichment_df == 0] = 1e-4

        pvalue_df.round(4)
        pvalue2_df.round(4)
        enrichment_df.round(3)

        if i == 0:
            fig = sns.clustermap(np.log2(enrichment_df), figsize=(15, 15),
               vmin= -1, vmax = 1, cmap=cm.bwr, annot=True)
            order_col = deepcopy(fig.dendrogram_col.reordered_ind)
            order_row = deepcopy(fig.dendrogram_row.reordered_ind)        
        else:
            df = enrichment_df.iloc[order_row, order_col]
            fig = sns.clustermap(np.log2(df), figsize=(15, 15),
               vmin=-1, vmax = 1, cmap=cm.bwr, row_cluster=False, col_cluster=False, annot=True)
        fig.cax.set_visible(False)
        # plt.tight_layout()
        fig.savefig(f'./plots/both_looping_heatmap/{loop_name}_enrichment.png')        

        pvalue_df = pvalue_df.iloc[order_row, order_col]
        pval_fig = sns.clustermap(pvalue_df, figsize=(15, 15),
           vmin=0, vmax = .1, cmap=cm.bwr, row_cluster=False, col_cluster=False, annot=True)
        pval_fig.cax.set_visible(False)
        plt.tight_layout()
        pval_fig.savefig(f'./plots/both_looping_heatmap/{loop_name}_pvalue.png')        

        pvalue2_df = pvalue2_df.iloc[order_row, order_col]
        pval_fig = sns.clustermap(pvalue2_df, figsize=(15, 15),
           vmin=0, vmax = .1, cmap=cm.bwr, row_cluster=False, col_cluster=False, annot=True)
        pval_fig.cax.set_visible(False)
        plt.tight_layout()
        pval_fig.savefig(f'./plots/both_looping_heatmap/{loop_name}_pvalue2.png')        

def one_sample_proportion_significance_test(p1_success, p1_tot, p2_success):
    zscore = (p1_success-p2_success)/(np.sqrt(p2_success*(1-p2_success))/p1_tot)
    pval = 2*(1-scipy.stats.norm.cdf(np.abs(zscore), 0, 1))    
    return pval




import bbi
def process_distance_mnase(pref):
    mnase_dict = {
        'peak_4_s10_positive' : bbi.open(pref + '/peak_4_s10_positive.bw'),
        'peak_4_s10_negative' : bbi.open(pref + '/peak_4_s10_negative.bw'),
        'peak_4_nc18_positive' : bbi.open(pref + '/peak_4_nc18_positive.bw'),
        'peak_4_nc18_negative' : bbi.open(pref + '/peak_4_nc18_negative.bw'),
        'peak_4_nc14_positive' : bbi.open(pref + '/peak_4_nc14_positive.bw'),
        'peak_4_nc14_negative' : bbi.open(pref + '/peak_4_nc14_negative.bw'),
        'peak_3_s10_positive' : bbi.open(pref + '/peak_3_s10_positive.bw'),
        'peak_3_s10_negative' : bbi.open(pref + '/peak_3_s10_negative.bw'),
        'peak_3_nc18_positive' : bbi.open(pref + '/peak_3_nc18_positive.bw'),
        'peak_3_nc18_negative' : bbi.open(pref + '/peak_3_nc18_negative.bw'),
        'peak_3_nc14_positive' : bbi.open(pref + '/peak_3_nc14_positive.bw'),
        'peak_3_nc14_negative' : bbi.open(pref + '/peak_3_nc14_negative.bw'),        
        'peak_2_s10_positive' : bbi.open(pref + '/peak_2_s10_positive.bw'),
        'peak_2_s10_negative' : bbi.open(pref + '/peak_2_s10_negative.bw'),
        'peak_2_nc18_positive' : bbi.open(pref + '/peak_2_nc18_positive.bw'),
        'peak_2_nc18_negative' : bbi.open(pref + '/peak_2_nc18_negative.bw'),
        'peak_2_nc14_positive' : bbi.open(pref + '/peak_2_nc14_positive.bw'),
        'peak_2_nc14_negative' : bbi.open(pref + '/peak_2_nc14_negative.bw'),
        'peak_1_s10_positive' : bbi.open(pref + '/peak_1_s10_positive.bw'),
        'peak_1_s10_negative' : bbi.open(pref + '/peak_1_s10_negative.bw'),
        'peak_1_nc18_positive' : bbi.open(pref + '/peak_1_nc18_positive.bw'),
        'peak_1_nc18_negative' : bbi.open(pref + '/peak_1_nc18_negative.bw'),
        'peak_1_nc14_positive' : bbi.open(pref + '/peak_1_nc14_positive.bw'),
        'peak_1_nc14_negative' : bbi.open(pref + '/peak_1_nc14_negative.bw'),
        'peak_0_s10_positive' : bbi.open(pref + '/peak_0_s10_positive.bw'),
        'peak_0_s10_negative' : bbi.open(pref + '/peak_0_s10_negative.bw'),
        'peak_0_nc18_positive' : bbi.open(pref + '/peak_0_nc18_positive.bw'),
        'peak_0_nc18_negative' : bbi.open(pref + '/peak_0_nc18_negative.bw'),
        'peak_0_nc14_positive' : bbi.open(pref + '/peak_0_nc14_positive.bw'),
        'peak_0_nc14_negative' : bbi.open(pref + '/peak_0_nc14_negative.bw'),
    }
    return mnase_dict
