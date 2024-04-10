


import scipy
import scipy.ndimage
import numpy as np 
import matplotlib.pyplot as plt

def make_resdict(df, mnase_bws, nearest_gene, strand='+', delta=3000):
	resdict = {}
	for cluster in range(4):
		subdf = df[df['Flybase_IDs'].isin(nearest_gene[cluster][1]['gene_id'])]
		subdf = subdf[subdf['strand'] == strand]
		chrom, s, e = subdf['seqname'], subdf['start'], subdf['end']
		for cond, bw in mnase_bws.items():
			v = bw.stackup(chrom, s-delta, s + delta)
			w = np.nanmean(v, axis=0)
			w = scipy.ndimage.gaussian_filter1d(w, sigma=10)
			resdict[f'{cluster}_{cond}'] = v
	return resdict

def make_and_plot_resdict(df, mnase_bws, nearest_gene, strand='+', delta=3000):
	resdict = {}
	for cluster in range(4):
		subdf = df[df['Flybase_IDs'].isin(nearest_gene[cluster][1]['gene_id'])]
		subdf = subdf[subdf['strand'] == strand]
		chrom, s, e = subdf['seqname'], subdf['start'], subdf['end']
		
		fig, ax = plt.subplots(1, 2, figsize=(8, 4))
		for cond, bw in mnase_bws.items():
			v = bw.stackup(chrom, s-delta, s + delta)
			w = np.nanmean(v, axis=0)
			w = scipy.ndimage.gaussian_filter1d(w, sigma=10)
			if 'negative' in cond:
				a = ax[0]
				a.plot(w, label=cond)
				a.set_title(f"TSS from Cluster: {cluster}, \n MNase_read=negative")
			else:
				a = ax[1]
				a.plot(w, label=cond)
				a.set_title(f"TSS from Cluster: {cluster}, \n MNase_read=positive")
			resdict[f'{cluster}_{cond}'] = v
		for c, a in enumerate(ax):
			a.legend(bbox_to_anchor=(.5, -.25), loc='center')
			a.set_ylim([.8, 1.8])
	return resdict, fig

def make_resdict_atac(df, bwdict, cond, nearest_gene, strand='+', delta=3000):
	resdict = {}
	bw = bwdict[cond]
	for cluster in range(4):
		subdf = df[df['Flybase_IDs'].isin(nearest_gene[cluster][1]['gene_id'])]
		subdf = subdf[subdf['strand'] == strand]
		chrom, s, e = subdf['seqname'], subdf['start'], subdf['end']
		v = bw.stackup(chrom, s-delta, s + delta)
		w = np.nanmean(v, axis=0)
		w = scipy.ndimage.gaussian_filter1d(w, sigma=10)
		resdict[f'{cluster}_{cond}'] = v
	return resdict

def make_resdict_atac_boundaries(df, bwdict, cond, delta=3000):
	resdict = {}
	bw = bwdict[cond]
	for cluster in df['itemRgb'].unique():
		subdf = df[df['itemRgb'] == cluster]
		chrom, s, e = subdf['chrom'], subdf['start'], subdf['end']
		v = bw.stackup(chrom, s-delta, s + delta)
		w = np.nanmean(v, axis=0)
		w = scipy.ndimage.gaussian_filter1d(w, sigma=10)
		resdict[f'{cluster}_{cond}'] = v
	return resdict


def make_and_plot_gaf(gaf_df, mnase_bws, strand='+', tf='gaf'):
	resdict = {}
	chrom, s, e = gaf_df['chrom'], gaf_df['start'], gaf_df['end']
	fig, ax = plt.subplots(2, 1, figsize=(8, 8))
	for cond, bw in mnase_bws.items():
		v = bw.stackup(chrom, s-100, s + 100)
		w = np.nanmean(v, axis=0)
		if 'negative' in cond:
			a = ax[0]
			a.plot(w, label=cond)
			a.set_title(f"Peaks from: {tf}, \n MNase_read=negative")
		else:
			a = ax[1]
			a.plot(w, label=cond)
			a.set_title(f"Peaks from: {tf}, \n MNase_read=positive")
		resdict[f'{cond}'] = v
	for c, a in enumerate(ax):
		a.legend(bbox_to_anchor=(.5, -.25), loc='center')
		a.set_ylim([0, 2.2])
	return resdict, fig