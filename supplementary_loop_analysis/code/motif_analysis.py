import scipy
import scipy.ndimage
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import pybedtools as pbt



def nochr_make_ctcf_peaks(anchors, name, shift=0):
    fasta = pbt.BedTool('./dm6/nochr.dm6.fa')
    coords = anchors.shift(g='./dm6/nochr.sizes', s=shift)
    seq = open(coords.sequence(fi=fasta).seqfn).read()
    with open(f'./meme/for_motifs/{name}', 'w+') as f:
        f.write(seq.upper())
        
def fimo_to_bed(file, out):
    with open(file) as f:
        print_lines = []
        for i, line in enumerate(f.readlines()):
            if i==0:
                continue
            elif line=='\n':
                break
            else:
                vals = line.split('\t')
                motif_name = vals[0]
                chrom, _ = vals[2].split(':')
                loop_start, loop_end = _.split('-')
                loop_start, loop_end = map(int, [loop_start, loop_end])
                motif_start, motif_end = map(int, [vals[3], vals[4]])
                start = loop_start+motif_start
                end = loop_start+motif_end
                strand = vals[5]
                pval = float(vals[6])
                qval = float(vals[7])
                print_lines.append([str(chrom), str(start), str(end), str(motif_name), str(qval), str(strand), str(pval)])
    with open(out, 'w+') as f:
        for line in print_lines:
            f.write("\t".join(line) + "\n")

arr = np.asarray
def get_col(bed, col):
    return arr(list(zip(*bed))[col])

def get_columns(bedtool, n=3):
    new_bedtool = pbt.BedTool.from_dataframe(bedtool.to_dataframe().iloc[:, :n])
    return new_bedtool

def make_motif_df(boundaries_as_bedtool, all_motif_names, motif_enc, boundary_enc, gpath, file='./peaks/motifs/all_anchor_motifs.bed', shift=0):
    n_motifs = len(all_motif_names); n_ancs = len(boundaries_as_bedtool)
    anc_motif_mat = np.zeros((n_ancs, n_motifs))
    anclist = []; motiflist = []
    sub_bedtool = get_columns(boundaries_as_bedtool, n=3)
    for i in sub_bedtool.intersect(pbt.BedTool(file).shift(s=shift, g=gpath), wo=True):
        anc = i[:3]
        motif = i[6]
        anclist.append('_'.join(anc))
        motiflist.append(motif)

    i1s = boundary_enc.transform(anclist)
    i2s = motif_enc.transform(motiflist)
    for i1, i2 in zip(i1s, i2s):
        anc_motif_mat[i1, i2] += 1

    motif_df = pd.DataFrame(anc_motif_mat)
    motif_df.columns = motif_enc.inverse_transform(motif_df.columns.values)
    motif_df.index = boundary_enc.inverse_transform(motif_df.index)
    return motif_df