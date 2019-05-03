import sys
import os
import glob
import argparse
import collections
import re

import pandas as pd


patt = re.compile(r'\s\(\d+\-\d+\)$')

def parse_lines(txt):
    header = txt.pop()
    seq_counts = collections.defaultdict(int)
    seq_anno = collections.defaultdict(set)
    anno_seq = collections.defaultdict(set)
    for line in txt:
        els = line.split('\t')
        seq = els.pop(0)
        count = int(els.pop(0))
        seq_counts[seq] += count
        for annotation in els:
            if annotation == '' or annotation == 'low complexity':
                a = 'no annotation'
            else:
                a = patt.sub('', annotation)
            seq_anno[seq].add(a)
            anno_seq[a].add(seq)
    return seq_counts, seq_anno, anno_seq

def make_feature_counts(seq_counts, seq_anno, anno_seq, frac_count=True):
    feature_counts = collections.defaultdict(float)
    for feature, seq_list in anno_seq.items():
        for seq in seq_list:
            count = seq_counts.get(seq, 0.0)
            if frac_count:
                n_assigned_anno = len(seq_anno.get(seq))
                if n_assigned_anno > 0:
                    count = count / n_assigned_anno
                else:
                    raise ValueError
            feature_counts[feature] += count
    return feature_counts

def argparser():
    parser = argparse.ArgumentParser(description='Aggregate unitas features from full annotation files')
    parser.add_argument('--sample-sheet', help='Optional sample sheet. Will subset aggregated table if needed', dest='samples')
    parser.add_argument('-o ', '--output', help='Output filename. Will default to stdout.')
    parser.add_argument('--frac-counts', action='store_true', dest='frac', help='Use fractional counts for sequences that map to multiple annotations')
    parser.add_argument('filenames', nargs='+')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = argparser()
    df_list = []
    for fn in args.filenames:
        sample_id = os.path.dirname(fn).split(os.path.sep)[-1]
        with open(fn) as fh:
            txt = fh.read().splitlines()
        seq_counts, seq_anno, anno_seq = parse_lines(txt)
        counts = make_feature_counts(seq_counts, seq_anno, anno_seq, frac_count=args.frac)
        df = pd.Dataframe.from_dict(counts, orient='index')
        df.columns = [sample_id]
        if args.samples:
            raise NotImplementedError
        df_list.append(df)
    DF = pd.concat(df_list, axis=1, join='outer', sort=False)
    if args.output is None:
        out_fn = sys.stdout
    else:
        out_fn = args.output
    DF.fillna(0, inplace=True)
    DF.index.name = 'feature_id'
    DF.to_csv(out_fn, sep='\t')

