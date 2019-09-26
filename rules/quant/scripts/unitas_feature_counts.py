import sys
import os
import glob
import argparse
import collections
import re

import pandas as pd


patt = re.compile(r'\s\(\d+\-\d+\)$')

def parse_lines(txt):
    header = txt.pop(0)
    if not 'annotation(s)' in header:
        raise ValueError
    seq_counts = collections.defaultdict(int)
    seq_anno = collections.defaultdict(set)
    anno_seq = collections.defaultdict(set)
    for line in txt:
        els = line.split('\t')
        seq = els[0]
        count = int(els[1])
        remainder = els[2:]
        seq_counts[seq] += count
        if len(els) == 3:
            # no annotation
            #seq_anno[seq].add(seq)
            #anno_seq[seq].add(seq)
            continue
        else:
            for annotation in remainder:
                if annotation == '':
                    continue
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
                anno = list(seq_anno.get(seq))
                n_assigned_anno = len(anno)
                if n_assigned_anno > 0:
                    #if n_assigned_anno == 2:
                    #    print(n_assigned_anno, count, seq_anno[seq], anno)
                    #    sys.exit()
                    count = count / n_assigned_anno
                else:
                    raise ValueError
            feature_counts[feature] += count
    return feature_counts


def parse_annotations(names):
    """find feature class from annotation names.
    """
    anno = []
    for n in names:
        if n.startswith('feature_class'):
            continue
        if 'miR' in n or n.startswith('let'):
            feature = 'miRNA'
            name = n
        elif 'piR' in n:
            feature = 'piRNA'
            name = n
        elif 'tRF' in n:
            feature = n.split()[0]
            name = n
        elif '|' in n:
            feature, name = n.split('|')
        else:
            feature = 'unknown'
            name = n
        anno.append([feature, name])
    df = pd.DataFrame(anno, index=names, columns=['feature_class', 'feature_name'])
    return df

def argparser():
    parser = argparse.ArgumentParser(description='Aggregate unitas features from full annotation files')
    parser.add_argument('filenames', nargs='+')
    parser.add_argument('--sample-sheet', help='Optional tab separated sample sheet with `Sample_ID` in header.\
    Will subset on samples in aggregated table if needed', dest='samples')
    parser.add_argument('--output', required=True, help='Output filename. Required')
    parser.add_argument('--frac-counts', action='store_true', dest='frac', help='Use fractional counts for sequences that map to multiple annotations')
    
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
        df = pd.DataFrame.from_dict(counts, orient='index')
        df.columns = [sample_id]
        if args.samples:
            raise NotImplementedError
        df_list.append(df)
    DF = pd.concat(df_list, axis=1, join='outer', sort=False)
    DF.fillna(0, inplace=True)
    DF.index.name = 'feature_id'
    DF = DF.reindex(DF.mean(axis=1).sort_values(ascending=False).index, axis=0)
    if not os.path.exists(os.path.dirname(out_fn)):
        os.makedirs(os.path.dirname(out_fn))
    DF.to_csv(out_fn, sep='\t')
    anno = parse_annotations(DF.index)
    anno.to_csv(out_fn + '.annotations', sep='\t', index=False)
