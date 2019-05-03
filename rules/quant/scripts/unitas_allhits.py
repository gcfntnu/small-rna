#!/usr/bin env python
"""

Parses the unitas *.hits_per_target.txt files.

Thsi file contains all hits in the cdna + refseq fasta files.

"""


import sys
import os
import glob
import argparse



def samplesheet_ids(fn, sep='\t'):
    sample_ids = []
    with open(fn) as fh:
        txt = fh.read().splitlines()
        header = txt.pop(0).split(sep)
        if not 'Sample_ID' in header:
            raise ValueError('`Sample_ID` column not found in samplesheet')
        for line in txt:
            sample_ids.append(line.split('\t')[0])
        return sample_ids

def argparser():
    parser = argparse.ArgumentParser(description='Aggregate unitas tables')
    parser.add_argument('--sample-sheet', help='Optional sample sheet. Will subset aggregated table if needed', dest='samples')
    parser.add_argument('-o ', '--output', help='Output filename. Will default to stdout.')
    parser.add_argument('filenames', nargs='+')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = argparser()

    import pandas as pd
    df_list = []
    for fn in args.filenames:
        sample_id = os.path.dirname(fn).split(os.path.sep)[-1]
        df = pd.read_csv(fn, sep='\t')
        df['Sample_ID'] = [sample_id] * df.shape[0]
        df_list.append(df)
    DF = pd.concat(df_list, axis=0, join='outer', sort=False)
    if args.output is None:
        out_fn = sys.stdout
    else:
        out_fn = args.output
    DF.fillna(0, inplace=True)
    DF.to_csv(out_fn, sep='\t', index=False)

    # pivot counts
    X = df.pivot(index='Sample_ID', columns='TRANSCRIPT_NAME', values="READ_COUNT")
    X = X.fillna(0)
    
    #class summary
    C = df.pivot_table(index='Sample_ID', columns='TRANSCRIPT_CLASS', values="READ_COUNT")
