#! /usr/bin/env python

import sys, os
from collections import OrderedDict
from scipy.stats import binom

def arg(argv):

    import argparse
    parser = argparse.ArgumentParser(description="For correcting the sites obtained according to reads number. (Gabey et al., Nat Commun, 2022)", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-D', '--input_dir', help="Directory for output of editing level files", type=str, required=True)
    parser.add_argument('-o', '--output', help="Output file", type=argparse.FileType('w'), required=True)
    parser.add_argument('-d', '--coverage', help="Coverage threshold [30]. Only record the sites with enough coverage", type=int, default=30)
    parser.add_argument('--seq_err', help="Possibility of sequencing error [0.001]", type=float, default=0.001)
    parser.add_argument('--hetero_err', help="Expected frequency of the alternative allele, for heterozygous SNP [0.5]", type=float, default=0.5)
    parser.add_argument('--p_hetero', help="A-priori probability for heterozygous [1e-4]", type=float, default=10**(-4))
    parser.add_argument('--homo_err', help="Expected frequency of the alternative allele, for homozygous SNP [1.0]", type=float, default=1.0)

    args = parser.parse_args()
    print(args)
    main(args.input_dir, args.output, coverage=args.coverage, seq_err=args.seq_err, hetero_err=args.hetero_err, homo_err=args.homo_err, p_hetero=args.p_hetero)

def main(input_dir, outputFile, coverage=30, seq_err=0.001, hetero_err=0.5, homo_err=1, p_hetero=10**(-4)):
    
    sample_list = os.walk(input_dir)[0][1]
    p_homo = p_hetero ** 2
    for sample in sample_list:
        inputFile = '%s/%s/%s.all.levels.txt' % (input_dir, sample, sample)
        editing_level(inputFile, coverage, seq_err, hetero_err, homo_err, p_hetero, p_homo)
    
def editing_level(inputFile, coverage, seq_err, hetero_err, homo_err, p_hetero, p_homo):

    ed_level_dict = OrderedDict()
    for line in open(inputFile, 'r'):
        if line.startswith('#'): continue
        info = line.strip().split('\t')
        ID = '%s_%s_%s' % (info[0], info[2], info[1])
        total = int(info[3]); mismatch = int(info[4])
        if total < coverage: continue
        seq_err_p_1 = binom.pmf(mismatch, total, seq_err) # due to sequencing error, if REF_num > ALT_num
        seq_err_p_2 = binom.pmf(total-mismatch, total, seq_err) # due to sequencing error, if ALT_num > REF_num
        seq_err_p = binom.pmf(mismatch, total, 1-seq_err) # due to sequencing error, according to Gabey et al., Nat Commun, 2022
        hetero_p = binom.pmf(mismatch, total, hetero_err) * p_hetero # represent a rare heterozygous SNP
        homo_p = binom.pmf(mismatch, total, homo_err) * p_homo # represent a rare homozygous SNP
        # ed_level_dict[ID] = min([seq_err_p, hetero_p, homo_p]) # seleted the most likely explanation of the edited reads
        # ed_level_dict[ID] = max([seq_err_p_1, seq_err_p_2]) # in case the REF is due to sequencing error

    return(ed_level_dict)

if __name__ == '__main__':
    arg(sys.argv)
