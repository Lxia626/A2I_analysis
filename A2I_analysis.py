#! /usr/bin/env python3

## the script is used for analyzing A2I editing
import sys, os
from python_analysis import base_pair, custom_sort
from pysam import FastaFile
from collections import OrderedDict

def arg(argv):

    import argparse
    parser = argparse.ArgumentParser(description="Analysis using Python", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help="Annotation file for filtered A2I sites [format: CHROM\tPOS\tREF_NUM,ALT_NUM\tREF\tALT\tRATIO\tSTRAND\tGENE_ID]", type=argparse.FileType('r'), required=True)
    parser.add_argument('--prefix', help="Prefix for output file", type=str)
    parser.add_argument('-R', help="Reference sequences", type=str, required=True)
    parser.add_argument('--depth', help="Threshold for depth [10, --EDITING_TYPE_COUNT]", type=int, default=10)
    parser.add_argument('--mis', help="Depth of mismatch[2, --EDITING_TYPE_COUNT]", type=int, default=2)
    parser.add_argument('--no_annotation', help="No annotation for each site[False, --EDITING_TYPE_COUNT]", action='store_true', default=False)
   
    args = parser.parse_args() 
    sys.stderr.write('args: '+str(args)+'\n')
    main(args.input, FastaFile(args.R), args.prefix, depth=args.depth, mismatch=args.mis, no_annotation=args.no_annotation)

def main(inputFile, ref_in, prefix, depth=10, mismatch=2, no_annotation=False):

    from collections import OrderedDict
    output_count = open(prefix+'.edit_type.txt', 'w')
    filter_output = open(prefix+'.filter.candsites', 'w')
    output_motif = open(prefix+'.motif.seq', 'w')
    output_count.write('TYPE\tCOUNT\tPERCENT\n')
    edit_type_dict = {}

    for i in ['A->T', 'A->C', 'A->G', 'C->G', 'C->A', 'C->T', 'T->A', 'T->G', 'T->C', 'G->C', 'G->T', 'G->A']:
        edit_type_dict[i] = 0

    for line in inputFile:
        if line.startswith('CHROM'): continue
        line_info = line.strip().split('\t')
        #sys.stderr.write(line)
        ref_alt = (int(line_info[2].split(',')[0]), int(line_info[2].split(',')[1]))
        if ref_alt[1] < mismatch: continue
        if sum(ref_alt) < depth: continue
        ref = line_info[3].upper(); alt = line_info[4].upper(); edit_type = ref+'->'+alt
        #ref_alt_ratio = float(line_info[5])
        if not no_annotation: trans_anno = line_info[-3].split(','); strand_set = frozenset(line_info[-2].split(',')); gene_symbol_list = line_info[-1].split(',')
        else: strand_set = frozenset(line_info[6].split(','))

        if len(strand_set) > 1: continue # exclude this site according to Song et al., 2020

        if strand_set == frozenset(['+']):
            edit_type_dict[edit_type] += 1
            if edit_type == 'A->G':
                output_motif.write('>%s\n%s\n' % (line.split('\t')[0]+"_"+line.split('\t')[1], ref_in.fetch(line.split('\t')[0],int(line.split('\t')[1])-3,int(line.split('\t')[1])-1).upper()+'A'+ref_in.fetch(line.split('\t')[0],int(line.split('\t')[1]),int(line.split('\t')[1])+2).upper()))
        else:
            edit_type_dict[base_pair(edit_type[0])+'->'+base_pair(edit_type[-1])] += 1
            if edit_type == 'T->C':
                up_1 = base_pair(ref_in.fetch(line.split('\t')[0],int(line.split('\t')[1])+1,int(line.split('\t')[1])+2)).upper()
                up_2 = base_pair(ref_in.fetch(line.split('\t')[0],int(line.split('\t')[1]),int(line.split('\t')[1])+1)).upper()
                down_1 = base_pair(ref_in.fetch(line.split('\t')[0],int(line.split('\t')[1])-2,int(line.split('\t')[1])-1)).upper()
                down_2 = base_pair(ref_in.fetch(line.split('\t')[0],int(line.split('\t')[1])-3,int(line.split('\t')[1])-2)).upper()
                output_motif.write('>%s\n%s\n' % (line.split('\t')[0]+"_"+line.split('\t')[1], up_1+up_2+'A'+down_1+down_2))
        filter_output.write(line)

    for _k, _v in edit_type_dict.items():
        output_count.write('%s\t%s\t%s\n' % (_k, _v, round(_v/sum(edit_type_dict.values()), 4)))

    output_count.close()
    filter_output.close()

def flank_seq(trans_id_list, chrom, site, strand):

    pass    

if __name__ == '__main__':
    arg(sys.argv)
