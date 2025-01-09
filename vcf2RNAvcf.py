#! /usr/bin/env python

import sys
from pysam import VariantFile, FastaFile
# from python_analysis import gtf_table_trans
from itertools import starmap
from time import time
from scipy.stats import binom

def arg(argv):

    import argparse
    parser = argparse.ArgumentParser(description="Discard SNP sites according to different database", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-v', help="BCFtools-indexed VCF file", type=str, required=True)
    parser.add_argument('-db', help="Candidate SNP databases, separated by ','", type=str, default=None) #! [ can be filtered previously using /public/usr/bin/bcftools filter in count_dbSNP_depth function in script Nanopore_A2I_editing_pipe.sh]
    parser.add_argument('-sr', '--simple_rep', help="Collection of simple repeats [not need if region=='ALU']", type=str, default=None)
    parser.add_argument('-R', '--ref', help="Reference genome [not need if region=='ALU']", type=str, default=None)
    parser.add_argument('-d', '--depth', help="Minimum depth of SNPs [10]", type=int, default=10) #! [can be filtered previously using /public/usr/bin/bcftools filter in count_dbSNP_depth function in script Nanopore_A2I_editing_pipe.sh]
    parser.add_argument('-m', '--mismatch', help="Minimum number of reads for alts [2]", type=int, default=2) 
    parser.add_argument('-gtf', help="GTF table file ['/public/home/Songlab/songyl/database/genome/human/GRCh37_release75/hg19_ENSEMBL_all_RNAs.table']", type=str, default='/public/home/Songlab/songyl/database/genome/human/GRCh37_release75/hg19_ENSEMBL_all_RNAs.table')
    parser.add_argument('--splice_nt', help="Distance from intron/exon boundary [not need if region=='ALU'] [4]", type=int, default=4)
    parser.add_argument('--buffer', help="Length for flanking sequence of editing candidate site for removing homopolymers [not need if region=='ALU'] [4]", type=int, default=4)
    parser.add_argument('-o', help="Output file", type=argparse.FileType('w'), required=True)
    parser.add_argument('-r', '--region', help="Region to analyse {'ALU', 'NONALU'}", default='ALU', type=str)
    parser.add_argument('--seq_p', help="Possibility for sequencing error", default=0.05, type=float)
    # parser.add_argument('--strand_correct', help="Corrected strand for each reads using 'correct_read_strand' from l-GIREMI", default=None, type=str)

    args = parser.parse_args()
    sys.stderr.write(str(args)+'\n')
    gtf_info = gtf_table_trans(args.gtf)

    if args.db:
        START_1 = time()
        snpdb = {}
        for FILE in args.db.split(','):
            for line in open(FILE, 'r'):
                if line.startswith('chr'): snpdb.setdefault(line.split('\t')[0], []).append(int(line.strip().split('\t')[1]))
                else: snpdb.setdefault('chr'+line.split('\t')[0], []).append(int(line.strip().split('\t')[1]))
        END_1 = time()
        sys.stderr.write('SNP databases have been loaded. Running time: %s s.\n' % (END_1-START_1))
    else: snpdb = {}
    if args.region != 'ALU':
        if not args.simple_rep or not args.ref:
            sys.stderr.write('Please specify -sr and -R when region !="ALU".\n'); sys.exit()
        ref_genome = FastaFile(args.ref)
        START_2 = time()
        simpleRep_dict = {}
        for line in open(args.simple_rep, 'r'):
            if line.split('\t')[5].startswith('chr'): simpleRep_dict.setdefault(line.split('\t')[5], []).append([int(line.split('\t')[6]), int(line.split('\t')[7])])
            else: simpleRep_dict.setdefault('chr'+line.split('\t')[5], []).append([int(line.split('\t')[6]), int(line.split('\t')[7])])
        END_2 = time()
        sys.stderr.write('Information of simple repeats has been loaded. Running time: %s s.\n' % (END_2-START_2))
        START_3 = time()
        sa_dict = {} # dictionary for splicing artefacts and transcript region
        for _t in gtf_info.keys():
            exon_s = [x[0]+1 for x in gtf_info[_t][2][1:-1]]; exon_e = [x[1] for x in gtf_info[_t][2][:-2]]
            sa_region = list(zip(exon_e, exon_s))
            for _sr in range(len(sa_region)):
                if _t in sa_dict: 
                    if [sa_region[_sr][0]-args.splice_nt, sa_region[_sr][0]+args.splice_nt] in sa_dict[_t]: continue
                    sa_dict[_t].append([sa_region[_sr][0]-args.splice_nt, sa_region[_sr][0]+args.splice_nt])
                else: sa_dict[_t] = [[sa_region[_sr][0]-args.splice_nt, sa_region[_sr][0]+args.splice_nt]]
        END_3 = time()
        sys.stderr.write('GTF information has been loaded. Running time: %s s.\n' % (END_3-START_3))
    else: 
        simpleRep_dict = None; sa_dict = None; ref_genome = None

    main(args.v, snpdb, simpleRep_dict, sa_dict, gtf_info, args.o, ref_genome, depth=args.depth, homo_buffer=args.buffer, region=args.region, mismatch=args.mismatch, seq_error_p=args.seq_p)

def gtf_table_trans(gtf_file):

    anno_dict = {}
    for line in open(gtf_file, 'r'):
        line_info = line.strip().split('\t')
        trans_strand = line_info[3]
        gene_region = [int(line_info[4]), int(line_info[5])]
        CDS_region = [int(line_info[6]), int(line_info[7])]
        exon_region = list(zip([int(x) for x in line_info[9].split(',')[:-1]], [int(x) for x in line_info[10].split(',')[:-1]]))
        gene_id = [line_info[-3], line_info[-2]] # (gene_id, gene_symbol)
        trans_type = line_info[-1]
        chrom = line_info[2]
        length = sum(starmap(lambda x,y:y-x, exon_region))
        anno_dict[line_info[1]] = [gene_region, CDS_region, exon_region, trans_strand, gene_id, chrom, trans_type, length]

    return(anno_dict)

def main(vcfFile, snpdb, simpleRep_dict, sa_dict, gtf_info, outFile, ref_genome, depth=10, homo_buffer=4, region='ALU', mismatch=2, seq_error_p=0.05):

    vcf_in = VariantFile(vcfFile, 'r')

    a = 0; START_x = time(); filtered_list = {}
    for rec in vcf_in.fetch():

        ##### basic filtration
        a += 1
        if a % 10000 == 0: 
            END_x = time()
            sys.stderr.write('%s records have been analyzed. Running time: %s s.\n' % (a, END_x-START_x))
            START_x = time()
        ref_for, ref_rev, alt_for, alt_rev = rec.info['DP4']
        # if alt_for > 0 and alt_rev > 0: continue # for editing site, the ALT should only on the same strand [modification 231002: maybe SNP and editing co-occured?]
        # if len([x for x in rec.info['DP4'] if x > 0]) > 2: continue # shoud on single strand [modification 231002: maybe ref can be on both strand?]
        if snpdb:
            if rec.chrom not in snpdb: continue # [modification 231002: SNPdb may have no such chromosome?]
        if region != 'ALU':
            if rec.chrom not in simpleRep_dict: continue
        if len(rec.alts) > 1: sys.stderr.write('%s:%s has more than one ALTs\n' % (rec.chrom, rec.pos)); continue # should be only one type of mutation
        ##### detailed filtration
        # not in SNP database
        if snpdb:
            if rec.chrom in snpdb:
                if rec.pos in snpdb[rec.chrom]: filtered_list.setdefault('snpdb', []).append(rec.pos); continue
        # enough coverage
        tmp_trans_id = []
        for x in gtf_info.keys(): # for each transcript
            if 'chr'+gtf_info[x][5] != rec.chrom and gtf_info[x][5] != rec.chrom: continue
            if gtf_info[x][0][0] < rec.pos <= gtf_info[x][0][1]: tmp_trans_id.append(x)
        # tmp_trans_id = [x for x in gtf_info.keys() if gtf_info[x][0][0] < rec.pos <= gtf_info[x][0][1]]
        if not tmp_trans_id: trans_id = 'NA' # in intergenic region
        else: trans_id = sorted(tmp_trans_id, key=lambda x:gtf_info[x][-1], reverse=True)[0] # the longest transcript covering target site
        if trans_id == 'NA': # in gene region
            if sum(rec.info['DP4']) < depth: filtered_list.setdefault('depth', []).append(rec.pos); continue
        else:
            strand_trans = gtf_info[trans_id][3]
            if strand_trans == '+':
                if ref_for + alt_for < depth: filtered_list.setdefault('depth', []).append(rec.pos); continue
            else:
                if ref_rev + alt_rev < depth: filtered_list.setdefault('depth', []).append(rec.pos); continue
        # sequence error, according to Gabay et al., Nat Commun, 2022
        if binom.pmf(k=alt_for+alt_rev, n=sum(rec.info['DP4']), p=seq_error_p) > 0.05: filtered_list.setdefault('seq_error', []).append(rec.pos); continue
        if region != 'ALU':
            # not in simple repeat regions
            if [x for x in simpleRep_dict[rec.chrom] if x[0] <= rec.pos <= x[1]]: continue
            # not containing homopolymer
            if homo_detect(ref_genome, homo_buffer, rec.chrom, rec.pos): continue
            # BLAT analysis

        ##### transcript filtration
        if region != 'ALU':
            # not in splicing artefects
            if trans_id in sa_dict: # otherwise maybe only one exon for target transcript or target site is in intergenic region
                if [x for x in sa_dict[trans_id] if x[0] <= rec.pos <= x[1]]: continue 
        if trans_id != 'NA': # in gene region
            strand_trans = gtf_info[trans_id][3]
            if strand_trans == '+':
                # sys.stderr.write('ref_for: '+str(ref_for)+'\n')    
                if alt_for < mismatch: filtered_list.setdefault('mismatch', []).append(rec.pos); continue
                if ref_for+alt_for == 0: filtered_list.setdefault('strand', []).append(rec.pos); continue # not on target strand
                if ref_for == 0: outFile.write('%s\t%s\t%s\t%s\t%s\t%s\t+\t%s\n' % (rec.chrom, rec.pos, str(ref_for)+','+str(alt_for), rec.ref, rec.alts[0], 1, gtf_info[trans_id][4][0]))
                else: outFile.write('%s\t%s\t%s\t%s\t%s\t%s\t+\t%s\n' % (rec.chrom, rec.pos, str(ref_for)+','+str(alt_for), rec.ref, rec.alts[0], round(alt_for/(ref_for+alt_for), 4), gtf_info[trans_id][4][0]))
                # if ref_for == 0: outFile.write('%s\t%s\t%s\t%s\t%s\t%s\t+\n' % (rec.chrom, rec.pos+1, rec.ref, rec.alts[0], str(ref_for)+','+str(alt_for), 1))
                # else: outFile.write('%s\t%s\t%s\t%s\t%s\t%s\t+\n' % (rec.chrom, rec.pos+1, rec.ref, rec.alts[0], str(ref_for)+','+str(alt_for), round(alt_for/(ref_for+alt_for), 4)))
            else:
                # sys.stderr.write('ref_rev: '+str(ref_rev)+'\n')
                if alt_rev < mismatch: filtered_list.setdefault('mismatch', []).append(rec.pos); continue
                if ref_rev+alt_rev == 0: filtered_list.setdefault('strand', []).append(rec.pos); continue # not on target strand
                if ref_rev == 0: outFile.write('%s\t%s\t%s\t%s\t%s\t%s\t-\t%s\n' % (rec.chrom, rec.pos, str(ref_rev)+','+str(alt_rev), rec.ref, rec.alts[0], 1, gtf_info[trans_id][4][0]))
                else: outFile.write('%s\t%s\t%s\t%s\t%s\t%s\t-\t%s\n' % (rec.chrom, rec.pos, str(ref_rev)+','+str(alt_rev), rec.ref, rec.alts[0], round(alt_rev/(ref_rev+alt_rev), 4), gtf_info[trans_id][4][0]))
                # if ref_rev == 0: outFile.write('%s\t%s\t%s\t%s\t%s\t%s\t-\n' % (rec.chrom, rec.pos+1, rec.ref, rec.alts[0], str(ref_rev)+','+str(alt_rev), 1))
                # outFile.write('%s\t%s\t%s\t%s\t%s\t%s\t-\n' % (rec.chrom, rec.pos+1, rec.ref, rec.alts[0], str(ref_rev)+','+str(alt_rev), round(alt_rev/(ref_rev+alt_rev), 4)))
        else: # in intergenic region and regard as strand_trans == '+'
            if alt_for < mismatch: filtered_list.setdefault('mismatch', []).append(rec.pos); continue
            if ref_for+alt_for == 0: filtered_list.setdefault('strand', []).append(rec.pos); continue # not on target strand
            if ref_for == 0: outFile.write('%s\t%s\t%s\t%s\t%s\t%s\t+\tINTERGENIC\n' % (rec.chrom, rec.pos, str(ref_for+ref_rev)+','+str(alt_for+alt_rev), rec.ref, rec.alts[0], 1))
            else: outFile.write('%s\t%s\t%s\t%s\t%s\t%s\t+\tINTERGENIC\n' % (rec.chrom, rec.pos, str(ref_for+ref_rev)+','+str(alt_for+alt_rev), rec.ref, rec.alts[0], round((alt_for+alt_rev)/sum(rec.info['DP4']), 4)))
            # if ref_for == 0: outFile.write('%s\t%s\t%s\t%s\t%s\t%s\t+\n' % (rec.chrom, rec.pos+1, rec.ref, rec.alts[0], str(ref_for+ref_rev)+','+str(alt_for+alt_rev), 1))
            # else: outFile.write('%s\t%s\t%s\t%s\t%s\t%s\t+\n' % (rec.chrom, rec.pos+1, rec.ref, rec.alts[0], str(ref_for+ref_rev)+','+str(alt_for+alt_rev), round((alt_for+alt_rev)/sum(rec.info['DP4']), 4)))
    
    for _k in filtered_list:
        sys.stderr.write('Sites filtered by %s: %s\n' % (_k, len(filtered_list[_k])))

    outFile.close()

def homo_detect(ref_genome, homo_buffer, chrom, pos):
    ### ATTENTION: for this part, further check is needed

    homo_len = homo_buffer
    seq = ref_genome.fetch(chrom, pos-homo_buffer, pos+homo_buffer+1)
    for x in range(len(seq)-homo_len+1):
        if len(set(list(seq[x:x+homo_len]))) == 1: return(True) # only one type of base, so in homopolymer
    return(False)

def sequence_error(ref_num, alt_num, p_value):

    seq_err_p = binom.cdf(k=alt_num, n=ref_num+alt_num, p=p_value)
    hetero_p = binom.cdf(k=alt_num, n=ref_num+alt_num, p=0.5) * 10^-4
    homo_p = binom.cdf(k=alt_num, n=ref_num+alt_num, p=1) * 10^-8
    p_list = [('seq_err_p', seq_err_p), ('hetero_p', hetero_p), ('homo_p', homo_p)]
    return(sorted(p_list, key=lambda x:x[1]))

if __name__ == '__main__':
    arg(sys.argv)
