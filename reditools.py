#! /usr/bin/env python

import sys
from math import ceil
from python_analysis import gtf_chrom_trans

def arg(argv):

    import argparse
    parser = argparse.ArgumentParser(description="Analyze REDItools results", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help="REDItools output", type=argparse.FileType('r'), required=True)
    parser.add_argument('--Alu_bed', help="BED file for ALU region [None]", type=str, default=None) # /public/home/Songlab/songyl/database/hg19_Alu.bed
    parser.add_argument('-gtf_table', help="Annotation table for current genome ['/public/home/Songlab/xial/Work/database/hg19/hg19_ENSEMBL_all_RNAs.table']", type=str, default='/public/home/Songlab/xial/Work/database/hg19/hg19_ENSEMBL_all_RNAs.table')
    parser.add_argument('-d', '--depth', help="Minimum depth of each candidate site [10]", type=int, default=10)
    parser.add_argument('-r', '--region', help="Alu and non-Alu region ['Both'(default), 'Alu', 'non-Alu']", type=str, default='Both')
    parser.add_argument('-p',  help="Number of multi-processes [1]", type=int, default=1)
    parser.add_argument('-o', '--output', help="Output file", type=str, required=True)

    args = parser.parse_args()
    sys.stderr.write('args: '+str(args)+'\n')
    if args.region != 'Both' and not args.Alu_bed: sys.stderr.write('Please specify ALU_bed.\n'); sys.exit()
    if args.region != 'Both': 
        Alu_region = {}
        for line in open(args.Alu_bed):
            Alu_region.setdefault(line.split('\t')[0], []).append((int(line.split('\t')[1]), int(line.strip().split('\t')[2])))
    else: Alu_region = None
    gtf_info = gtf_chrom_trans(args.gtf_table)
   
    inputFile = args.input.readlines() 
    sys.stderr.write('Start candidate site filtration...\n')
    if args.p == 1:
        sys.stderr.write('...Single process: Candidate site number: %s\n' % (len(inputFile)))
        out_dict = main(inputFile, depth=args.depth, Alu_region=Alu_region, region=args.region, gtf_info=gtf_info)
    else:
        from multiprocessing import Pool
        pool = Pool(processes=args.p); results = []
        separate = int(ceil(len(inputFile)/args.p))
        for i in range(args.p):
            inputSite = inputFile[separate*i:separate*(i+1)]
            sys.stderr.write('...Simultaneously: subprocess %s start running. Candidate site number: %s\n' % (i, len(inputSite)))
            result = pool.apply_async(main, (inputSite, args.depth, Alu_region, args.region, i, gtf_info, ))
            results.append(result)
        pool.close()
        pool.join()
        out_dict = {}
        for _rlt in results:
            rlt = _rlt.get()
            for _k, _v in rlt.items():
                out_dict.setdefault(_k, []).extend(_v)

    sys.stderr.write('Start output final candidate sites.....\n')
    output_file(out_dict, args.output)
    sys.stderr.write('Finish running successfully!\n')
    

def main(inFile, depth=10, Alu_region=None, region='Both', subpro='Single', gtf_info=None):
    
    out_dict = {}; a = 0
    for line in inFile:
        a += 1
        if a % 100000 == 0: sys.stderr.write('%s records has analyzed for subprocess %s.\n' % (a, subpro))
        if line.startswith('Region') or line.split('\t')[7] == '-': continue
        info = line.strip().split('\t')
        if info[0] == 'chrM': continue
        pos = int(info[1]); chrom = info[0]; ref = info[2]
        if len(info[7]) > 2: continue
        TYPE = info[7][0]+'->'+info[7][1]; alt = info[7][1]; ratio=float(info[8]); strand = info[3]
        base_index = ['A', 'C', 'G', 'T']
        base_count_info = info[6].replace('[', '').replace(']', '').split(', ')
        base_count = dict([(base_index[x], int(base_count_info[x])) for x in range(len(base_count_info))])
        if sum(base_count.values()) < depth: continue
        record = False
        if region != 'Both':
            if chrom not in Alu_region: continue
            if [x for x in Alu_region[chrom] if x[0] < pos <= x[1]] and region == 'Alu': record = True # in Alu region
            elif not [x for x in Alu_region[chrom] if x[0] < pos <= x[1]] and region == 'non-Alu': record = True # in non-Alu region
        else:
            record = True
        if record:
            _g = gtf_info_collect(chrom, pos, gtf_info)
            ratio = round(base_count[alt]/sum(base_count.values()), 4)
            out_dict.setdefault(chrom, []).append('\t'.join([str(pos), ref, alt, str(base_count[ref])+','+str(base_count[alt]), str(ratio), _g]))
    
    return(out_dict)    

def output_file(out_dict, out_file):

    outFile = open(out_file, 'w')
    outFile.write('CHROM\tPOS\tREF\tALT\tDEPTH\tRATIO\tTRANS_ANNO\tSTRAND\tGENE_SYMBOL\n')
    summary_out = open(out_file+'.summary', 'w')
    summary_out.write('TYPE\tCOUNT\tPERCENT\n')
    record_dict = {}

    for _k in sorted(out_dict.keys()):
        for _r in out_dict[_k]:
            TYPE = _r.split('\t')[1]
            if TYPE in record_dict: record_dict[TYPE] += 1
            else: record_dict[TYPE] = 1
            outFile.write('%s\t%s\n' % (_k, _r))
    outFile.close()

    for _k in record_dict:
        summary_out.write('%s\t%s\t%s\n' % (_k, record_dict[_k], round(record_dict[_k]/sum(record_dict.values()),4)))
    summary_out.close()

def gtf_info_collect(chrom, pos, gtf_info):

    chrom = chrom.replace('chr', '')
    cand_info = [x for x in gtf_info[chrom] if x[0][0] < pos <= x[0][1]]
    trans_id = ','.join([x[5] for x in cand_info])
    strand = ','.join([x[3] for x in cand_info])
    gene_symbol = ','.join(x[4][1] for x in cand_info)

    # in case no transcript has been annotated
    if trans_id == '': return('\t'.join(['intergenic', 'intergenic', 'intergenic']))
    else: return('\t'.join([trans_id, strand, gene_symbol]))

if __name__ == '__main__':
    arg(sys.argv)
