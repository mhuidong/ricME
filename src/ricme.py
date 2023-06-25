# -*- coding: UTF-8 -*-
import os.path
import sys
import time
import logging
import argparse
import pysam
from Bio import SeqIO
from multiprocessing import Pool
from ricme_extraction import parse_read
from ricme_combfa import comb_ins_fa
from ricme_callme import *

def parseArgs(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', type=str)
    parser.add_argument('reference', type=str)
    parser.add_argument('vcf', type=str)
    parser.add_argument('rmetlfa', type=str)
    parser.add_argument('--concensus', type=str, default='concensus.fa')
    parser.add_argument('--batchsize', type=int, default=10)
    parser.add_argument('-t', '--thread', type=int, default=16)
    parser.add_argument('-s', '--searchrange', type=int, default=500)
    parser.add_argument('-T', '--tempdir', type=str, default='ricme')
    parser.add_argument('--preset', type=str, default='pacbio')
    parser.add_argument('--minscore', type=float, default='0.7')

    # 多参数并且扩增
    # parser.add_argument('-i','--input',action='append',nargs=2, metavar=('url','name'))
    args = parser.parse_args(argv)
    return args

def parse_base_info(seq):
    info = {'SVLEN': 0, 'END': 0, "SVTYPE": '', "RE": 0}
    for i in seq.split(';'):
        if i.split('=')[0] in ["SVLEN", "END", "RE"]:
            try:
                info[i.split('=')[0]] = abs(int(i.split('=')[1]))
            except:
                pass
        if i.split('=')[0] == "SVTYPE":
            info[i.split('=')[0]] = i.split('=')[1][0:3]
    return info

def is_overlap(b1, b2, c1, c2):
    if min(b2, c2) - max(b1, c1) / (b2-b1) > 0:
        return True
    return False

def find_me(bamfp, reference, taskname, melist, tempdir, searchrange):
    ref = SeqIO.to_dict(SeqIO.parse(reference, "fasta"))
    bam = pysam.AlignmentFile(bamfp)
    with open(tempdir + '/' + taskname + '.fa', 'w') as f:
        for me in melist:
            chrname, pos, svl, svt, line = me
            readname = '*'.join([chrname, str(pos), str(svl), svt, str(line)])
            if svt == 'DEL':
                seq = str(ref[chrname].seq[pos: pos+svl])
                f.write('>'+readname + '\n' + seq + '\n')
            else:
                bp1, bp2 = pos - searchrange, pos + searchrange
                ins = {'INS': dict()}
                for read in bam.fetch(chrname, bp1, bp2):
                    parse_read(read, Chr_name=chrname, SV_size=50, min_mapq=20, max_split_parts=7, min_read_len=500, candidate=ins, min_siglength=10, merge_del_threshold=0, merge_ins_threshold=100, MaxSize=100000)
                index = 0
                for i in ins['INS'][chrname]:
                    if is_overlap(pos, pos+svl, i[0], i[0] + i[1]):
                        index += 1
                        f.write('>' + readname + '*' + str(index) + '\n')
                        f.write(i[3] + '\n')
    f.close()


def multi_find_me(args):
    return find_me(*args)

def main_ctrl(args):
    logging.info('1.get candidate me sequence, input: bam, reference, vcf   output: fa')
    batchsize = args.batchsize
    reference = args.reference
    vcf = args.vcf
    count = 0
    task = dict()
    tempdir = args.tempdir
    if not os.path.exists(tempdir):
        os.mkdir(tempdir)
    mes_tempdir = tempdir + '//' + 'me_seq'
    met_tempdir = tempdir + '//' + 'me_type'
    if not os.path.exists(mes_tempdir):
        os.mkdir(mes_tempdir)
    if not os.path.exists(met_tempdir):
        os.mkdir(met_tempdir)
    ricmefa = tempdir + '//' + 'ricme.fa'
    ricmesam = tempdir + '//' + 'ricme.sam'
    combinefa = tempdir + '//' + 'combine.fa'


    for i, line in enumerate(open(vcf, 'r').readlines()):
        if '#' in line:
            continue
        sv = line.strip().split('\t')
        info = parse_base_info(sv[7])
        svt = info['SVTYPE']
        if svt not in ['INS', 'DEL']:
            continue
        chrname, pos, svl = sv[0], int(sv[1]), info['SVLEN']
        if count % batchsize == 0:
            task['task' + str(count // batchsize)] = []
        task['task' + str(count // batchsize)].append([chrname, pos, svl, svt, i])
        count += 1

    pools = Pool(processes=args.thread)
    for key, value in task.items():
        para = [(args.bam, reference, key, value, mes_tempdir, args.searchrange)]
        pools.map_async(multi_find_me, para)
    pools.close()
    pools.join()

    # combine xxx.fa into one file
    

    with open(ricmefa, 'w') as f:
        filelist = os.listdir(mes_tempdir)
        for file in filelist:
            for line in open(mes_tempdir+'//'+file, 'r'):
                f.write(line)
    f.close()
    
    logging.info('2 Start combine rmetlfa and ricme fa.')
    comb_ins_fa(args.vcf, args.rmetlfa, ricmefa, combinefa)
    
    logging.info('3 Realignment output fa to the reference, input: fa, reference, output: sam')
    os.system('ngmlr -r {} -q {} -o {} -t {} -x {} --subread-length {} --subread-corridor {}'.format(args.concensus, ricmefa, ricmesam, 8, args.preset, 128, 20))
    
    logging.info('3 Judge the me type, input: sam, vcf, output fa')

    callme(ricmesam, vcf, tempdir, batchsize, args.thread, met_tempdir, args.minscore)



def setupLogging(debug=False):
    logLevel = logging.DEBUG if debug else logging.INFO
    logFormat = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig(stream=sys.stderr, level=logLevel, format=logFormat)
    logging.info("Running %s" % " ".join(sys.argv))

def run(argv):
    setupLogging()
    args = parseArgs(argv)
    starttime = time.time()
    main_ctrl(args)
    logging.info("Finished in %0.2f seconds." % (time.time() - starttime))

if __name__ == '__main__':
    run(sys.argv[1:])