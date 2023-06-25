from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from multiprocessing import Pool
import os

def get_met(sam):
    me = dict()
    for line in open(sam, 'r'):
        if '@' in line:
            continue
        info = line.split('\t')
        qst = info[2]
        if qst == '*':
            continue
        index = info[0].split('*')[4]
        if index not in me:
            me[index] = [qst]
        else:
            me[index].append(qst)
    return me

def judge_met(me, minscore):
    final = dict()
    for k, v in me.items():     # score=1, allnum, allscore         score=1表示被比对到了cs上
        st = {'Alu': [0, 0, 0, 'Alu'], 'L1': [0, 0, 0, 'L1'], 'SVA': [0, 0, 0, 'SVA']}  # count, score
        for i in v:
            if i[1] < minscore:
                continue
            if int(i[1]) == 1:
                st[i[0]][0] += 1
            st[i[0]][1] += 1
            st[i[0]][2] += float(i[1])
        stl = list(st.values())
        stl = sorted(stl, key=lambda x: [x[0], x[1], x[2]])
        if stl[-1][0:3] == [0, 0, 0]:
            final[k] = 'NULL'
        else:
            final[k] = stl[-1][-1]
    return final

def write_vcf(orivcf, output, me):
    vcf = open(orivcf, 'r').readlines()
    for k, v in me.items():
        sv = vcf[int(k)].strip().split('\t')
        sv[2] = '<' + ':'.join(sv[2].split('.')[0:2] + [v]) + '>'
        vcf[int(k)] = '\t'.join(sv) + '\n'
    with open(output, 'w') as f:
        for i in vcf:
            f.write(i)
    f.close()

standard = dict()
concensus = open('concensus.fa').readlines()
for i, line in enumerate(concensus):
    if '>' in line:
        standard[line[1:].strip()] = concensus[i + 1].strip()
aligner = PairwiseAligner(mode='local')
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5
aligner.match_score = 5
aligner.mismatch_score = -4

def comp_seq(taskname, melist, tempdir):
    met = list()
    for me in melist:
        index, seq = me
        most_sim = ''
        max_score = -1e9
        seq = Seq(seq)
        for type, s in standard.items():
            seq2 = Seq(s)
            align = aligner.align(seq2, seq)[0]
            align = str(align)
            match = align.count('|')
            mismatch = align.count('.')
            gap = align.count('-') // 2
            similarity = match / (mismatch + match + gap)
            if similarity > max_score:
                max_score = similarity
                most_sim = type
        met.append([index, most_sim, max_score])
    with open(tempdir + '//' + taskname + '.txt', 'w') as f:
        f.write(str(met))
    f.close()

def multi_comp_seq(args):
    return comp_seq(*args)

def callme(sam, vcf, output, batchsize, thread, tempdir, minscore):
    task = dict()
    st_dic = dict()
    count = 0

    for line in open(sam, 'r'):
        if '@' in line:
            continue
        info = line.strip().split('\t')
        st = info[2]
        index = info[0].split('*')[4]
        if st == '*':
            seq = info[9]
            if count % batchsize == 0:
                task['task' + str(count // batchsize)] = []
            task['task' + str(count // batchsize)].append([index, seq])
            count += 1
        else:
            if index not in st_dic:
                st_dic[index] = [[st, 1]]
            else:
                st_dic[index].append([st, 1])
    '''
    finalme = judge_met(st_dic, minscore)
    write_vcf(vcf, output+'//deepme_align.vcf', finalme)      #只用到了比对
    pools = Pool(processes=thread)
    for key, value in task.items():
        para = [(key, value, tempdir)]
        pools.map_async(multi_comp_seq, para)
    pools.close()
    pools.join()
    '''
    parse_me = list()
    for f in os.listdir(tempdir):
        parse_me += eval(open(tempdir + '//' + f).read())
    for i in parse_me:
        if i[0] not in st_dic:
            st_dic[i[0]] = [[i[1], i[2]]]
        else:
            st_dic[i[0]].append([i[1], i[2]])

    finalme = judge_met(st_dic, minscore)
    write_vcf(vcf, output+'//deepme_sim.vcf', finalme)     # 比对加序列相似性