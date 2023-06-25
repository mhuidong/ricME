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


def judge_pos(pos1, svl1, pos2, svl2):
    b1, b2 = pos1, pos1 + svl1
    c1, c2 = pos2, pos2 + svl2
    if min(b2, c2) - max(b1, c1) > 0:
        return 1
    elif c1 >= b2:
        return 2
    elif c2 <= b1:
        return 3


def comb_ins_fa(orivcf, rmetlfa, deepmefa, output):
    vcf = list()
    index = dict()
    for i, line in enumerate(open(orivcf, 'r').readlines()):
        if '#' in line:
            continue
        info = line.strip().split('\t')
        otherinfo = parse_base_info(info[7])
        svt = otherinfo['SVTYPE']
        if svt == 'INS':
            chrname, pos, svl = info[0], int(info[1]), otherinfo['SVLEN']
            vcf.append([chrname, pos, svl, str(i)])
    rmetl = list()
    rmetl_rd = dict()
    for line in open(rmetlfa, 'r'):
        if '>' in line and 'INS' in line:
            rmetl.append(line)
    i = j = 0
    while i < len(vcf) and j < len(rmetl):

        chr1, pos1, svl1, id = vcf[i]
        chr2, pos2, svl2 = rmetl[j][1:].strip().split('*')[1:4]
        value = judge_pos(int(pos1), int(svl1), int(pos2), int(svl2))
        if value == 1:
            if id not in index:
                index[id] = 0
            else:
                index[id] += 1
            rmetl_rd[rmetl[j]] = '>' + '*'.join([chr2, pos2, svl2, 'INS', id, str(index[id])]) + '\n'
            j += 1
        elif value == 2:
            i += 1
        elif value == 3:
            j += 1
        # if (i + 1) % 10000 == 0:
        #     print('i', i)
        # if (j + 1) % 10000 == 0:
        #     print('j', j)
    deepme_rd = dict()
    for line in open(deepmefa, 'r'):
        if '>' in line and 'INS' in line:
            chrname, pos, svl, svt, id = line[1:].strip().split('*')[:-1]
            if id not in index:
                index[id] = 0
            else:
                index[id] += 1
            deepme_rd[line] = '>' + '*'.join([chrname, pos, svl, 'INS', id, str(index[id])]) + '\n'
    with open(output, 'w') as f:
        flag = 0
        for line in open(rmetlfa, 'r'):
            if '>' in line:
                if line in rmetl_rd:
                    f.write(rmetl_rd[line])
                    flag = 1
            else:
                if flag == 1:
                    f.write(line)
                    flag = 0
        for line in open(deepmefa, 'r'):
            if '>' in line:
                if line in deepme_rd:
                    f.write(deepme_rd[line])
                    flag = 1
            else:
                if flag == 1:
                    f.write(line)
                    flag = 0
    f.close()
