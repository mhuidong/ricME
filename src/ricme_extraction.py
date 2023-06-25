import cigar

dic_starnd = {1: '+', 2: '-'}
signal = {1 << 2: 0,
          1 >> 1: 1,
          1 << 4: 2,
          1 << 11: 3,
          1 << 4 | 1 << 11: 4}

def detect_flag(Flag):
    back_sig = signal[Flag] if Flag in signal else 0
    return back_sig


def parse_read(read, Chr_name, SV_size, min_mapq, max_split_parts, min_read_len, candidate, min_siglength,
               merge_del_threshold, merge_ins_threshold, MaxSize):
    if read.query_length < min_read_len:
        return 0
    Combine_sig_in_same_read_ins = list()
    Combine_sig_in_same_read_del = list()
    process_signal = detect_flag(read.flag)
    if read.query_length < 500:
        return 0
    if read.mapq >= min_mapq:
        pos_start = read.reference_start  # 0-based
        pos_end = read.reference_end
        shift_del = 0
        shift_ins = 0
        softclip_left = 0
        softclip_right = 0
        hardclip_left = 0
        hardclip_right = 0
        shift_ins_read = 0
        if read.cigar[0][0] == 4:
            softclip_left = read.cigar[0][1]
        if read.cigar[0][0] == 5:
            hardclip_left = read.cigar[0][1]
        for element in read.cigar:
            if element[0] in [0, 7, 8]:
                shift_del += element[1]
            if element[0] == 2 and element[1] < 50:
                shift_del += element[1]
            if element[0] == 2 and element[1] >= 50:
                Combine_sig_in_same_read_del.append([pos_start + shift_del, element[1]])
                shift_del += element[1]

            # calculate offset of an ins sig in read
            if element[0] != 2:
                shift_ins_read += element[1]

            if element[0] in [0, 2, 7, 8]:
                shift_ins += element[1]
            if element[0] == 1 and element[1] >= 50:
                Combine_sig_in_same_read_ins.append([pos_start + shift_ins, element[1],
                                                     str(read.query_sequence[shift_ins_read - element[
                                                         1] - hardclip_left:shift_ins_read - hardclip_left])])

        if read.cigar[-1][0] == 4:
            softclip_right = read.cigar[-1][1]
        if read.cigar[-1][0] == 5:
            hardclip_right = read.cigar[-1][1]

        if hardclip_left != 0:
            softclip_left = hardclip_left
        if hardclip_right != 0:
            softclip_right = hardclip_right

    generate_combine_sigs(Combine_sig_in_same_read_ins, Chr_name, read.query_name, "INS", candidate,
                          merge_ins_threshold)
    generate_combine_sigs(Combine_sig_in_same_read_del, Chr_name, read.query_name, "DEL", candidate,
                          merge_del_threshold)
    if process_signal == 1 or process_signal == 2:
        Tags = read.get_tags()
        if read.mapq >= min_mapq:
            if process_signal == 1:
                primary_info = [softclip_left, read.query_length - softclip_right, pos_start,
                                pos_end, Chr_name, dic_starnd[process_signal]]
            else:
                primary_info = [softclip_right, read.query_length - softclip_left, pos_start,
                                pos_end, Chr_name, dic_starnd[process_signal]]
        else:
            primary_info = []

        for i in Tags:
            if i[0] == 'SA':
                Supplementary_info = i[1].split(';')[:-1]
                organize_split_signal(Chr_name, primary_info, Supplementary_info, read.query_length,
                                      SV_size, min_mapq, max_split_parts, read.query_name, candidate, MaxSize,
                                      read.query_sequence)

def acquire_clip_pos(deal_cigar):
    seq = list(cigar.Cigar(deal_cigar).items())
    if seq[0][1] == 'S':
        first_pos = seq[0][0]
    else:
        first_pos = 0
    if seq[-1][1] == 'S':
        last_pos = seq[-1][0]
    else:
        last_pos = 0

    bias = 0
    for i in seq:
        if i[1] == 'M' or i[1] == 'D' or i[1] == '=' or i[1] == 'X':
            bias += i[0]
    return [first_pos, last_pos, bias]


def organize_split_signal(chr, primary_info, Supplementary_info, total_L, SV_size,
                          min_mapq, max_split_parts, read_name, candidate, MaxSize, query):
    split_read = list()
    if len(primary_info) > 0:
        split_read.append(primary_info)
        min_mapq = 0
    for i in Supplementary_info:
        seq = i.split(',')
        local_chr = seq[0]
        local_start = int(seq[1])
        local_cigar = seq[3]
        local_strand = seq[2]
        local_mapq = int(seq[4])
        if local_mapq >= min_mapq:
            # if local_mapq >= 0:
            local_set = acquire_clip_pos(local_cigar)
            if local_strand == '+':
                split_read.append([local_set[0], total_L - local_set[1], local_start,
                                   local_start + local_set[2], local_chr, local_strand])
            else:
                try:
                    split_read.append([local_set[1], total_L - local_set[0], local_start,
                                       local_start + local_set[2], local_chr, local_strand])
                except:
                    pass
    if len(split_read) <= max_split_parts or max_split_parts == -1:
        analysis_split_read(split_read, SV_size, total_L, read_name, candidate, MaxSize, query)


def generate_combine_sigs(sigs, Chr_name, read_name, svtype, candidate, merge_dis):
    # for i in sigs:
    # 	print(svtype,i, len(sigs))
    if len(sigs) == 0:
        pass
    elif len(sigs) == 1:
        if Chr_name not in candidate[svtype]:
            candidate[svtype][Chr_name] = list()
        if svtype == 'INS':
            candidate[svtype][Chr_name].append([sigs[0][0],
                                                sigs[0][1],
                                                read_name,
                                                sigs[0][2]])
        else:
            candidate[svtype][Chr_name].append([sigs[0][0],
                                                sigs[0][1],
                                                read_name])
    else:
        temp_sig = sigs[0]
        if svtype == "INS":
            temp_sig += [sigs[0][0]]
            for i in sigs[1:]:
                if i[0] - temp_sig[3] <= merge_dis:
                    temp_sig[1] += i[1]
                    temp_sig[2] += i[2]
                    temp_sig[3] = i[0]
                else:
                    if Chr_name not in candidate[svtype]:
                        candidate[svtype][Chr_name] = list()
                    candidate[svtype][Chr_name].append([temp_sig[0],
                                                        temp_sig[1],
                                                        read_name,
                                                        temp_sig[2]])
                    temp_sig = i
                    temp_sig.append(i[0])
            if Chr_name not in candidate[svtype]:
                candidate[svtype][Chr_name] = list()
            candidate[svtype][Chr_name].append([temp_sig[0],
                                                temp_sig[1],
                                                read_name,
                                                temp_sig[2]])
        else:
            temp_sig += [sum(sigs[0])]
            # merge_dis_bias = max([i[1]] for i in sigs)
            for i in sigs[1:]:
                if i[0] - temp_sig[2] <= merge_dis:
                    temp_sig[1] += i[1]
                    temp_sig[2] = sum(i)
                else:
                    if Chr_name not in candidate[svtype]:
                        candidate[svtype][Chr_name] = list()
                    candidate[svtype][Chr_name].append([temp_sig[0],
                                                        temp_sig[1],
                                                        read_name])
                    temp_sig = i
                    temp_sig.append(i[0])
            if Chr_name not in candidate[svtype]:
                candidate[svtype][Chr_name] = list()
            candidate[svtype][Chr_name].append([temp_sig[0],
                                                temp_sig[1],
                                                read_name])


def analysis_split_read(split_read, SV_size, RLength, read_name, candidate, MaxSize, query):
    '''
    read_start	read_end	ref_start	ref_end	chr	strand
    #0			#1			#2			#3		#4	#5
    '''
    SP_list = sorted(split_read, key=lambda x: x[0])
    # print(read_name)
    # for i in SP_list:
    # 	print(i)
    # detect INS involoved in a translocation
    trigger_INS_TRA = 0
    # Store Strands of INV
    if len(SP_list) == 2:
        ele_1 = SP_list[0]
        ele_2 = SP_list[1]
        if ele_1[4] == ele_2[4]:
            if ele_1[5] == ele_2[5]:
                # dup & ins & del
                a = 0
                if ele_1[5] == '-':
                    ele_1 = [RLength - SP_list[a + 1][1], RLength - SP_list[a + 1][0]] + SP_list[a + 1][2:]
                    ele_2 = [RLength - SP_list[a][1], RLength - SP_list[a][0]] + SP_list[a][2:]
                    query = query[::-1]
                if ele_1[3] - ele_2[2] < SV_size:
                    if ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1] >= SV_size:
                        if ele_2[4] not in candidate["INS"]:
                            candidate["INS"][ele_2[4]] = list()
                        if ele_2[2] - ele_1[3] <= 100 and (
                                ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1] <= MaxSize or MaxSize == -1):
                            candidate["INS"][ele_2[4]].append([(ele_2[2] + ele_1[3]) / 2,
                                                               ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1],
                                                               read_name,
                                                               str(query[
                                                                   ele_1[1] + int((ele_2[2] - ele_1[3]) / 2):ele_2[
                                                                                                                 0] - int(
                                                                       (ele_2[2] - ele_1[3]) / 2)])])
                    if ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3] >= SV_size:
                        if ele_2[4] not in candidate["DEL"]:
                            candidate["DEL"][ele_2[4]] = list()

                        if ele_2[0] - ele_1[1] <= 100 and (
                                ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3] <= MaxSize or MaxSize == -1):
                            candidate["DEL"][ele_2[4]].append([ele_1[3],
                                                               ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3],
                                                               read_name])
        else:
            trigger_INS_TRA = 1
    else:
        # over three splits
        for a in range(len(SP_list[1:-1])):
            ele_1 = SP_list[a]
            ele_2 = SP_list[a + 1]
            ele_3 = SP_list[a + 2]

            if ele_1[4] == ele_2[4]:
                if ele_2[4] == ele_3[4]:
                    if ele_1[5] == ele_3[5] and ele_1[5] == ele_2[5]:
                        # dup & ins & del
                        if ele_1[5] == '-':
                            ele_1 = [RLength - SP_list[a + 2][1], RLength - SP_list[a + 2][0]] + SP_list[a + 2][2:]
                            ele_2 = [RLength - SP_list[a + 1][1], RLength - SP_list[a + 1][0]] + SP_list[a + 1][2:]
                            ele_3 = [RLength - SP_list[a][1], RLength - SP_list[a][0]] + SP_list[a][2:]
                            query = query[::-1]

                        if ele_1[3] - ele_2[2] < SV_size:
                            if ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1] >= SV_size:
                                if ele_2[4] not in candidate["INS"]:
                                    candidate["INS"][ele_2[4]] = list()

                                if ele_2[2] - ele_1[3] <= 100 and (
                                        ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1] <= MaxSize or MaxSize == -1):
                                    if ele_3[2] >= ele_2[3]:
                                        candidate["INS"][ele_2[4]].append([(ele_2[2] + ele_1[3]) / 2,
                                                                           ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1],
                                                                           read_name,
                                                                           str(query[ele_1[1] + int(
                                                                               (ele_2[2] - ele_1[3]) / 2):ele_2[
                                                                                                              0] - int(
                                                                               (ele_2[2] - ele_1[3]) / 2)])])
                            if ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3] >= SV_size:
                                if ele_2[4] not in candidate["DEL"]:
                                    candidate["DEL"][ele_2[4]] = list()

                                if ele_2[0] - ele_1[1] <= 100 and (
                                        ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3] <= MaxSize or MaxSize == -1):
                                    if ele_3[2] >= ele_2[3]:
                                        candidate["DEL"][ele_2[4]].append([ele_1[3],
                                                                           ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3],
                                                                           read_name])

                        if len(SP_list) - 3 == a:
                            ele_1 = ele_2
                            ele_2 = ele_3
                            if ele_1[3] - ele_2[2] < SV_size <= ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1]:
                                if ele_2[4] not in candidate["INS"]:
                                    candidate["INS"][ele_2[4]] = list()

                                if ele_2[2] - ele_1[3] <= 100 and (
                                        ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1] <= MaxSize or MaxSize == -1):
                                    candidate["INS"][ele_2[4]].append([(ele_2[2] + ele_1[3]) / 2,
                                                                       ele_2[0] + ele_1[3] - ele_2[2] - ele_1[1],
                                                                       read_name,
                                                                       str(query[
                                                                           ele_1[1] + int((ele_2[2] - ele_1[3]) / 2):
                                                                           ele_2[0] - int((ele_2[2] - ele_1[3]) / 2)])])

                            if ele_1[3] - ele_2[2] < SV_size <= ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3]:
                                if ele_2[4] not in candidate["DEL"]:
                                    candidate["DEL"][ele_2[4]] = list()

                                if ele_2[0] - ele_1[1] <= 100 and (
                                        ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3] <= MaxSize or MaxSize == -1):
                                    candidate["DEL"][ele_2[4]].append([ele_1[3],
                                                                       ele_2[2] - ele_2[0] + ele_1[1] - ele_1[3],
                                                                       read_name])

    if len(SP_list) >= 3 and trigger_INS_TRA == 1:
        if SP_list[0][4] == SP_list[-1][4]:
            # print(SP_list[0])
            # print(SP_list[-1])
            if SP_list[0][5] != SP_list[-1][5]:
                pass
            else:
                if SP_list[0][5] == '+':
                    ele_1 = SP_list[0]
                    ele_2 = SP_list[-1]
                else:
                    ele_1 = [RLength - SP_list[-1][1], RLength - SP_list[-1][0]] + SP_list[-1][2:]
                    ele_2 = [RLength - SP_list[0][1], RLength - SP_list[0][0]] + SP_list[0][2:]
                    query = query[::-1]
                # print(ele_1)
                # print(ele_2)
                dis_ref = ele_2[2] - ele_1[3]
                dis_read = ele_2[0] - ele_1[1]
                if dis_ref < 100 and dis_read - dis_ref >= SV_size and (dis_read - dis_ref <= MaxSize or MaxSize == -1):
                    # print(min(ele_2[2], ele_1[3]), dis_read - dis_ref, read_name)
                    if ele_1[4] not in candidate['INS']:
                        candidate['INS'][ele_1[4]] = list()
                    candidate["INS"][ele_2[4]].append([min(ele_2[2],
                                                           ele_1[3]),
                                                       dis_read - dis_ref,
                                                       read_name,
                                                       str(query[
                                                           ele_1[1] + int(dis_ref / 2):ele_2[0] - int(dis_ref / 2)])])