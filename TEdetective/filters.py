
import sys, os

def exec_filter(args):
    #
    log_FH=open(args.log_file, 'w')
    preprocess_dir_realpath = os.path.realpath(args.preprocess_dir)
    log_FH.write('preprocessing/intermediate file directory: '+ preprocess_dir_realpath +'\n')
    #
    dir_path = os.getcwd()
    log_FH.write('working directory: '+ dir_path +'\n')
    #
    fofn_bed_file = os.path.realpath(args.fofn_bed)
    log_FH.write('fofn_bed file: '+ fofn_bed_file +'\n')
    #
    input_file_name = args.ofa_inp
    log_FH.write('input file name: '+ input_file_name +'\n')
    #
    insert_size = args.isz_inp 
    read_length = args.rdl_inp
    qual_lim = args.qlm_inp
    tcr = args.tcr_inp
    trd = args.trd_inp    
    rp = args.rp_inp
#    trd = tcr + tdr

    #Remove reads falling into existing TE regions
    
    chr_track = 'chr00'
    dict_bed = {}
    with open(fofn_bed_file, 'r') as bed_file:
        bed_file_lines = bed_file.readlines()
    bed_file.close()

    for line in bed_file_lines:
        if line.split()[0] != chr_track:
            line_num = 0
            chr_track = line.split()[0]
            dict_bed[chr_track] = []
            dict_bed[chr_track+'_idx'] = {}
            pos_track = int((int(line.split()[1]))/10000)
            dict_bed[chr_track+'_idx'][pos_track] = 1 #line_num to handle -1 for very first line
            dict_bed[chr_track].append([int(line.split()[1]), int(line.split()[2]), line.split()[3]])
            continue    
        dict_bed[chr_track].append([int(line.split()[1]), int(line.split()[2]), line.split()[3]])
        line_num += 1
        if pos_track != int((int(line.split()[1]))/10000):
            for i in range(pos_track, int((int(line.split()[1]))/10000)):
                if pos_track + 1 == int((int(line.split()[1]))/10000):
                    dict_bed[chr_track+'_idx'][int((int(line.split()[1]))/10000)] = line_num
                    pos_track = int((int(line.split()[1]))/10000)
                if pos_track + 1 < int((int(line.split()[1]))/10000):
                    dict_bed[chr_track+'_idx'][pos_track+1] = dict_bed[chr_track+'_idx'][pos_track]
                    pos_track = pos_track + 1
    del bed_file_lines[:]
##
    left_clipped_rd, right_clipped_rd, left_discord_rd, right_discord_rd, test_class_score = 0, 0, 0, 0, 0

    with open(input_file_name, 'r') as input_file_file:
        input_file_lines = input_file_file.readlines()
    input_file_file.close()

    for line in input_file_lines:
        te_loc = 'null'

        if line.startswith('#'):
            continue
        try:
            for info in dict_bed[line.split()[1]][(dict_bed[line.split()[1]+'_idx'].get(int((int(line.split()[3]))/10000),1))-1:]:
                if int(line.split()[3]) in range( int(info[0]), int(info[1])+1 ): #and int(line.split()[4]) in range (int(info[0]), int(info[1])+1):
                    te_loc = info[2]
                    break
#                elif int(line.split()[3]) in range( int(info[0])-100, int(info[1])+101 ):
#                    te_loc = info[2]+'_100'
#                    break
#                elif int(line.split()[3]) < int(info[0]):  ### <<---- Fix this
#                    te_loc = 'noTE'
#                    break
        except (ValueError, KeyError):
            continue

        if te_loc == 'null':
            try:
                for info in dict_bed[line.split()[1]][(dict_bed[line.split()[1]+'_idx'].get(int((int(line.split()[3]))/10000),1))-1:]:
                    if int(line.split()[3]) in range( int(info[0])-100, int(info[1])+101 ):
                        te_loc = info[2]+'_100'
                        break
                    elif int(line.split()[3]) < int(info[0]):
                        break
            except (ValueError, KeyError):
                continue

        if te_loc == 'null':
            te_loc = 'noTE'

        words = line.split()

#        samfile = pysam.AlignmentFile(bam_full, 'rb')
#        cnt_rd = 0
#        for read in samfile.fetch(words[1], int(words[3])-insert_size, int(words[3])+insert_size):
#            if read.cigarstring != None:
#                if read.get_reference_positions()[0] >= int(words[3])-insert_size and read.get_reference_positions()[-1] <= int(words[3])+insert_size:
#                    cnt_rd += 1
#        samfile.close()

        test_class = words[0]

        if (words[7] == test_class and float(words[8]) > qual_lim) or words[7] == 'NA':    #left clipped read
            left_clipped_rd = int(words[9])
#        if words[7] == test_class or words[7] == 'NA':
#            test_class_score += 1
        if (words[11] == test_class and float(words[12]) > qual_lim) or words[11] == 'NA':    # right clipped reads
            right_clipped_rd = int(words[13])
        # Left poly-A/T
        p_pat_num = int(words[19]) 
        # Right poly-A/T
        n_pat_num = int(words[20])
#        if words[11] == test_class or words[11] == 'NA':
#            test_class_score += 1
        if (words[21] == test_class and float(words[22]) > qual_lim) or words[21] == 'NA':    # left discordant reads
            left_discord_rd = int(words[23])
#        if words[19] == test_class or words[19] == 'NA':
#            test_class_score += 1
        if (words[25] == test_class and float(words[26]) > qual_lim) or words[25] == 'NA':    # right discordant reads
            right_discord_rd = int(words[27])
#        if words[23] == test_class or words[23] == 'NA':
#            test_class_score += 1

        total_clipped_rd = left_clipped_rd + right_clipped_rd
        total_clipped_rd_wpat = total_clipped_rd + p_pat_num + n_pat_num
        #
        total_discord_rd = left_discord_rd + right_discord_rd

        total_rd_left = left_clipped_rd + left_discord_rd
        total_rd_right = right_clipped_rd + right_discord_rd

        filter_result = 'FAIL'

#        if total_clipped_rd >= 3 or ( (total_clipped_rd >= 1) and ( (total_clipped_rd_wpat+total_discord_rd) >= 5) ) or ( total_discord_rd >= 10 ): # L1base Sim
#        if total_clipped_rd >= 5 or ( total_clipped_rd >=3 and (total_clipped_rd_wpat + total_discord_rd >= 7)) or total_discord_rd >= 10: # and test_class_score == 4 \ CEU
#        if total_clipped_rd >= 5 or ( total_clipped_rd >=3 and (total_clipped_rd + total_discord_rd >= 10)) or total_discord_rd >= 15: # and test_class_score == 4 \ ALU
#        if total_clipped_rd >= 5 or ( total_clipped_rd >=3 and (total_clipped_rd + total_discord_rd >= 10)) or total_discord_rd >= 25: # and test_class_score == 4 \ BL6NJ
#        if total_clipped_rd >= 5 or ( total_clipped_rd >=3 and (total_clipped_rd + total_discord_rd >= 10)) or total_discord_rd >= 15: # and test_class_score == 4 \ LTR
#        if total_clipped_rd >= 3 or ( total_clipped_rd >=2 and (total_clipped_rd + total_discord_rd >= 5)) or total_discord_rd >= 10: # and test_class_score == 4 \ LTR-SUB
#        if total_clipped_rd >= 3 or ( total_clipped_rd >=2 and (total_clipped_rd + total_discord_rd >= 5)): # or total_discord_rd >= 25: # and test_class_score == 4 \ Ecat11
#            and ((total_clipped_rd + total_discord_rd)*100/cnt_rd >= rp): # and total_rd_left > 0 and total_rd_right > 0:
        if total_clipped_rd >= 3 or ( (total_clipped_rd >= 1) and ( (total_clipped_rd_wpat+total_discord_rd) >= 5) ):
            filter_result = 'PASS'
        elif total_discord_rd >= 10:
            filter_result = 'PASS_D'

        log_FH.write(filter_result+'\t'+te_loc+'\t'+line)

        left_clipped_rd, right_clipped_rd, left_discord_rd, right_discord_rd, test_class_score = 0, 0, 0, 0, 0
    log_FH.close()


def exec_filter_p(args):
    #
    log_FH=open(args.log_file, 'w')
    preprocess_dir_realpath = os.path.realpath(args.preprocess_dir)
    log_FH.write('preprocessing/intermediate file directory: '+ preprocess_dir_realpath +'\n')
    #
    dir_path = os.getcwd()
    log_FH.write('working directory: '+ dir_path +'\n')
    #
    input_file_name = args.ofa_inp
    log_FH.write('input file: '+ input_file_name +'\n')
    #
    fofn_bed_file = os.path.realpath(args.fofn_bed)
    log_FH.write('fofn_bed file: '+ fofn_bed_file +'\n')
    #
    insert_size = args.isz_inp
    read_length = args.rdl_inp
    qual_lim = args.qlm_inp
    tcr = args.tcr_inp
    trd = args.trd_inp
    rp = args.rp_inp

    #Remove reads falling into existing TE regions

    chr_track = 'chr00'
    dict_bed = {}
    with open(fofn_bed_file, 'r') as bed_file:
        bed_file_lines = bed_file.readlines()
    bed_file.close()

    for line in bed_file_lines:
        if line.split()[0] != chr_track:
            line_num = 0
            chr_track = line.split()[0]
            dict_bed[chr_track] = []
            dict_bed[chr_track+'_idx'] = {}
            pos_track = int((int(line.split()[1]))/10000)
            dict_bed[chr_track+'_idx'][pos_track] = 1 #line_num to handle -1 for very first line
            dict_bed[chr_track].append([int(line.split()[1]), int(line.split()[2]), line.split()[3]])
            continue
        dict_bed[chr_track].append([int(line.split()[1]), int(line.split()[2]), line.split()[3]])
        line_num += 1
        if pos_track != int((int(line.split()[1]))/10000):
            for i in range(pos_track, int((int(line.split()[1]))/10000)):
                if pos_track + 1 == int((int(line.split()[1]))/10000):
                    dict_bed[chr_track+'_idx'][int((int(line.split()[1]))/10000)] = line_num
                    pos_track = int((int(line.split()[1]))/10000)
                if pos_track + 1 < int((int(line.split()[1]))/10000):
                    dict_bed[chr_track+'_idx'][pos_track+1] = dict_bed[chr_track+'_idx'][pos_track]
                    pos_track = pos_track + 1
    del bed_file_lines[:]
##
    left_clipped_rd, right_clipped_rd, left_discord_rd, right_discord_rd, test_class_score = 0, 0, 0, 0, 0

    with open(input_file_name, 'r') as input_file_file:
        input_file_lines = input_file_file.readlines()
    input_file_file.close()

    for line in input_file_lines:

        te_loc = 'null'

        if line.startswith('#'):
            continue
#
        try:
            for info in dict_bed[line.split()[1]][(dict_bed[line.split()[1]+'_idx'].get(int((int(line.split()[3]))/10000),1))-1:]:
                if int(line.split()[3]) in range( int(info[0]), int(info[1])+1 ): #and int(line.split()[4]) in range (int(info[0]), int(info[1])+1):
                    te_loc = info[2]
                    break
#                elif int(line.split()[3]) in range( int(info[0])-100, int(info[1])+101 ):
#                    te_loc = info[2]+'_100'
#                    break
#                elif int(line.split()[3]) < int(info[0]):
#                    te_loc = 'noTE'
#                    break
        except (ValueError, KeyError):
            continue
#
        if te_loc == 'null':
            try:
                for info in dict_bed[line.split()[1]][(dict_bed[line.split()[1]+'_idx'].get(int((int(line.split()[3]))/10000),1))-1:]:
                    if int(line.split()[3]) in range( int(info[0])-100, int(info[1])+101 ):
                        te_loc = info[2]+'_100'
                        break
                    elif int(line.split()[3]) < int(info[0]):
                        break
            except (ValueError, KeyError):
                continue

        if te_loc == 'null':
            te_loc = 'noTE' 

        words = line.split()

#        samfile = pysam.AlignmentFile(bam_full, 'rb')
#        cnt_rd = 0
#        for read in samfile.fetch(words[1], int(words[3])-insert_size, int(words[3])+insert_size):
#            if read.cigarstring != None:
#                if read.get_reference_positions()[0] >= int(words[3])-insert_size and read.get_reference_positions()[-1] <= int(words[3])+insert_size:
#                    cnt_rd += 1
#        samfile.close()

        test_class = words[0]

        if (words[7] == test_class and float(words[8]) > qual_lim) or words[7] == 'NA':    #left clipped read
            left_clipped_rd = int(words[9])
#        if words[7] == test_class or words[7] == 'NA':
#            test_class_score += 1
        if (words[11] == test_class and float(words[12]) > qual_lim) or words[11] == 'NA':    # right clipped reads
            right_clipped_rd = int(words[13])
        # Left poly-A/T
            p_pat_num = int(words[19])
        # Right poly-A/T
            n_pat_num = int(words[20])
#        if words[11] == test_class or words[11] == 'NA':
#            test_class_score += 1
        if (words[21] == test_class and float(words[22]) > qual_lim) or words[21] == 'NA':    # left discordant reads
            left_discord_rd = int(words[23])
#        if words[19] == test_class or words[19] == 'NA':
#            test_class_score += 1
        if (words[25] == test_class and float(words[26]) > qual_lim) or words[25] == 'NA':    # right discordant reads
            right_discord_rd = int(words[27])
#        if words[23] == test_class or words[23] == 'NA':
#            test_class_score += 1

        total_clipped_rd = left_clipped_rd + right_clipped_rd
        total_clipped_rd_wpat = total_clipped_rd + p_pat_num + n_pat_num
        #
        total_discord_rd = left_discord_rd + right_discord_rd

        total_rd_left = left_clipped_rd + left_discord_rd
        total_rd_right = right_clipped_rd + right_discord_rd

        filter_result = 'FAIL'

#        if ((total_clipped_rd + total_discord_rd >= 3) or ( (total_clipped_rd + total_discord_rd >= 2) and (total_clipped_rd_wpat + total_discord_rd >= 5)) ): 
#        if total_clipped_rd >= 1 and ( (total_clipped_rd + total_discord_rd >= 3) or (total_clipped_rd_wpat + total_discord_rd >= 5) ) : # and test_class_score == 4 \ CEU
#        if total_clipped_rd >= 2 or ( total_clipped_rd >=1 and (total_clipped_rd + total_discord_rd >= 5)) or total_discord_rd >= 10:
#        if total_clipped_rd >=1 and (total_clipped_rd_wpat + total_discord_rd >= 3):
#        if total_clipped_rd >= tcr and (total_clipped_rd + total_discord_rd) >= trd: # and ((total_clipped_rd + total_discord_rd)*100/cnt_rd > 0):
        if total_clipped_rd >= 3:
            filter_result = 'PASS_C3'
        elif total_clipped_rd >= 1 and ( (total_clipped_rd + total_discord_rd >= 3) or (total_clipped_rd_wpat + total_discord_rd >= 5) ) :
            filter_result = 'PASS_C1'
        elif ( total_discord_rd >= 3 ):
            filter_result = 'PASS_D'

        #log_FH.write(filter_result+'\t'+te_loc+'\t'+line) #+'\t'+str(cnt_rd))
        sys.stdout.write(filter_result+'\t'+te_loc+'\t'+line) #+'\t'+str(cnt_rd))
        

        left_clipped_rd, right_clipped_rd, left_discord_rd, right_discord_rd, test_class_score = 0, 0, 0, 0, 0
    log_FH.close()
