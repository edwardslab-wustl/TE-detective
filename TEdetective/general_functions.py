import re

import numpy as np

from TEdetective.io_functions import check_file


def cigar_to_tup(cigar_string, map_dir):
    cigar_flag = {'M':0, 'I':1, 'D':2, 'N':3, 'S':4, 'H':5, 'P':6, '=':7, 'X':8, 'B':9}
    cgr_list = re.findall(r'[A-Za-z]+|\d+', cigar_string)
    cgr_list_grp = [cgr_list[i:i + 2] for i in range(0, len(cgr_list), 2)]
    if map_dir == 'n':
        cgr_list_grp = cgr_list_grp[::-1]
    cigar_tup = []
    for data in cgr_list_grp:
        cigar_tup.append((cigar_flag[data[1]], int(data[0])))
    return(cigar_tup)


def cgr_to_mpb(cigar_tup):
    i= 0
    mpb_lst = []
    for info in cigar_tup:
        for j in range(0, info[1]):
            i += 1
            if info[0] == 0:
                mpb_lst.append(i)
    return(mpb_lst)


def check_uniq_mapping( read, args ):

    write_flag = 'y'
    
    if read.mapping_quality < args.mpq_inp:
        #
        read_cgr_tup = read.cigartuples
        if read.is_reverse:
            read_cgr_tup = read_cgr_tup[::-1]
        maped_bases = cgr_to_mpb(read_cgr_tup)
        if read.has_tag('XA'): # Secondary alignment
            tag_line = read.get_tag('XA')
            for i in range(0, len(tag_line.split(';'))-1): # why -1 ? -> ends with ; .
                align_info = tag_line.split(';')[i]
                if int(align_info.split(',')[1]) > 0:
                    secondary_maped_bases = cgr_to_mpb(cigar_to_tup(align_info.split(',')[2], 'p'))
                if int(align_info.split(',')[1]) < 0:
                    secondary_maped_bases = cgr_to_mpb(cigar_to_tup(align_info.split(',')[2], 'n'))
                if len(set(maped_bases).intersection(secondary_maped_bases)) > 5:
                    write_flag = 'n'

        if read.has_tag('SA'): # Chimeric alignment
            sa_tag_line = read.get_tag('SA')
            for i in range(0, len(sa_tag_line.split(';'))-1):
                sa_align_info = sa_tag_line.split(';')[i]
                if sa_align_info.split(',')[2] == '+':
                    sa_secondary_maped_bases = cgr_to_mpb(cigar_to_tup(sa_align_info.split(',')[3], 'p'))
                if sa_align_info.split(',')[2] == '-':
                    sa_secondary_maped_bases = cgr_to_mpb(cigar_to_tup(sa_align_info.split(',')[3], 'n'))
                if len(set(maped_bases).intersection(sa_secondary_maped_bases)) > 5:
                    write_flag = 'n'

    return( write_flag )


def break_points_2d(data_set, clust_denst, end_range, offset_value):
    clusters_data=[]
    ref_chrom = 'None'
    tnum_dat = []
    for items in data_set:
        if items[0] != ref_chrom:
            if len(tnum_dat[:-1]) >= clust_denst:
                clusters_data.append([ ref_chrom,int(round(clust_ref_point+offset_value)),len(tnum_dat[:-1]), tnum_dat[:-1] ])
            ref_chrom = items[0]
            del tnum_dat[:-1]
        tnum_dat.append(items[1])
        ref_point = np.mean(np.array(tnum_dat))

        if (items[1] >= (ref_point - end_range)) and (items[1] <= (ref_point + end_range)):
            clust_ref_point = ref_point
        elif len(tnum_dat) >= 3 and ( (items[1] >= (np.mean(np.array(tnum_dat[1:])) - end_range)) \
            and  (items[1] <= (np.mean(np.array(tnum_dat[1:])) + end_range)) ) \
            and np.std(np.array(tnum_dat[1:])) < np.std(np.array(tnum_dat[:-1])):
            del tnum_dat[0]
            clust_ref_point = np.mean(np.array(tnum_dat))
        else:
            if len(tnum_dat[:-1]) >= clust_denst:
                clusters_data.append([ ref_chrom,int(round(clust_ref_point+offset_value)),len(tnum_dat[:-1]), tnum_dat[:-1] ])
            del tnum_dat[:-1]
    if len(tnum_dat) >= clust_denst:
#        clusters_data.append((ref_chrom,int(round(clust_ref_point)),len(tnum_dat)))
        clusters_data.append([ ref_chrom, int(round(clust_ref_point+offset_value)), len(tnum_dat), tnum_dat ])

    return(clusters_data)


def break_points(data_set, clust_denst, end_range):
    clusters_data=[]
    tnum_dat = []
    for tnum in data_set:
        tnum_dat.append(tnum)
        ref_point = np.mean(np.array(tnum_dat))
        if (tnum >= (ref_point - end_range)) and (tnum <= (ref_point + end_range)):
            clust_ref_point = ref_point # we need this to print
        elif len(tnum_dat) >= 3 and ( (tnum >= (np.mean(np.array(tnum_dat[1:])) - end_range)) \
            and  (tnum <= (np.mean(np.array(tnum_dat[1:])) + end_range)) ) \
            and np.std(np.array(tnum_dat[1:])) < np.std(np.array(tnum_dat[:-1])):
            del tnum_dat[0]
            clust_ref_point = np.mean(np.array(tnum_dat))
        else:
            if len(tnum_dat[:-1]) >= clust_denst:
                #if out_flag == 0:
                clusters_data.append( int(round(clust_ref_point)) )
                #if out_flag == 1:
                #    clusters_data.append( ( int(round(clust_ref_point)), len(tnum_dat)  ) )
            del tnum_dat[:-1]

    if len(tnum_dat) >= clust_denst:

        #if out_flag == 0:
        clusters_data.append(int(round(clust_ref_point)))
        #if out_flag == 1:
        #    clusters_data.append( ( int(round(clust_ref_point), len(tnum_dat)  ) )

    return(clusters_data)


def calc_length(inp_tup_a, inp_tup_b):
    map_index = []
    map_index.extend((int(inp_tup_a[0][1]), int(inp_tup_a[0][2]), int(inp_tup_b[0][1]), int(inp_tup_b[0][2])))
    length = np.amax(map_index) - np.amin(map_index)
    return(length)


def pat_check(inp_seq, query_len, mis_match):
    inp_seq = inp_seq.upper()
    match_len = query_len - mis_match
    polyAT_test_flag = 0
    polyAT_type = 'X'
    # 1. Check poly-A
    start_idx  = 0
    end_idx = query_len
    while ( end_idx <= len(inp_seq) ):
        if ( inp_seq[start_idx:end_idx].count('A') >= match_len ):
            polyAT_test_flag  = 1
            polyAT_type = 'A'
            break
        elif ( inp_seq[start_idx:end_idx].count('T') >= match_len ):
            polyAT_test_flag  = 1
            polyAT_type = 'T'
            break
        start_idx += 1
        end_idx += 1
    return( polyAT_test_flag, polyAT_type )
