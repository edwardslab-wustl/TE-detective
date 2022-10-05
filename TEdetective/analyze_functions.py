import os
import subprocess

from Bio.Seq import Seq

from TEdetective.general_functions import pat_check, check_uniq_mapping

def run_censor (fasta_file, library_file, log_FH, message):
    log_FH.write(message + "\n")
    log_FH.flush()
    command = ['censor.ncbi', fasta_file, '-lib', library_file]
    result = subprocess.run(command, stderr = log_FH, stdout = log_FH)
    return

def fix_te_sequence_file_path (te_class_file, dir_path, fofn_ref_realpath, log_FH):
    tmp_file1 = os.path.realpath(te_class_file)
    tmp_file2 = dir_path + "/" + te_class_file
    tmp_file3 = os.path.dirname(fofn_ref_realpath) + "/" + te_class_file
    if os.path.exists(tmp_file1):
        te_class_file = tmp_file1
    elif os.path.exists(tmp_file2):
        te_class_file = tmp_file2
    elif os.path.exists(tmp_file3):
        te_class_file = tmp_file3
    elif te_class_file != 'none':
        log_FH.write("can't find file " + te_class_file+ " from input file "+ fofn_ref_realpath + "\n")
        log_FH.write(tmp_file1 + "\n")
        log_FH.write(tmp_file2 + "\n")
        log_FH.write(tmp_file3 + "\n")
    log_FH.write("Using rep file: " + te_class_file + "\n")
    return te_class_file

def analysis_pat_check( inp_fa_file, inp_fa_map_file ):
    with open(inp_fa_file, 'r') as fa_file:
        fa_file_lines = fa_file.readlines()
    fa_file.close()
    mapeed_seq_idx = []
    try:
        with open(inp_fa_map_file, 'r') as map_file:
            map_file_lines = map_file.readlines()
        map_file.close()
        for line in  map_file_lines:
            mapeed_seq_idx.append( line.strip().split()[0] )
    except FileNotFoundError:
        pass
    fa_file_lines_itr = iter( fa_file_lines )
    pat_test_out_lines = []
    for line in fa_file_lines_itr:
        info_line = line
        seq_line = next( fa_file_lines_itr )
        seq_idx = (info_line.strip().split()[0]).strip('>')
        if seq_idx in mapeed_seq_idx:
            continue
        pat_flag, pat_type = pat_check( seq_line, 9, 1 )
        if pat_flag == 1:
            pat_test_out_lines.append([seq_idx, pat_type ])
    return ( pat_test_out_lines )


def flt_discord( discord_mate_file ):
    with open(discord_mate_file, 'r') as inp_file:
        inp_file_lines = inp_file.readlines()
    inp_file.close()
    inp_file_lines_itr = iter( inp_file_lines )
    dict_info = {}
    for line in inp_file_lines_itr:
        info_line = line
        seq_id = line.strip().split()[1]
        seq_line = next( inp_file_lines_itr )
        seq = Seq( seq_line.strip() )
        try:
            for item in dict_info[seq_id]:
                if ( ( item == seq ) or ( item == seq.complement() ) \
                        or ( item == seq.reverse_complement() ) ):
                    continue
                dict_info[seq_id].append(seq)
        except KeyError:
            dict_info[seq_id] = []
            dict_info[seq_id].append(seq)
    #backup file
    os.rename(discord_mate_file, discord_mate_file+'.prefltr')
    with open(discord_mate_file, 'w') as out_file:    
        idx_cnt = 0
        for key in dict_info:
            for item in dict_info[key]:
                idx_cnt += 1
                out_file.write('>'+str(idx_cnt)+' '+key+'\n'+str(item)+'\n')
    out_file.close()


def find_clipped_ends( iterator_reads, insert_guess, args ):
    array_p_ps = []
    array_n_ps = []
    insert_size = args.isz_inp
    clipped_length = args.cll_inp
    anchor_length = args.ahl_inp
    insert_range = args.rdl_inp
    min_mapq = args.mpq_inp
    
    # 1. Read is uniqly mapped.
    # 2. Soft clipped at either ends. Clipped len > min clipped len, Anchor len > min anchor len
    # 3. For reads p end of TE: insert_guess-read_len <= ref_end < insert_guess+read_len
    # 4. For reads n end of TE: insert_guess-read_len <= ref_start < insert_guess+read_len
             
    for read in iterator_reads:
        if read.reference_start in range(insert_guess-insert_size, insert_guess+insert_size) and \
            read.reference_end in range(insert_guess-insert_size, insert_guess+insert_size) and \
            (read.cigartuples != None) and (read.is_supplementary != True) and (read.has_tag('XA') or \
            read.has_tag('SA') or read.mapping_quality >= min_mapq) and (read.mapping_quality >= args.mpqu_inp):

            write_flag = check_uniq_mapping( read, args ) 

            if read.cigartuples[-1][0] == 4 and read.cigartuples[-1][1] > clipped_length and write_flag == 'y': # soft clipped at p-end of TE
                if read.reference_end-read.reference_start > anchor_length and read.reference_end \
                    in range(insert_guess-insert_range,insert_guess+insert_range): # |------>|-read_length|*|___
                    array_p_ps.append(read.reference_end)

            if read.cigartuples[0][0] == 4 and read.cigartuples[0][1] > clipped_length and write_flag == 'y': # soft clipped at n-end of TE
                if read.reference_end-read.reference_start > anchor_length and read.reference_start \
                    in range(insert_guess-insert_range,insert_guess+insert_range): #__|*|+read_length|<------|
                    array_n_ps.append(read.reference_start)

    array_p_s = sorted(array_p_ps, key=int)
    array_n_s = sorted(array_n_ps, key=int)
    return(array_p_s, array_n_s)