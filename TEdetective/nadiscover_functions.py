import pysam

from TEdetective.general_functions import cgr_to_mpb, cigar_to_tup
from TEdetective.general_functions import check_uniq_mapping,pat_check
from TEdetective.io_functions import eprint
#from TEdetective.io_functions import read_type_info, check_file, read_class_info

def alt_mapped_pos( read, args):
    if read.has_tag('XA'): # Secondary alignment
        tag_line = read.get_tag('XA')
#        for i in range(0, len(tag_line.split(';'))-1): # why -1 ? -> ends with ; .
        align_info = tag_line.split(';')[0]
        if int(align_info.split(',')[1]) > 0:
            secondary_maped_bases = cgr_to_mpb(cigar_to_tup(align_info.split(',')[2], 'p'))
            chrom = align_info.split(',')[0]
            reg_start = abs(int(align_info.split(',')[1]))
            reg_end = reg_start + secondary_maped_bases[-1]
        if int(align_info.split(',')[1]) < 0:
            secondary_maped_bases = cgr_to_mpb(cigar_to_tup(align_info.split(',')[2], 'n'))
            chrom = align_info.split(',')[0]
            reg_end = abs(int(align_info.split(',')[1]))
            reg_start = reg_end - secondary_maped_bases[-1]

    elif read.has_tag('SA'): # Chimeric alignment Ex:
        sa_tag_line = read.get_tag('SA')
#        for i in range(0, len(sa_tag_line.split(';'))-1):
        sa_align_info = sa_tag_line.split(';')[0]
        if sa_align_info.split(',')[2] == '+':
            sa_secondary_maped_bases = cgr_to_mpb(cigar_to_tup(sa_align_info.split(',')[3], 'p'))
            chrom = sa_align_info.split(',')[0]
            reg_start = int(sa_align_info.split(',')[1])
            reg_end = reg_start + sa_secondary_maped_bases[-1]
        if sa_align_info.split(',')[2] == '-':
            sa_secondary_maped_bases = cgr_to_mpb(cigar_to_tup(sa_align_info.split(',')[3], 'n'))
            chrom = sa_align_info.split(',')[0]
            reg_end = int(sa_align_info.split(',')[1])
            reg_start = reg_end - sa_secondary_maped_bases[-1]

    return( chrom, reg_start, reg_end )


def write_pat_clipped_dat(read_bam_clipped, min_mapq_uniq, clipped_length, pat_out_file, pat_query_len, pat_mis_match, args ):
        pat_out_file_lines = []
        samfile_clipped = pysam.AlignmentFile(read_bam_clipped, "rb")
        for read in samfile_clipped.fetch():
            #eprint(read.is_supplementary, read.has_tag('XA'), read.has_tag('SA'), read.mapping_quality)
            if (read.is_supplementary != True) and ( read.has_tag('XA') or read.has_tag('SA') ) \
                and ( read.mapping_quality >= min_mapq_uniq ):
                
                write_flag = check_uniq_mapping( read, args )
                #eprint(write_flag, read.is_supplementary, read.has_tag('XA'), read.has_tag('SA'), read.mapping_quality)
                eprint(read.cigartuples)
                if write_flag == 'y':
                    clipped_side = 'X'
                    query_sequence = 'tmp'
                    if  ( read.cigartuples[-1][0] == 4  
                         and read.cigartuples[-1][1] > clipped_length 
                         and (read.cigartuples[0][0] == 4 and read.cigartuples[0][1] > 5) != True ):
                        clipped_side = 'R'
                        query_sequence=read.query_sequence[read.infer_query_length() \
                                    - read.cigartuples[-1][1]-1:read.infer_query_length()]
                        eprint('top')
                    elif ( read.cigartuples[0][0] == 4 ) and ( read.cigartuples[0][1] > clipped_length ) and \
                        ( (read.cigartuples[-1][0] == 4 and read.cigartuples[-1][1] > 5) != True ):
                        clipped_side = 'L'
                        query_sequence=read.query_sequence[0:read.cigartuples[0][1]-1]
                        eprint('bot')
                        
                    pat_flag, pat_type = pat_check(query_sequence, pat_query_len, pat_mis_match)
                    eprint (pat_flag, pat_type, pat_query_len, pat_mis_match, query_sequence)
                    if pat_flag == 1:
                        pat_out_file_lines.append( read.query_name +' '+ str(read.flag) + ' ' \
                            + str(read.reference_name) + ' ' + str(read.reference_start) + ' ' + \
                                         str(read.reference_end) + ' ' + clipped_side )
                            # +'\t'+query_sequence+'\t'+pat_type )
                        #pat_dat_lines.append( (read.reference_name, read.reference_start) )
        samfile_clipped.close()
        #pat_out_file = open(preprocess_dir_realpath+'/pat_clipped_read-bam_id_flag.dat','w')
        pat_out_file.write('\n'.join(pat_out_file_lines))
        #del pat_out_file_lines[:]
        pat_out_file.close()
        return pat_out_file_lines
    
def write_pat_clipped_dat_new(full_bam, read_bam_clipped, min_mapq_uniq, clipped_length, pat_out_file, pat_query_len, pat_mis_match, args ):
        pat_out_file_lines = []
        samfile_clipped = pysam.AlignmentFile(read_bam_clipped, "rb")
        bamfile = pysam.AlignmentFile(full_bam, 'rb')
        full_bam_index =  pysam.IndexedReads(bamfile)
        full_bam_index.build()
        for orig_read in samfile_clipped.fetch():
            try:
                full_bam_index.find(orig_read.query_name)
            except KeyError:
                #eprint("1, can't find: ", orig_read.query_name)
                pass
            else:
                iterator = full_bam_index.find(orig_read.query_name)
                if iterator == None:
                    eprint("Can't find: ", orig_read.query_name)
                else:
                    for read in iterator:
                        #eprint(read.is_supplementary, read.has_tag('XA'), read.has_tag('SA'), read.mapping_quality)
                        if ( read.is_supplementary != True  
                             and ( read.has_tag('XA') or read.has_tag('SA') ) 
                             and  read.mapping_quality >= min_mapq_uniq ):
                            write_flag = check_uniq_mapping( read, args )
                            #eprint(write_flag, read.is_supplementary, read.has_tag('XA'), read.has_tag('SA'), read.mapping_quality)
                            #eprint(read.query_name, read.cigartuples)
                            if write_flag == 'y':
                                clipped_side = 'X'
                                query_sequence = 'tmp'
                                if ( read.cigartuples[-1][0] == 4  
                                    and read.cigartuples[-1][1] > clipped_length 
                                    and (read.cigartuples[0][0] == 4 and read.cigartuples[0][1] > 5) != True ):
                                    clipped_side = 'R'
                                    start = read.infer_query_length() - read.cigartuples[-1][1]-1
                                    end = read.infer_query_length() 
                                    query_sequence=read.query_sequence[start:end]
                                    #eprint('top')
                                elif ( read.cigartuples[0][0] == 4 
                                      and read.cigartuples[0][1] > clipped_length 
                                      and ( (read.cigartuples[-1][0] == 4 and read.cigartuples[-1][1] > 5) != True ) ):
                                    clipped_side = 'L'
                                    start = 0
                                    end = read.cigartuples[0][1]-1
                                    query_sequence=read.query_sequence[start:end]
                                    #eprint('bot')
                                pat_flag, pat_type = pat_check(query_sequence, pat_query_len, pat_mis_match)
                                eprint (pat_flag, pat_type, pat_query_len, pat_mis_match, query_sequence)
                                if pat_flag == 1:
                                    outLine = ' '.join([ read.query_name,
                                                        str(read.flag),
                                                        str(read.reference_name),
                                                        str(read.reference_start),
                                                        str(read.reference_end), 
                                                        clipped_side ])
                                    pat_out_file_lines.append(outLine)
                                    #pat_out_file_lines.append( read.query_name +' '+ str(read.flag) + ' ' \
                                    #    + str(read.reference_name) + ' ' + str(read.reference_start) + ' ' + \
                                    #                 str(read.reference_end) + ' ' + clipped_side )
                                    #OLD     # +'\t'+query_sequence+'\t'+pat_type )
                                    #OLD #pat_dat_lines.append( (read.reference_name, read.reference_start) )
        samfile_clipped.close()
        bamfile.close()
        pat_out_file.write('\n'.join(pat_out_file_lines))
        pat_out_file.close()
        return pat_out_file_lines
    
def write_pat_clipped_dat_new2(full_bam, read_bam_clipped_cmpl, min_mapq_uniq, clipped_length, pat_out_file, pat_query_len, pat_mis_match, args ):
        pat_out_file_lines = []
        samfile_clipped_cmpl = pysam.AlignmentFile(read_bam_clipped_cmpl, "rb")
        #bamfile = pysam.AlignmentFile(full_bam, 'rb')
        #full_bam_index =  pysam.IndexedReads(bamfile)
        #full_bam_index.build()
        for read in samfile_clipped_cmpl.fetch():
            #eprint(read.is_supplementary, read.has_tag('XA'), read.has_tag('SA'), read.mapping_quality)
            if ( read.is_supplementary != True  
                 and ( read.has_tag('XA') or read.has_tag('SA') ) 
                 and  read.mapping_quality >= min_mapq_uniq ):
                write_flag = check_uniq_mapping( read, args )
                #eprint(write_flag, read.is_supplementary, read.has_tag('XA'), read.has_tag('SA'), read.mapping_quality)
                #eprint(read.query_name, read.cigartuples)
                if write_flag == 'y':
                    clipped_side = 'X'
                    query_sequence = 'tmp'
                    if ( read.cigartuples[-1][0] == 4  
                        and read.cigartuples[-1][1] > clipped_length 
                        and (read.cigartuples[0][0] == 4 and read.cigartuples[0][1] > 5) != True ):
                        clipped_side = 'R'
                        start = read.infer_query_length() - read.cigartuples[-1][1]-1
                        end = read.infer_query_length() 
                        query_sequence=read.query_sequence[start:end]
                        #eprint('top')
                    elif ( read.cigartuples[0][0] == 4 
                          and read.cigartuples[0][1] > clipped_length 
                          and ( (read.cigartuples[-1][0] == 4 and read.cigartuples[-1][1] > 5) != True ) ):
                        clipped_side = 'L'
                        start = 0
                        end = read.cigartuples[0][1]-1
                        query_sequence=read.query_sequence[start:end]
                        #eprint('bot')
                    pat_flag, pat_type = pat_check(query_sequence, pat_query_len, pat_mis_match)
                    #eprint (pat_flag, pat_type, pat_query_len, pat_mis_match, query_sequence)
                    if pat_flag == 1:
                        outLine = ' '.join([ read.query_name,
                                            str(read.flag),
                                            str(read.reference_name),
                                            str(read.reference_start),
                                            str(read.reference_end), 
                                            clipped_side ])
                        pat_out_file_lines.append(outLine)
                        #pat_out_file_lines.append( read.query_name +' '+ str(read.flag) + ' ' \
                        #    + str(read.reference_name) + ' ' + str(read.reference_start) + ' ' + \
                        #                 str(read.reference_end) + ' ' + clipped_side )
                        #OLD     # +'\t'+query_sequence+'\t'+pat_type )
                        #OLD #pat_dat_lines.append( (read.reference_name, read.reference_start) )
        samfile_clipped_cmpl.close()
        pat_out_file.write('\n'.join(pat_out_file_lines))
        pat_out_file.close()
        return pat_out_file_lines