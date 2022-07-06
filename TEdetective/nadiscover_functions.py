from TEdetective.general_functions import cgr_to_mpb, cigar_to_tup
#from TEdetective.io_functions import read_type_info, check_file, read_class_info

def alt_mapped_pos( read, args):
    #
    #
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