
import sys, os
import numpy as np

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def check_file(file_name):
    if os.path.isfile(file_name) and os.path.getsize(file_name) > 0:
        return(True)
    else:
        return(False) 


def num_of_lines(fname):
    i = -1 # in case file in empty
    with open(fname) as fl:
        for i, l in enumerate(fl):
            pass
    fl.close()
    return(i+1)


def extract_line(line_number, filename):
    line_value = []
    for i,file_line in enumerate(open(filename,'r')):
        if i+1==line_number: #i is 0 based
            line_value.append([file_line.split()[0], file_line.split()[1]])
            break
    return(line_value)


def print_tup(inp_tup, noftf, seperation, end_chr):
    for typ_dat in inp_tup:
        noftf.write(seperation.join(str(dat) for dat in typ_dat) + end_chr)
 
        
def get_class(file_name, args):
    #
    qual_interval_inp = args.qii_inp
    #
    te_class_utop = None
    qual_bound = 1
    #
    with open(file_name,'r') as input_file:
        input_file_lines = input_file.readlines()
    input_file.close()

    while (te_class_utop == None):
        te_class = []
        qual_bound = qual_bound - qual_interval_inp
        for line in input_file_lines:
            words = line.strip().split()
            if float(words[7]) >= qual_bound:
                te_class.append(words[3])
        if len(te_class) > 0:
            te_class_u, te_counts = np.unique(te_class, return_counts=True)
            te_class_utop = (te_class_u[np.argmax(te_counts)])
    #
    out_arr = []
    ref_cnt_max = 0
    ref_cnt_min = sys.maxsize

    for text in input_file_lines:
        words = text.strip().split()
        if te_class_utop == words[3]:
            if int(words[5]) > int(ref_cnt_max):
                ref_cnt_max = words[5]
            if int(words[4]) < int(ref_cnt_min):
                ref_cnt_min = words[4]
    out_arr.append([te_class_utop, ref_cnt_min, ref_cnt_max, np.amax(te_counts)])
    #
    return(tuple(out_arr))

def read_type_info(rd_type, end, file_name, args):
    # Decides TE type and some flags.
    #
    output_file_lines = []
    flag_value = 'n'
    #
    output_file_lines.append( 'Mapping of ('+ end +')end '+ rd_type +' reads.. \n' )    
    #
    if check_file(file_name): # file_name = .fa.map
        #map_flname = (file_name)
        with open(file_name,'r') as fi: 
            output_file_lines.append(fi.read())
        fi.close()
        output_file_lines.append('TE type: ')
        flag_value = 'y'

    type_info = get_type(file_name, args)
    output_file_lines.append(';'.join(str(dat) for dat in type_info)+' ')        
    output_file_lines.append('\n\n')
    #

    return(flag_value, type_info, output_file_lines)


def read_class_info(rd_type, end, file_name, args):

    # Finds TE class
    #
    output_file_lines = []
    #
    if check_file(file_name):
        output_file_lines.append( 'Mapping of ('+ end +')end '+ rd_type +' reads.. \n' )
        #
        #map_flname = (file_name)
        with open(file_name,'r') as fi:
            output_file_lines.append(fi.read())
        fi.close()
        output_file_lines.append('TE class: ')
        #
        class_info = get_class(file_name, args)
        output_file_lines.append(';'.join(str(dat) for dat in class_info)+' ')
        output_file_lines.append('\n\n')
        #
    if not check_file(file_name):
        output_file_lines.append( file_name.split('_')[0] +' '+ file_name.split('_')[1] +' '+ \
                    rd_type +'reads/mates of ('+ end +')end do not align to any type\n')

    return(class_info, output_file_lines)

def te_type_setup(file_name_p, file_name_n, type_p, type_n, cnt_p, cnt_n, fofn_ref_file):
    #
    #fofn_ref_file = args.fofn_ref    
    te_class_file = 'none'
    flag_value = 'n'
    output_file_lines = []
    #
    if check_file(file_name_p) == True and check_file(file_name_n) == True:
        ratio_p ='na'
        ratio_n = 'na'
        if cnt_p > 0:
            ratio_p = str(float("{0:.2f}".format((type_p[0][2]/cnt_p)*100)))
        if cnt_n > 0:
            ratio_n = str(float("{0:.2f}".format((type_n[0][2]/cnt_n)*100))) 
        if type_p[0][0] == 'NA' or type_n[0][0] == 'NA' or (type_p[0][0] != type_n[0][0]):
            output_file_lines.append("%s (%s percent) reads at 5'-end support %s with quality >=%s, \
                and %s (%s percent) reads at 3'-end support %s with quality >=%s.\
                Will not go for class and length determination.\n" \
                % (str(type_p[0][2]), ratio_p, \
                str(type_p[0][0]), str(type_n[0][1]), str(type_n[0][2]), ratio_n, \
                str(type_n[0][0]), str(type_n[0][1])))
            flag_value = 'n'
        else:
            output_file_lines.append("%s (%s percent) reads at 5'-end support %s with quality >=%s, \
                and %s (%s percent) reads at 3'-end support %s with quality >=%s.\
                Going for class and length(for clipped only) determination.\n" \
                % (str(type_p[0][2]), ratio_p, \
                str(type_p[0][0]), str(type_n[0][1]), str(type_n[0][2]), ratio_n, \
                str(type_n[0][0]), str(type_n[0][1])))
            flag_value = 'y'
            #
            with open(fofn_ref_file, 'r') as ref_type_file_file:
                ref_type_file_lines = ref_type_file_file.readlines()
            ref_type_file_file.close()
            #
            for line in ref_type_file_lines:
                if type_p[0][0] == line.split()[0]:
                    te_class_file = line.split()[1]
                    break
    
    return(flag_value, te_class_file, output_file_lines)


def get_type(file_name, args):

    # Scans through *.fa.map in various quality intervals, 
    # finds type with highest occurance in each interval.  
    start_qual = 1.0
    te_type_name = []
    qual_interval = args.qii_inp #qual_interval_inp
    num_interval = args.nii_inp #num_interval_inp
    te_type_u = []
    te_counts = []
    #
    if check_file(file_name):
        #
        with open(file_name, 'r') as file_inp:
            file_inp_lines = file_inp.readlines()
        file_inp.close()
        #
        while ( num_interval != 0 ):
            te_type = []
            for line in file_inp_lines:
                #
                words = line.strip().split()
                type_id = words[3].split(':')[1]
                #
                if ( float(words[7]) > start_qual-qual_interval ) and ( float(words[7]) <= start_qual ):
                    te_type.append(type_id)
            #
            if len(te_type) ==0:
                start_qual = start_qual - qual_interval
                num_interval -= 1
                continue
            #
            te_type_u, te_counts = np.unique(te_type, return_counts=True)
            te_type_name.append([te_type_u[np.argmax(te_counts)], \
                    float("{0:.2f}".format(start_qual-qual_interval)), np.amax(te_counts)])
            #
            start_qual = start_qual - qual_interval
            num_interval -= 1
            try:
                del te_type_u[:]
                del te_counts[:]
            except ValueError:
                pass

    else:
        te_type_name.append(['NA', 0, 0])

    if len(te_type_name) == 0: #for very low quality 
        te_type_name.append(['NA', 0, 0])

    return(te_type_name)


def output_header ():
    header = ''
    header = "\t".join( "#Type",
                        "Chromosome",
                        "Initial_Guess",
                        "Estimated_insertion_point",
                        "Hetrozygous/Homozygous",
                        "num_reads_forHet",
                        "TSD_length",
                        "(+)clipped_type",
                        "(+)clipped_type_quality",
                        "num_reads_supporting_(+)clipped_type",
                        "frac_reads_supporting_(+)clipped_type",
                        "(-)clipped_type",
                        "(-)clipped_type_quality",
                        "num_reads_supporting_(-)clipped_type",
                        "frac_reads_supporting_(-)clipped_type",
                        "Estimated_TE_class",
                        "Both_clipped_end_support",
                        "Estimated_TE_Length(bp)",
                        "gap_bw_ends",
                        "num_pat_p",
                        "num_pat_n",
                        "(+)discord_mate_type",
                        "(+)discord_mate_type_quality",
                        "num_reads_supporting_(+)discord_mate_type",
                        "frac_reads_supporting_(+)discord_mate_type",
                        "(-)discord_mate_type",
                        "(-)discord_mate_type_quality",
                        "num_reads_supporting_(-)discord_mate_type",
                        "frac_reads_supporting_(-)discord_mate_type",
                        "Estimated_TE_class-discord",
                        "Both_end_discord_mate_support",
                        "special_comment" )
    return header