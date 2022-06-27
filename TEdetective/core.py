################################################################
# Python based tool to identify novel transposable element insertions from WGS data.
# Contributors: Manoj Kumar Singh (manoj@wustl.edu),
#        and John Edwards (jredwards@wustl.edu)
#
# External dependencies: censor, NCBI blast (provided with package)
#
################################################################

import sys
import os
import copy
#from tkinter import W
import pysam
import numpy as np
import argparse
import re
import subprocess
from subprocess import Popen, PIPE
from Bio.Seq import Seq
from Bio.Sequencing.Applications import BwaIndexCommandline
from Bio.Sequencing.Applications import BwaAlignCommandline
from Bio.Sequencing.Applications import BwaSamseCommandline
from TEdetective.io_functions import eprint


#--------------------------------------------------------------
def check_file(file_name):
    if os.path.isfile(file_name) and os.path.getsize(file_name) > 0:
        return(True)
    else:
        return(False)    
#--------------------------------------------------------------
def num_of_lines(fname):
    i = -1 # in case file in empty
    with open(fname) as fl:
        for i, l in enumerate(fl):
            pass
    fl.close()
    return(i+1)
#--------------------------------------------------------------
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
#--------------------------------------------------------------
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
#--------------------------------------------------------------
def extract_line(line_number, filename):
    line_value = []
    for i,file_line in enumerate(open(filename,'r')):
        if i+1==line_number: #i is 0 based
            line_value.append([file_line.split()[0], file_line.split()[1]])
            break
    return(line_value)
#--------------------------------------------------------------
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

#--------------------------------------------------------------

def print_tup(inp_tup, noftf, seperation, end_chr):
    for typ_dat in inp_tup:
        noftf.write(seperation.join(str(dat) for dat in typ_dat) + end_chr)

#--------------------------------------------------------------

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

#--------------------------------------------------------------

def calc_length(inp_tup_a, inp_tup_b):
    map_index = []
    map_index.extend((int(inp_tup_a[0][1]), int(inp_tup_a[0][2]), int(inp_tup_b[0][1]), int(inp_tup_b[0][2])))
    length = np.amax(map_index) - np.amin(map_index)
    return(length)

#--------------------------------------------------------------

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

#--------------------------------------------------------------

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

#--------------------------------------------------------------

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

#--------------------------------------------------------------

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

#--------------------------------------------------------------

def cgr_to_mpb(cigar_tup):
    i= 0
    mpb_lst = []
    for info in cigar_tup:
        for j in range(0, info[1]):
            i += 1
            if info[0] == 0:
                mpb_lst.append(i)
    return(mpb_lst)

#--------------------------------------------------------------

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

#--------------------------------------------------------------

def exec_preprocess(args):
    #
    dir_path = os.getcwd()
    sys.stdout.write('working directory: '+ dir_path +'\n')
    #
    preprocess_dir_realpath = os.path.realpath(args.preprocess_dir)
    sys.stdout.write('preprocessing/intermediate file directory: '+ preprocess_dir_realpath +'\n')
    #
    fofn_ref_realpath = os.path.realpath(args.fofn_ref)
    sys.stdout.write('fofn_ref file: '+ fofn_ref_realpath +'\n')
    #
    bam_full = os.path.realpath(args.bam_inp)
    sys.stdout.write('Input bam file: '+str(bam_full)+'\n')
    #
    bam_short_name = bam_full.split('/')[-1][:-4]
    sys.stdout.write('bam short name: '+ bam_short_name +'\n')
    #
    clipped_length = args.cll_inp
    sys.stdout.write('Minimum clipped length: '+str(clipped_length)+'\n')
    #
    #fofn_ref_file = args.fofn_ref
    #sys.stdout.write('FoFn reference sequence: '+str(fofn_ref_file)+'\n')
    #
    sys.stdout.write('working directory: '+os.getcwd()+'\n')

    #discord_bam = bam_full.split('/')[-1][:-4]+'_discord.bam'
    #clipped_bam = bam_full.split('/')[-1][:-4]+'_clipped.bam'
    discord_bam = bam_short_name+'_discord.bam'
    clipped_bam = bam_short_name+'_clipped.bam'
    #subprocess.run(['mkdir' , '-p' , args.preprocess_outdir ])
    subprocess.run(['mkdir' , '-p' , preprocess_dir_realpath ])

    samfile = pysam.AlignmentFile(bam_full, "rb")
    #newsam_d = pysam.AlignmentFile('./preprocessed_files/'+discord_bam, 'wb', template=samfile)
    #newsam_c = pysam.AlignmentFile('./preprocessed_files/'+clipped_bam, 'wb', template=samfile)
    #newsam_d = pysam.AlignmentFile(args.preprocess_outdir + '/' + discord_bam, 'wb', template=samfile)
    #newsam_c = pysam.AlignmentFile(args.preprocess_outdir + '/' + clipped_bam, 'wb', template=samfile)
    newsam_d = pysam.AlignmentFile(preprocess_dir_realpath + '/' + discord_bam, 'wb', template=samfile)
    newsam_c = pysam.AlignmentFile(preprocess_dir_realpath + '/' + clipped_bam, 'wb', template=samfile)

    for read in samfile.fetch():
        if read.is_paired == True and read.is_proper_pair != True:
            newsam_d.write(read)
        if read.cigartuples != None:
            if ( read.cigartuples[-1][0] == 4 ) and ( read.cigartuples[-1][1] > clipped_length ) and \
                    ( (read.cigartuples[0][0] == 4 and read.cigartuples[0][1] > 5) != True ):
                a = pysam.AlignedSegment()
                a.query_name = read.query_name
                a.query_sequence=read.query_sequence[read.infer_query_length() \
                            - read.cigartuples[-1][1]-1:read.infer_query_length()]
                a.flag = read.flag
                a.reference_id = read.reference_id
                a.reference_start = read.reference_start
                a.mapping_quality = 0 # read.mapping_quality
                a.cigar = []
                a.cigar.append((read.cigartuples[-1][0], read.cigartuples[-1][1]))
                a.next_reference_id = read.next_reference_id
                a.next_reference_start=read.next_reference_start
                a.template_length=read.template_length
                a.query_qualities = read.query_qualities[read.infer_query_length() \
                            - read.cigartuples[-1][1]-1:read.infer_query_length()]
                a.tags = read.tags
                newsam_c.write(a)

            elif ( read.cigartuples[0][0] == 4 ) and ( read.cigartuples[0][1] > clipped_length ) and \
                    ( (read.cigartuples[-1][0] == 4 and read.cigartuples[-1][1] > 5) != True ):
                a = pysam.AlignedSegment()
                a.query_name = read.query_name
                a.query_sequence=read.query_sequence[0:read.cigartuples[0][1]-1]
                a.flag = read.flag
                a.reference_id = read.reference_id
                a.reference_start = read.reference_start
                a.mapping_quality = 0 #read.mapping_quality
                a.cigar = []
                a.cigar.append((read.cigartuples[0][0], read.cigartuples[0][1])) 
                a.next_reference_id = read.next_reference_id
                a.next_reference_start=read.next_reference_start
                a.template_length=read.template_length
                a.query_qualities=read.query_qualities[0:read.cigartuples[0][1]-1]
                a.tags = read.tags
                newsam_c.write(a)
    newsam_d.close()
    newsam_c.close()
    samfile.close()

    #pysam.index('./preprocessed_files/'+discord_bam)
    #pysam.index('./preprocessed_files/'+clipped_bam)
    #sys.stdout.write('Created ./preprocessed_files/%s\n' % discord_bam)
    #sys.stdout.write('Created ./preprocessed_files/%s\n' % clipped_bam)
    #pysam.index(args.preprocess_outdir + '/' + discord_bam)
    #pysam.index(args.preprocess_outdir + '/' + clipped_bam)
    #sys.stdout.write('Created %s/%s\n' % (args.preprocess_outdir,discord_bam))
    #sys.stdout.write('Created %s/%s\n' % (args.preprocess_outdir,clipped_bam))
    pysam.index(preprocess_dir_realpath + '/' + discord_bam)
    pysam.index(preprocess_dir_realpath + '/' + clipped_bam)
    sys.stdout.write('Created %s/%s\n' % (preprocess_dir_realpath,discord_bam))
    sys.stdout.write('Created %s/%s\n' % (preprocess_dir_realpath,clipped_bam))

    ref_type_file_name = []
#    with open(fofn_ref_file, 'r') as ref_type_file_file:
    with open(fofn_ref_realpath, 'r') as ref_type_file_file:
        for line in ref_type_file_file:
            ref_type_file_name.append([line.split()[0], line.split()[1]])

    #type_file_all = open('./preprocessed_files/te_ref_type.fa', 'w')
    #type_file_all_bwa = open('./preprocessed_files/te_ref_type_bwa.fa', 'w')
#    type_file_all = open(args.preprocess_outdir + '/te_ref_type.fa', 'w')
#    type_file_all_bwa = open(args.preprocess_outdir + '/te_ref_type_bwa.fa', 'w')
    type_file_all = open(preprocess_dir_realpath + '/te_ref_type.fa', 'w')
    type_file_all_bwa = open(preprocess_dir_realpath + '/te_ref_type_bwa.fa', 'w')
    for cnt in range(0, len(ref_type_file_name)):
        i=0
        with open(ref_type_file_name[cnt][1], 'r') as ref_file:
            for line in ref_file:
                line_bwa = line
                if line.startswith('>'):
                    i+=1
                    line = line.split()[0]+':'+ref_type_file_name[cnt][0]+'\n'
                    line_bwa = '>'+ref_type_file_name[cnt][0]+'_'+str(i)+'\n'
                type_file_all.write(line)
                type_file_all_bwa.write(line_bwa)
        ref_file.close()
        cnt += 1
    type_file_all.close()
    type_file_all_bwa.close()

    #sys.stdout.write('Created ./preprocessed_files/te_ref_type.fa\n')
    #sys.stdout.write('Created ./preprocessed_files/te_ref_type_bwa.fa\n')
    #sys.stdout.write('Created ' + args.preprocess_outdir + '/te_ref_type.fa\n')
    #sys.stdout.write('Created ' + args.preprocess_outdir + '/te_ref_type_bwa.fa\n')
    sys.stdout.write('Created ' + preprocess_dir_realpath + '/te_ref_type.fa\n')
    sys.stdout.write('Created ' + preprocess_dir_realpath + '/te_ref_type_bwa.fa\n')

#--------------------------------------------------------------

def exec_discover(args):
    #
    dir_path = os.getcwd() #os.path.dirname(os.path.realpath(__file__))
    sys.stdout.write('Working directory: '+str(dir_path)+'\n')
    #
    preprocess_dir_realpath = os.path.realpath(args.preprocess_dir)
    sys.stdout.write('preprocessing/intermediate file directory: '+ preprocess_dir_realpath +'\n')
    #
    fofn_ref_realpath = os.path.realpath(args.fofn_ref)
    sys.stdout.write('fofn_ref file: '+ fofn_ref_realpath +'\n')
    #
    bam_full = os.path.realpath(args.bam_inp)
    sys.stdout.write('bam file name: '+str(bam_full)+'\n')
    #
    bam_short_name = bam_full.split('/')[-1][:-4]
    sys.stdout.write('bam short name: '+ bam_short_name +'\n')
    #
    insert_size = args.isz_inp
    sys.stdout.write('Insert size estimate: '+str(insert_size)+'\n')
    #
    discord_rd_clust_denst = args.drd_inp # <--change this variable name, confusing.
    sys.stdout.write('Number of reads in a cluster to call it insertion: '+str(discord_rd_clust_denst)+'\n')
    #
    read_length = args.rdl_inp
    sys.stdout.write('Average read length: '+str(read_length)+'\n')
    #
    coverage_cutoff    = args.cct_inp
    sys.stdout.write('Coverage cutoff to skip a region: '+str(coverage_cutoff)+'\n')
    #
    min_mapq = args.mpq_inp
    sys.stdout.write('minimum mapping quality: '+str(min_mapq)+'\n')
    #
    min_mapq_uniq = args.mpqu_inp
    sys.stdout.write('minimum mapping quality: '+str(min_mapq_uniq)+'\n')
    #    
    #fofn_ref_file = args.fofn_ref
    #sys.stdout.write('FoFn reference sequence: '+str(fofn_ref_file)+'\n\n')
    #
    clipped_length = args.cll_inp
    sys.stdout.write('Minimum clipped length: '+str(clipped_length)+'\n')
    #
    #subprocess.run(['mkdir' , '-p' , 'intermediate_files'])
    subprocess.run(['mkdir' , '-p' , preprocess_dir_realpath])
    #
    #reference_genome = dir_path+'/preprocessed_files/te_ref_type_bwa.fa'
    #reference_genome = args.preprocess_dir+'/te_ref_type_bwa.fa'
    reference_genome = preprocess_dir_realpath+'/te_ref_type_bwa.fa'
    index_cmd = BwaIndexCommandline(infile=reference_genome, algorithm='bwtsw')
    index_cmd()
    #
    #read_bam = dir_path+'/preprocessed_files/'+bam_full.split('/')[-1][:-4]+'_discord.bam'
    #read_bam = args.preprocess_dir + '/'+bam_full.split('/')[-1][:-4]+'_discord.bam'
    read_bam = preprocess_dir_realpath + '/'+bam_short_name+'_discord.bam'
    #output_sai_file = dir_path+'/intermediate_files/aligned_reads.sai'
    output_sai_file = preprocess_dir_realpath + '/aligned_reads.sai'
    read_group='@RG ID:foo  SM:bar'
    align_cmd = BwaAlignCommandline(reference=reference_genome, b='b', read_file=read_bam)
    align_cmd(stdout=output_sai_file)
    #
    #output_sam_file = dir_path+'/intermediate_files/aligned_reads.sam'
    output_sam_file = preprocess_dir_realpath + '/aligned_reads.sam'
    samse_cmd = BwaSamseCommandline(reference=reference_genome, read_file=read_bam, sai_file=output_sai_file)
    samse_cmd(stdout=output_sam_file)
    #
    #read_bam_clipped = dir_path+'/preprocessed_files/'+bam_full.split('/')[-1][:-4]+'_clipped.bam'
    #read_bam_clipped = args.preprocess_dir+'/'+bam_short_name+'_clipped.bam'
    read_bam_clipped = preprocess_dir_realpath+'/'+bam_short_name+'_clipped.bam'
    #output_sai_file_clipped = dir_path+'/intermediate_files/aligned_reads_clipped.sai'
    output_sai_file_clipped = preprocess_dir_realpath+ '/aligned_reads_clipped.sai'
    read_group='@RG ID:foo  SM:bar'
    align_cmd_clipped = BwaAlignCommandline(reference=reference_genome, b='b', read_file=read_bam_clipped)
    align_cmd_clipped(stdout=output_sai_file_clipped)
    #
    #output_sam_file_clipped = dir_path+'/intermediate_files/aligned_reads_clipped.sam'
    output_sam_file_clipped = preprocess_dir_realpath + '/aligned_reads_clipped.sam'
    samse_cmd_clipped = BwaSamseCommandline(reference=reference_genome, read_file=read_bam_clipped, sai_file=output_sai_file_clipped)
    samse_cmd_clipped(stdout=output_sam_file_clipped)    
    #
    ref_type_file_name = []
    #with open(fofn_ref_file, 'r') as ref_type_file_file:
    with open(fofn_ref_realpath, 'r') as ref_type_file_file:
        for line in ref_type_file_file:
            ref_type_file_name.append([line.split()[0], line.split()[1]])
    #
    # n is 0 based *********
    dict_aof = {}
    for cnt_1 in range( len(ref_type_file_name) ):
        dict_aof[ref_type_file_name[cnt_1][0]] = []
    #
    with open(output_sam_file,'r') as output_sam_file_fl:
        output_sam_file_fl_lines = output_sam_file_fl.readlines()
    output_sam_file_fl.close
    #
    count_header_lines = 0
    for n,line in enumerate(output_sam_file_fl_lines):  #<<--- read it in ram
        if line.startswith('@'):
            count_header_lines += 1
        if line.startswith('@') == False: # and (int(line.split()[1]) & 0x0004) != 0:
            try:
                dict_aof[line.split()[2].split('_')[0]].append(str(n-count_header_lines+1)+' '+str(line.split()[0]))
            except KeyError:
                continue
    del output_sam_file_fl_lines[:]
    #
    for cnt_1 in range( len(ref_type_file_name) ):
        #with open(dir_path+'/intermediate_files/'+ref_type_file_name[cnt_1][0]+'_ref-aligned_line-num_id.dat', 'w') as temp_fl:
        with open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_ref-aligned_line-num_id.dat', 'w') as temp_fl:
            temp_fl.write('\n'.join(dict_aof[ref_type_file_name[cnt_1][0]]))    
        temp_fl.close
    dict_aof.clear()
    #
    ######
    dict_aofc = {}
    for cnt_1 in range( len(ref_type_file_name) ):
        dict_aofc[ref_type_file_name[cnt_1][0]] = []

    with open(output_sam_file_clipped,'r') as output_sam_file_clipped_fl:
        output_sam_file_clipped_fl_lines = output_sam_file_clipped_fl.readlines()
    output_sam_file_clipped_fl.close()
    #
    count_header_lines = 0
    for n,line in enumerate(output_sam_file_clipped_fl_lines):
        if line.startswith('@'):
            count_header_lines += 1
        if line.startswith('@') == False: # and line.split()[2] != '*':
            try:
                dict_aofc[line.split()[2].split('_')[0]].append(str(n-count_header_lines+1)+' '+str(line.split()[0]))
            except KeyError:
                continue
    del output_sam_file_clipped_fl_lines[:]
    #
    for cnt_1 in range( len(ref_type_file_name) ):
        #with open(dir_path+'/intermediate_files/'+ref_type_file_name[cnt_1][0]+'_clipped_ref-aligned_line-num_id.dat', 'w') as temp_fl:
        with open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_clipped_ref-aligned_line-num_id.dat', 'w') as temp_fl:
            temp_fl.write('\n'.join(dict_aofc[ref_type_file_name[cnt_1][0]]))
        temp_fl.close
    dict_aofc.clear()
    #
    # Extract ids of mapped discordant reads.
    for cnt_1 in range( len(ref_type_file_name) ):
        line_number = 0
        #with open(dir_path+'/intermediate_files/'+ref_type_file_name[cnt_1][0]+'_ref-aligned_line-num_id.dat', 'r') as rali:
        with open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_ref-aligned_line-num_id.dat', 'r') as rali:
            rali_lines = rali.readlines()
        rali.close()
        dict_rali = {}
        for item in rali_lines:
            dict_rali[int(item.split()[0])] = item.split()[1]    

        raw_out_file_line = []
        samfile = pysam.AlignmentFile(read_bam, 'rb')
        for read in samfile.fetch():
            line_number += 1
            try:
                if dict_rali[line_number] == read.query_name:
                    raw_out_file_line.append(read.query_name +' '+ str(read.flag))
            except KeyError:
                continue
        samfile.close()
        #raw_out_file = open(dir_path+'/intermediate_files/'+ref_type_file_name[cnt_1][0]+'_read-bam_id_flag.dat','w')
        raw_out_file = open(preprocess_dir_realpath +'/'+ref_type_file_name[cnt_1][0]+'_read-bam_id_flag.dat','w')
        raw_out_file.write('\n'.join(raw_out_file_line))
        del raw_out_file_line[:]
        raw_out_file.close()
        dict_rali.clear()

#    Extract ids of uniquly mapped clipped reads
    for cnt_1 in range( len(ref_type_file_name) ):
        line_number = 0
        #with open(dir_path+'/intermediate_files/'+ref_type_file_name[cnt_1][0]+'_clipped_ref-aligned_line-num_id.dat', 'r') as crali:
        with open(preprocess_dir_realpath +'/'+ref_type_file_name[cnt_1][0]+'_clipped_ref-aligned_line-num_id.dat', 'r') as crali:
            crali_lines = crali.readlines()
        crali.close()
        dict_crali = {}
        for item in crali_lines:
            dict_crali[int(item.split()[0])] = item.split()[1]

        raw_out_file_line = []
        samfile = pysam.AlignmentFile(read_bam_clipped, 'rb')
        for read in samfile.fetch():
            line_number += 1
            try:
                if (dict_crali[line_number] == read.query_name) and (read.is_supplementary != True) and \
                    (read.cigartuples != None) and (read.has_tag('XA') or read.has_tag('SA') or \
                    read.mapping_quality >= min_mapq) and (read.mapping_quality >= min_mapq_uniq):
        
                    write_flag = check_uniq_mapping( read, args )

                    if write_flag == 'y':
                        #
                        clipped_side = 'X'
                        if  ( read.cigartuples[-1][0] == 4 ) and ( read.cigartuples[-1][1] > clipped_length ):
                            clipped_side = 'R'
                        elif ( read.cigartuples[0][0] == 4 ) and ( read.cigartuples[0][1] > clipped_length ):
                            clipped_side = 'L'
                        #
                        raw_out_file_line.append(read.query_name +' '+ str(read.flag) + ' ' \
                            + str(read.reference_name) + ' ' + str(read.reference_start) + ' ' + str(read.reference_end) + ' ' + clipped_side)
            except KeyError:
                continue
        samfile.close()
        #raw_out_file = open(dir_path+'/intermediate_files/'+ref_type_file_name[cnt_1][0]+'_clipped_read-bam_id_flag.dat','w')
        raw_out_file = open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_clipped_read-bam_id_flag.dat','w')
        raw_out_file.write('\n'.join(raw_out_file_line))
        del raw_out_file_line[:]
        raw_out_file.close()
        dict_crali.clear()

    # Search for mate and their position
    samfile_idx = pysam.AlignmentFile(read_bam, 'rb')
    id_index = pysam.IndexedReads(samfile_idx)
    id_index.build()

    for cnt_1 in range( len(ref_type_file_name) ):
        ##
        #with open(dir_path+'/intermediate_files/'+ref_type_file_name[cnt_1][0]+'_read-bam_id_flag.dat' ,'r') as read_bam_dat:
        with open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_read-bam_id_flag.dat' ,'r') as read_bam_dat:
            read_bam_dat_lines = read_bam_dat.readlines()
        read_bam_dat.close

        mate_out_file_line = []
        for bam_line in read_bam_dat_lines:
            iterator = id_index.find(bam_line.split()[0])
            for read in iterator:
                if ((int(bam_line.split()[1]) & 0x40) != (read.flag & 0x40)) and (read.is_supplementary != True) and \
                    (read.cigartuples != None) and (read.has_tag('XA') or read.has_tag('SA') or \
                    read.mapping_quality >= min_mapq) and (read.mapping_quality >= min_mapq_uniq): #!= (int(bam_line.split()[1]) & 0x40):
                    
                    write_flag = check_uniq_mapping( read, args )

                    if write_flag == 'y':
                        mate_out_file_line.append(read.query_name + ' ' + str(read.flag) + ' ' + \
                        str(read.reference_name) + ' ' + str(read.reference_start) + ' ' + str(read.reference_end))
        del read_bam_dat_lines[:]
        read_bam_dat.close()

        #mate_out_file = open(dir_path+'/intermediate_files/'+ref_type_file_name[cnt_1][0]+'_read-bam_mate_id_pos.dat', 'w')
        mate_out_file = open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_read-bam_mate_id_pos.dat', 'w')
        mate_out_file.write('\n'.join(mate_out_file_line))
        del mate_out_file_line[:]
        mate_out_file.close()
    samfile_idx.close()
    
    read_positions_clusters_file = open('initial_predictions.txt', 'w')
    for cnt_1 in range( len(ref_type_file_name) ):
        flag_read_position = 'y'
        read_positions = []
        #with open(dir_path+'/intermediate_files/'+ref_type_file_name[cnt_1][0]+'_read-bam_mate_id_pos.dat', 'r') as mate_id_dat:
        with open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_read-bam_mate_id_pos.dat', 'r') as mate_id_dat:
            mate_id_dat_lines = mate_id_dat.readlines()
        mate_id_dat.close()

        for line in mate_id_dat_lines:
            read_positions.append((line.split()[2], int(line.split()[3])))
        del mate_id_dat_lines[:]

        #with open(dir_path+'/intermediate_files/'+ref_type_file_name[cnt_1][0]+'_clipped_read-bam_id_flag.dat', 'r') as read_id_dat:
        with open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_clipped_read-bam_id_flag.dat', 'r') as read_id_dat:
            read_id_dat_lines = read_id_dat.readlines() 
        read_id_dat.close()
        for line in read_id_dat_lines:
            read_positions.append((line.split()[2], int(line.split()[3])))
        del read_id_dat_lines[:]

        read_positions_sorted = sorted(read_positions, key=lambda x: (x[0], x[1]))
        read_positions_clusters = break_points_2d(read_positions_sorted, discord_rd_clust_denst, \
                        insert_size, read_length/2) #<<<<< change it to insertsize-read_length

        #Remove predictions from very high coverage areas
        samfile = pysam.AlignmentFile(bam_full, 'rb')
        read_positions_clusters_nohicov = []
        coverage_values = []
        for clust_pos in read_positions_clusters:
            try:
                for pileupcolumn in samfile.pileup(clust_pos[0], int(clust_pos[1])-(insert_size+read_length), \
                                            int(clust_pos[1])+(insert_size+read_length)):    
                    coverage_values.append(int(pileupcolumn.n))            
                coverage_values_np = np.array(coverage_values)
                print(np.mean(coverage_values_np[read_length:-1*read_length]))
                if np.mean(coverage_values_np[read_length:-1*read_length]) < coverage_cutoff: #cct
                    read_positions_clusters_nohicov.append(clust_pos)
                del coverage_values[:]
                del coverage_values_np
            except ValueError:
                read_positions_clusters_nohicov.append(clust_pos)
        samfile.close()

        #Write final cluster to a file    
        for info in read_positions_clusters_nohicov:
            read_positions_clusters_file.write(ref_type_file_name[cnt_1][0]+'\t'+'\t'.join(str(dat) for dat in info[:-1])+'\n') 

        del read_positions[:]
        del read_positions_sorted[:]
        del read_positions_clusters[:]
        del read_positions_clusters_nohicov[:]

    read_positions_clusters_file.close()
#--------------------------------------------------------------
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

#--------------------------------------------------------------
def pat_check(inp_seq, query_len, mis_match):
    #
    inp_seq = inp_seq.upper()
    #
    match_len = query_len - mis_match
    #
    polyAT_test_flag = 0
    polyAT_type = 'X'
    # 1. Check poly-A
    start_idx  = 0
    end_idx = query_len
    while ( end_idx <= len(inp_seq) ):
        #
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
        #
    return( polyAT_test_flag, polyAT_type )
#--------------------------------------------------------------
def exec_nadiscover(args):
    #
    dir_path = os.getcwd() #os.path.dirname(os.path.realpath(__file__))
    sys.stdout.write('Working directory: '+str(dir_path)+'\n')
    #
    preprocess_dir_realpath = os.path.realpath(args.preprocess_dir)
    sys.stdout.write('preprocessing/intermediate file directory: '+ preprocess_dir_realpath +'\n')
    #
    fofn_ref_realpath = os.path.realpath(args.fofn_ref)
    sys.stdout.write('fofn_ref file: '+ fofn_ref_realpath +'\n')
    #
    bam_full = args.bam_inp
    sys.stdout.write('bam file name: '+str(bam_full)+'\n')
    #
    bam_short_name = bam_full.split('/')[-1][:-4]
    sys.stdout.write('bam short name: '+ bam_short_name +'\n')
    #
    insert_size = args.isz_inp
    sys.stdout.write('Insert size estimate: '+str(insert_size)+'\n')
    #
    discord_rd_clust_denst = args.drd_inp
    sys.stdout.write('Number of reads in a cluster to call it insertion: '+str(discord_rd_clust_denst)+'\n')
    #
    read_length = args.rdl_inp
    sys.stdout.write('Average read length: '+str(read_length)+'\n')
    #
    coverage_cutoff = args.cct_inp
    sys.stdout.write('Coverage cutoff to skip a region: '+str(coverage_cutoff)+'\n')
    #
    #fofn_ref_file = args.fofn_ref
    #sys.stdout.write('FoFn reference sequence: '+str(fofn_ref_file)+'\n\n')
    #
    min_mapq = args.mpq_inp
    sys.stdout.write('minimum mapping quality: '+str(min_mapq)+'\n')
    #
    min_mapq_uniq = args.mpqu_inp
    sys.stdout.write('minimum mapping quality: '+str(min_mapq_uniq)+'\n')
    #
    clipped_length = args.cll_inp
    sys.stdout.write('Minimum clipped length: '+str(clipped_length)+'\n')
    #
    pat_query_len = args.pql_inp
    
    #
    pat_mis_match = args.pmm_inp

    #
    rmsk_bed = args.rmsk_bed
    #
    #read_bam_clipped = dir_path+'/preprocessed_files/'+bam_full.split('/')[-1][:-4]+'_clipped.bam'
    #read_bam_clipped = args.preprocess_dir +'/'+bam_full.split('/')[-1][:-4]+'_clipped.bam'
    read_bam_clipped = preprocess_dir_realpath +'/'+bam_short_name+'_clipped.bam'
    #
    #
    
    if args.pat_inp:
        #
        pat_out_file_lines = []
        #
        samfile_clipped = pysam.AlignmentFile(read_bam_clipped, "rb")
        #    
        for read in samfile_clipped.fetch():
            if (read.is_supplementary != True) and ( read.has_tag('XA') or read.has_tag('SA') ) \
                and ( read.mapping_quality >= min_mapq_uniq ):
                #
                write_flag = check_uniq_mapping( read, args )
                #
                if write_flag == 'y':
                    clipped_side = 'X'
                    query_sequence = 'tmp'
                    if ( read.cigartuples[-1][0] == 4 ) and ( read.cigartuples[-1][1] > clipped_length ) and \
                        ( (read.cigartuples[0][0] == 4 and read.cigartuples[0][1] > 5) != True ):
                        #    
                        clipped_side = 'R'
                        query_sequence=read.query_sequence[read.infer_query_length() \
                                    - read.cigartuples[-1][1]-1:read.infer_query_length()]
                        #
                    elif ( read.cigartuples[0][0] == 4 ) and ( read.cigartuples[0][1] > clipped_length ) and \
                        ( (read.cigartuples[-1][0] == 4 and read.cigartuples[-1][1] > 5) != True ):
                        #
                        clipped_side = 'L'
                        query_sequence=read.query_sequence[0:read.cigartuples[0][1]-1]
                        #
                    pat_flag, pat_type = pat_check(query_sequence, pat_query_len, pat_mis_match)
                    if pat_flag == 1:
                        #
                        pat_out_file_lines.append( read.query_name +' '+ str(read.flag) + ' ' \
                            + str(read.reference_name) + ' ' + str(read.reference_start) + ' ' + \
                                         str(read.reference_end) + ' ' + clipped_side )
                            # +'\t'+query_sequence+'\t'+pat_type )
                        #pat_dat_lines.append( (read.reference_name, read.reference_start) )

        samfile_clipped.close()
        #
        #    
        #pat_out_file = open(dir_path+'/intermediate_files/pat_clipped_read-bam_id_flag.dat','w')
        pat_out_file = open(preprocess_dir_realpath+'/pat_clipped_read-bam_id_flag.dat','w')
        pat_out_file.write('\n'.join(pat_out_file_lines))
        #del pat_out_file_lines[:]
        pat_out_file.close()
        #
        #
    #sys.exit()        
    #
    ref_type_file_name = []
    #with open(fofn_ref_file, 'r') as ref_type_file_file:
    with open(fofn_ref_realpath, 'r') as ref_type_file_file:
        for line in ref_type_file_file:
            ref_type_file_name.append([line.split()[0], line.split()[1]])
    #
    #
    if args.nas_inp:
        # Process RMSK file.
        dict_ref_bed = {}
        for items in ref_type_file_name:
            dict_ref_bed[items[0]] = []
        #
        with open(rmsk_bed, 'r') as rmsk_bed_file:
            rmsk_bed_file_lines = rmsk_bed_file.readlines()
        rmsk_bed_file.close()
        #
        # Consolidate rmsk bed items.
        rmsk_bed_items = []
        prev_chr, prev_start, prev_end, prev_type = 'chr00', 0, 0, 'none'
        for line in rmsk_bed_file_lines:
            items = line.strip().split()
            if (items[0] == prev_chr) and (items[3] == prev_type) and (int(items[1]) <= prev_end):
                prev_end = int(items[2])
                continue
            else:
                rmsk_bed_items.append( prev_chr +'\t'+ str(prev_start) +'\t'+ str(prev_end) +'\t'+ prev_type +'\n' )
                #rmsk_bed_items.append( items[0] +'\t'+ items[1] +'\t'+ items[2] +'\t'+ items[3] +'\n' )
                prev_chr, prev_start, prev_end, prev_type = items[0], int(items[1]), int(items[2]), items[3]
        #open(dir_path+'/intermediate_files/rmsk_cons.bed', 'w').write(''.join(rmsk_bed_items))
        open(preprocess_dir_realpath+'/rmsk_cons.bed', 'w').write(''.join(rmsk_bed_items))
        #
        # Indexing
        dict_bed = {}
        chr_track = 'chr00'
        for line in rmsk_bed_items:
            if line.startswith('chr00'):
                continue
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
        #
        discord_pos_list = []
        #discord_bam = dir_path+'/preprocessed_files/'+bam_full.split('/')[-1][:-4]+'_discord.bam'
        #discord_bam = args.preprocess_dir + '/'+bam_full.split('/')[-1][:-4]+'_discord.bam'
        discord_bam = preprocess_dir_realpath + '/'+bam_short_name+'_discord.bam'
        sys.stdout.write('Discordant bam file: '+discord_bam+'\n')
        #
        samfile_discord = pysam.AlignmentFile(discord_bam, 'rb')
        #
        for read in samfile_discord.fetch():
            discord_pos_list.append( read.query_name +'\t'+ str(read.flag) +'\t'+ read.reference_name +'\t'+ str(read.reference_start) +'\t'+ str(read.reference_end)+'\n')
        #
        samfile_discord.close()
        #
        #open(dir_path+'/intermediate_files/discord_pos_list.txt', 'w').write(''.join(discord_pos_list))
        open(preprocess_dir_realpath+'/discord_pos_list.txt', 'w').write(''.join(discord_pos_list))
        #
        if args.flg_all:
            clipped_pos_list = []
            #clipped_bam = dir_path+'/preprocessed_files/'+bam_full.split('/')[-1][:-4]+'_clipped.bam'
            #clipped_bam = args.preprocess_dir + '/'+bam_full.split('/')[-1][:-4]+'_clipped.bam'
            clipped_bam = preprocess_dir_realpath + '/'+bam_short_name+'_clipped.bam'
            sys.stdout.write('Clipped bam file: '+clipped_bam+'\n')
            #
            samfile_clipped = pysam.AlignmentFile(clipped_bam, 'rb')
            for read in samfile_clipped.fetch():
                #
                if (read.is_supplementary != True) and ( read.has_tag('XA') or read.has_tag('SA') ) \
                    and ( read.mapping_quality >= min_mapq_uniq ):
                    #
                    write_flag = check_uniq_mapping( read, args )
                    #
                    if write_flag == 'y':
                        chrom, clp_start, clp_end = alt_mapped_pos( read, args )
                        clipped_pos_list.append( read.query_name +'\t'+ str(read.flag) +'\t'+ read.reference_name +'\t'+ \
                            str(read.reference_start) +'\t'+ str(read.reference_end) +'\t'+ chrom +'\t'+ str(clp_start) +'\t'+ str(clp_end) +'\n' )

            #open(dir_path+'/intermediate_files/clipped_pos_list.txt', 'w').write(''.join(clipped_pos_list))
            open(preprocess_dir_realpath+'/clipped_pos_list.txt', 'w').write(''.join(clipped_pos_list))
            samfile_clipped.close()
        #
        #
        dict_discord = {}
        for line in discord_pos_list:
            items = line.strip().split()
            try:
                for info in dict_bed[items[2]][(dict_bed[items[2]+'_idx'].get(int((int(items[3]))/10000),1))-1:]:
                    if len( range(max(int(items[3]), int(info[0])), min(int(items[4]), int(info[1]))+1) ) >=25: # atleast 25bp overlap
                        try:
                            dict_discord[info[2]].append(items[0] +'\t'+ items[1] +'\n')

                        except KeyError:
                            dict_discord[info[2]] = []
                            dict_discord[info[2]].append(items[0] +'\t'+ items[1] +'\n')
                        break
                    elif (int(items[3]) > int(info[0])):
                        break

            except (ValueError, KeyError):
                continue
        #
        for cnt_1 in range(0, len(ref_type_file_name)):
            try:
                #open(dir_path+'/intermediate_files/'+ref_type_file_name[cnt_1][0]+'_noalign_read-bam_id_flag.dat', 'w').write(''.join(dict_discord[ref_type_file_name[cnt_1][0]]))
                open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_noalign_read-bam_id_flag.dat', 'w').write(''.join(dict_discord[ref_type_file_name[cnt_1][0]]))
            except KeyError:
                #open(dir_path+'/intermediate_files/'+ref_type_file_name[cnt_1][0]+'_noalign_read-bam_id_flag.dat', 'w').write('')
                open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_noalign_read-bam_id_flag.dat', 'w').write('')
        #
        #
        if args.flg_all:
            dict_clipped = {}
            for line in clipped_pos_list:
                items = line.strip().split()
                try:
                    for info in dict_bed[items[5]][(dict_bed[items[5]+'_idx'].get(int((int(items[6]))/10000),1))-1:]:
                        if len( range(max(int(items[6]), int(info[0])), min(int(items[7]), int(info[1]))+1) ) >=20: # atleast 25bp overlap
                            try:
                                dict_clipped[info[2]].append( line )

                            except KeyError:
                                dict_clipped[info[2]] = []
                                dict_clipped[info[2]].append( line )
                            break
                        elif (int(items[6]) > int(info[0])):
                            break

                except (ValueError, KeyError):
                    continue
            #
            for cnt_1 in range(0, len(ref_type_file_name)):
                try:
                    #open(dir_path+'/intermediate_files/'+ref_type_file_name[cnt_1][0]+'_noalign_clipped_read-bam_id_flag.dat', 'w').write(''.join(dict_clipped[ref_type_file_name[cnt_1][0]]))
                    open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_noalign_clipped_read-bam_id_flag.dat', 'w').write(''.join(dict_clipped[ref_type_file_name[cnt_1][0]]))
                except KeyError:
                    #open(dir_path+'/intermediate_files/'+ref_type_file_name[cnt_1][0]+'_noalign_clipped_read-bam_id_flag.dat', 'w').write('')
                    open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_noalign_clipped_read-bam_id_flag.dat', 'w').write('')
            #
        #
        # Search for mate and their position
        samfile_idx = pysam.AlignmentFile(discord_bam, 'rb')
        id_index = pysam.IndexedReads(samfile_idx)
        id_index.build()
        #
        for cnt_1 in range( len(ref_type_file_name) ):
            ##
            #with open(dir_path+'/intermediate_files/'+ref_type_file_name[cnt_1][0]+'_noalign_read-bam_id_flag.dat' ,'r') as read_bam_dat:
            with open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_noalign_read-bam_id_flag.dat' ,'r') as read_bam_dat:
                read_bam_dat_lines = read_bam_dat.readlines()
            read_bam_dat.close
            #
            mate_out_file_line = []
            for bam_line in read_bam_dat_lines:
                iterator = id_index.find(bam_line.split()[0])
                for read in iterator:
                    if ((int(bam_line.split()[1]) & 0x40) != (read.flag & 0x40)) and (read.is_supplementary != True) and \
                        (read.cigartuples != None) and (read.has_tag('XA') or read.has_tag('SA') or \
                        read.mapping_quality >= min_mapq) and (read.mapping_quality >= min_mapq_uniq): #!= (int(bam_line.split()[1]) & 0x40):

                        write_flag = check_uniq_mapping( read, args )

                        if write_flag == 'y':
                            mate_out_file_line.append(read.query_name + ' ' + str(read.flag) + ' ' + \
                            str(read.reference_name) + ' ' + str(read.reference_start) + ' ' + str(read.reference_end))
            del read_bam_dat_lines[:]
            read_bam_dat.close()
            #
            #mate_out_file = open(dir_path+'/intermediate_files/'+ref_type_file_name[cnt_1][0]+'_noalign_read-bam_mate_id_pos.dat', 'w')
            mate_out_file = open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_noalign_read-bam_mate_id_pos.dat', 'w')
            mate_out_file.write('\n'.join(mate_out_file_line))
            del mate_out_file_line[:]
            mate_out_file.close()
        samfile_idx.close()
        #
        # nas condition done
        #
    #read_positions_clusters_file = open('initial_predictions_noalign.txt', 'w')
    read_positions_clusters_file = open(args.ouput_file, 'w')
    #
    for cnt_1 in range( len(ref_type_file_name) ):
        #
        read_positions = []
        dict_pat_test = {}
        #
        # nas condition
        if args.nas_inp:
            #
            #with open(dir_path+'/intermediate_files/'+ref_type_file_name[cnt_1][0]+'_noalign_read-bam_mate_id_pos.dat', 'r') as noalign_mate_id_dat:
            with open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_noalign_read-bam_mate_id_pos.dat', 'r') as noalign_mate_id_dat:
                noalign_mate_id_dat_lines = noalign_mate_id_dat.readlines()
            noalign_mate_id_dat.close()
            #
            if args.merged:
                #
                #with open(dir_path+'/intermediate_files/'+ref_type_file_name[cnt_1][0]+'_read-bam_mate_id_pos.dat', 'r') as align_mate_id_dat:
                with open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_read-bam_mate_id_pos.dat', 'r') as align_mate_id_dat:
                    align_mate_id_dat_lines = align_mate_id_dat.readlines()
                align_mate_id_dat.close()
                #
                dict_align = {}
                align_mate_read_positions = []
                #
                for line in align_mate_id_dat_lines:
                    #
                    align_mate_read_positions.append( ( line.split()[2], int(line.split()[3]) ) )
                    #
                    item = line.strip().split()[0]+'_'+line.strip().split()[1]+'_'+line.strip().split()[2]+'_'+line.strip().split()[3]
                    dict_align[item] = 1
                    #
                del align_mate_id_dat_lines[:]
                #
                noalign_mate_read_positions = []
                #
                for line in noalign_mate_id_dat_lines:
                    #
                    test_item = line.strip().split()[0]+'_'+line.strip().split()[1]+'_'+line.strip().split()[2]+'_'+line.strip().split()[3]
                    #
                    try:
                        if dict_align[test_item] == 1:
                            continue
                    except KeyError:
                        noalign_mate_read_positions.append( ( line.split()[2], int(line.split()[3]) ) )
                    #
                # del noalign_mate_id_dat_lines[:] need this if no merged
                #
                mate_read_positions = align_mate_read_positions + noalign_mate_read_positions
                #
                #print(mate_read_positions)
                #print(align_mate_read_positions)
                #print(noalign_mate_read_positions)
                #
                del align_mate_read_positions[:]
                del noalign_mate_read_positions[:]
                dict_align.clear()
                #
                # Clipped
                #
                #with open(dir_path+'/intermediate_files/'+ref_type_file_name[cnt_1][0]+'_clipped_read-bam_id_flag.dat', 'r') as align_id_dat:
                with open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_clipped_read-bam_id_flag.dat', 'r') as align_id_dat:
                    align_id_dat_lines = align_id_dat.readlines()
                align_id_dat.close()
                #
                dict_clp_align = {}
                align_read_positions = []
                #
                for line in align_id_dat_lines:
                    #
                    align_read_positions.append( ( line.split()[2], int(line.split()[3]) ) )
                    #
                    item = line.strip().split()[0]+'_'+line.strip().split()[1]+'_'+line.strip().split()[2]+'_'+line.strip().split()[3]
                    #
                    dict_clp_align[item] = 1
                    dict_pat_test[line.strip().split()[0]] = line.strip().split()[1]
                    #
                del align_id_dat_lines[:]

                noalign_read_positions = []
                if args.flg_all:
                    #
                    #with open(dir_path+'/intermediate_files/'+ref_type_file_name[cnt_1][0]+'_noalign_clipped_read-bam_id_flag.dat', 'r') as noalign_id_dat:
                    with open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_noalign_clipped_read-bam_id_flag.dat', 'r') as noalign_id_dat:
                        noalign_id_dat_lines = noalign_id_dat.readlines()
                    noalign_id_dat.close()
                    #
                    for line in noalign_id_dat_lines:
                        #
                        test_item = line.strip().split()[0]+'_'+line.strip().split()[1]+'_'+line.strip().split()[2]+'_'+line.strip().split()[3]
                        #
                        try:
                            if dict_clp_align[test_item] == 1:
                                continue
                        except KeyError:
                            noalign_read_positions.append( ( line.split()[2], int(line.split()[3]) ) )
                            dict_pat_test[line.strip().split()[0]] = line.strip().split()[1]
                        #
                    del noalign_id_dat_lines[:]
                    #
                clpd_read_positions = align_read_positions + noalign_read_positions
                #
                del align_read_positions[:]
                #
                del noalign_read_positions[:]
                #
                dict_clp_align.clear()
                #
                #
                read_positions = mate_read_positions + clpd_read_positions
                del mate_read_positions[:]
                del clpd_read_positions[:]
                #
            if not args.merged:
                #
                #
                mate_read_positions = []
                for line in noalign_mate_id_dat_lines:
                    mate_read_positions.append( ( line.split()[2], int(line.split()[3]) ) )

                clpd_read_positions = []
                if args.flg_all:
                    #
                    #with open(dir_path+'/intermediate_files/'+ref_type_file_name[cnt_1][0]+'_noalign_clipped_read-bam_id_flag.dat', 'r') as noalign_id_dat:
                    with open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_noalign_clipped_read-bam_id_flag.dat', 'r') as noalign_id_dat:
                        noalign_id_dat_lines = noalign_id_dat.readlines()
                    noalign_id_dat.close()
                    #
                    #
                    for line in noalign_id_dat_lines:
                        clpd_read_positions.append( ( line.strip().split()[2], int(line.strip().split()[3]) ) )
                        dict_pat_test[line.strip().split()[0]] = line.strip().split()[1]
                    del noalign_id_dat_lines[:]
                    #
                read_positions = mate_read_positions + clpd_read_positions
                del mate_read_positions[:]
                del clpd_read_positions[:]
                #
            del noalign_mate_id_dat_lines[:]
            # nas consition ends
            #
        if not args.nas_inp:
            print(args.nas_inp)
            #    
            #with open(dir_path+'/intermediate_files/'+ref_type_file_name[cnt_1][0]+'_read-bam_mate_id_pos.dat', 'r') as mate_id_dat:
            with open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_read-bam_mate_id_pos.dat', 'r') as mate_id_dat:
                mate_id_dat_lines = mate_id_dat.readlines()
            mate_id_dat.close()

            for line in mate_id_dat_lines:
                read_positions.append((line.strip().split()[2], int(line.strip().split()[3])))
            del mate_id_dat_lines[:]
            #
            #with open(dir_path+'/intermediate_files/'+ref_type_file_name[cnt_1][0]+'_clipped_read-bam_id_flag.dat', 'r') as clpd_read_id_dat:
            with open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_clipped_read-bam_id_flag.dat', 'r') as clpd_read_id_dat:
                clpd_read_id_dat_lines = clpd_read_id_dat.readlines()
            clpd_read_id_dat.close()
            #
            for line in clpd_read_id_dat_lines:
                read_positions.append((line.strip().split()[2], int(line.strip().split()[3])))
                dict_pat_test[line.strip().split()[0]] = line.strip().split()[1]            
    
            del clpd_read_id_dat_lines[:]
    #
    #
        read_positions_sorted = sorted(read_positions, key=lambda x: (x[0], x[1]))
        #
        if args.pat_inp:
            #
            #
            #if args.pat_inp: # polyA/T search
            min_non_pat_rd_denst = 3 #min_non_pat_rd_denst
            #
            pat_dat_lines = []
            for line in pat_out_file_lines:
                try:
                    if dict_pat_test[ line.strip().split()[0] ] ==  line.strip().split()[1]:
                        continue
                except (KeyError, ValueError):
                    pat_dat_lines.append( ( line.strip().split()[2], int(line.strip().split()[3]) ) )
            #
            dict_pat_test.clear()
            pat_dat_lines_sorted = sorted(pat_dat_lines, key=lambda x: (x[0], x[1]))
            #
            # Indexing
            dict_pat = {}
            chr_track = 'chr00'
            for items in pat_dat_lines_sorted:
                if items[0] != chr_track:
                    line_num = 0
                    chr_track = items[0]
                    dict_pat[chr_track] = []
                    dict_pat[chr_track+'_idx'] = {}
                    pos_track = int( int(items[1])/10000 )
                    dict_pat[chr_track+'_idx'][pos_track] = 1 #line_num to handle -1 for very first line
                    dict_pat[chr_track].append( int(items[1]) )
                    continue
                dict_pat[chr_track].append( int(items[1]) )
                #
                line_num += 1
                if pos_track != int( int(items[1]) /10000 ):
                    for i in range(pos_track, int( int(items[1])/10000)):
                        if pos_track + 1 == int( int(items[1])/10000):
                            dict_pat[chr_track+'_idx'][ int( int(items[1])/10000) ] = line_num
                            pos_track = int( int(items[1])/10000 )
                        if pos_track + 1 < int( int(items[1])/10000):
                            dict_pat[chr_track+'_idx'][pos_track+1] = dict_pat[chr_track+'_idx'][pos_track]
                            pos_track = pos_track + 1

        #    
            read_positions_clusters = break_points_2d(read_positions_sorted, min_non_pat_rd_denst, \
                        insert_size, read_length/2) #<<<<< change it to insertsize-read_length
            #
            #    if args.pat_inp:
            #
            read_positions_clusters_pat = []
            for clust_pos in read_positions_clusters:
                if clust_pos[2] >= discord_rd_clust_denst:
                    read_positions_clusters_pat.append( clust_pos )
                    continue
                # append pat data
                test_lower_lim = clust_pos[3][0] - insert_size 
                test_upper_lim = clust_pos[3][-1] + insert_size
                try:
                    for info in dict_pat[ clust_pos[0] ][ ( dict_pat[ clust_pos[0]+'_idx' ].get( int( clust_pos[1]/10000 ),1 ) )-1: ]:
                        if info in range( test_lower_lim, test_upper_lim ):
                            clust_pos[3].append( info )    
                except KeyError:
                    continue
                # break point
                tmp_position_sorted = sorted( clust_pos[3], key=lambda x:x )
                tmp_position_clstrs = break_points( tmp_position_sorted, discord_rd_clust_denst, insert_size )
                if not tmp_position_clstrs:
                    continue    # not enough reads in cluster        
                # closest is within isz range
                print( min(tmp_position_clstrs, key=lambda x:abs(x-clust_pos[1])) )
                if abs(clust_pos[1] - (min(tmp_position_clstrs, key=lambda x:abs(x-clust_pos[1])))) < insert_size:
                    clust_pos[1] = min(tmp_position_clstrs, key=lambda x:abs(x-clust_pos[1]))
                    read_positions_clusters_pat.append( clust_pos )
                #
                del tmp_position_clstrs[:]
                #
            del read_positions_clusters[:]
            #
            read_positions_clusters = read_positions_clusters_pat.copy()
            #
            del read_positions_clusters_pat[:]
            #
            dict_pat.clear()
            #
        if not args.pat_inp:
            read_positions_clusters = break_points_2d(read_positions_sorted, discord_rd_clust_denst, \
                        insert_size, read_length/2)    
        #Remove predictions from very high coverage areas
        samfile = pysam.AlignmentFile(bam_full, 'rb')
        read_positions_clusters_nohicov = []
        coverage_values = []
        for clust_pos in read_positions_clusters:
            try:
                for pileupcolumn in samfile.pileup(clust_pos[0], int(clust_pos[1])-(insert_size+read_length), \
                                            int(clust_pos[1])+(insert_size+read_length)):
                    coverage_values.append(int(pileupcolumn.n))
                coverage_values_np = np.array(coverage_values)
                if np.mean(coverage_values_np[read_length:-1*read_length]) < coverage_cutoff: #cct
                    read_positions_clusters_nohicov.append(clust_pos)
                del coverage_values[:]
                del coverage_values_np
            except ValueError:
                read_positions_clusters_nohicov.append(clust_pos)
        samfile.close()

        #Write final cluster to a file
        for info in read_positions_clusters_nohicov:
            read_positions_clusters_file.write(ref_type_file_name[cnt_1][0]+'\t'+'\t'.join(str(dat) \
                                    for dat in info[:-1])+'\t'+str(len(info[-1]))+'\n')

        del read_positions[:]
        del read_positions_sorted[:]
        del read_positions_clusters[:]
        del read_positions_clusters_nohicov[:]

        read_positions_clusters_file.close()
#--------------------------------------------------------------
def exec_cluster2d(args):    
    #
    preprocess_dir_realpath = os.path.realpath(args.preprocess_dir)
    sys.stdout.write('preprocessing/intermediate file directory: '+ preprocess_dir_realpath +'\n')
    #
    fofn_ref_realpath = os.path.realpath(args.fofn_ref)
    sys.stdout.write('fofn_ref file: '+ fofn_ref_realpath +'\n')
    #
    dir_path = os.getcwd()
    sys.stdout.write('working directory: '+ dir_path +'\n')
    #
    bam_full = args.bam_inp
    insert_size = args.isz_inp
    discord_rd_clust_denst = args.drd_inp
    read_length = args.rdl_inp
    coverage_cutoff    = args.cct_inp
    #fofn_ref_file = args.fofn_ref
    ref_type_file_name = []
    #with open(fofn_ref_file, 'r') as ref_type_file_file:
    with open(fofn_ref_realpath, 'r') as ref_type_file_file:
        for line in ref_type_file_file:
            ref_type_file_name.append([line.split()[0], line.split()[1]])

    read_positions_clusters_file = open('recluster_initial_predictions.txt', 'w')
    for cnt_1 in range(0, len(ref_type_file_name)):
        flag_read_position = 'y'
        read_positions = []
        #
        if args.flg_all:
            #with open(dir_path+'/intermediate_files/'+ref_type_file_name[cnt_1][0]+'_read-bam_mate_id_pos.dat', 'r') as mate_id_dat:
            with open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_read-bam_mate_id_pos.dat', 'r') as mate_id_dat:
                mate_id_dat_lines = mate_id_dat.readlines()
            mate_id_dat.close()

            for line in mate_id_dat_lines:
                read_positions.append((line.split()[2], int(line.split()[3])))
            del mate_id_dat_lines[:]
        #
        #with open(dir_path+'/intermediate_files/'+ref_type_file_name[cnt_1][0]+'_clipped_read-bam_id_flag.dat', 'r') as read_id_dat:
        with open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_clipped_read-bam_id_flag.dat', 'r') as read_id_dat:
            read_id_dat_lines = read_id_dat.readlines()
        read_id_dat.close()
        for line in read_id_dat_lines:
            read_positions.append((line.split()[2], int(line.split()[3])))
        del read_id_dat_lines[:]
        #
        read_positions_sorted = sorted(read_positions, key=lambda x: (x[0], x[1]))
        read_positions_clusters = break_points_2d(read_positions_sorted, discord_rd_clust_denst, insert_size, read_length/2) #<<<<< change it to insertsize-read_length
        #
        #Remove predictions from very high coverage areas
        samfile = pysam.AlignmentFile(bam_full, 'rb')
        read_positions_clusters_nohicov = []
        coverage_values = []
        for clust_pos in read_positions_clusters:
            try:
                for pileupcolumn in samfile.pileup(clust_pos[0], int(clust_pos[1])-(insert_size+read_length), int(clust_pos[1])+(insert_size+read_length)):
                    coverage_values.append(int(pileupcolumn.n))
                coverage_values_np = np.array(coverage_values)
                if np.mean(coverage_values_np[read_length:-1*read_length]) < coverage_cutoff: #cct
                    read_positions_clusters_nohicov.append(clust_pos)
                del coverage_values[:]
                del coverage_values_np
            except ValueError:
                read_positions_clusters_nohicov.append(clust_pos)
        samfile.close()

        #Write final cluster to a file
        for info in read_positions_clusters_nohicov:
            read_positions_clusters_file.write(ref_type_file_name[cnt_1][0]+'\t'+'\t'.join(str(dat) for dat in info[:-1])+'\n')

        del read_positions[:]
        del read_positions_sorted[:]
        del read_positions_clusters[:]
        del read_positions_clusters_nohicov[:]
    read_positions_clusters_file.close()
#--------------------------------------------------------------
def exec_filter(args):
    #
    preprocess_dir_realpath = os.path.realpath(args.preprocess_dir)
    sys.stdout.write('preprocessing/intermediate file directory: '+ preprocess_dir_realpath +'\n')
    #
    #fofn_ref_realpath = os.path.realpath(args.fofn_ref)
    #sys.stdout.write('fofn_ref file: '+ fofn_ref_realpath +'\n')
    #
    dir_path = os.getcwd()
    sys.stdout.write('working directory: '+ dir_path +'\n')
    #
    fofn_bed_file = os.path.realpath(args.fofn_bed)
    sys.stdout.write('fofn_bed file: '+ fofn_bed_file +'\n')
    #
    input_file_name = args.ofa_inp
    sys.stdout.write('input file name: '+ input_file_name +'\n')
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

        sys.stdout.write(filter_result+'\t'+te_loc+'\t'+line)

        left_clipped_rd, right_clipped_rd, left_discord_rd, right_discord_rd, test_class_score = 0, 0, 0, 0, 0
##    
#--------------------------------------------------------------
def exec_filter_p(args):
    #
    preprocess_dir_realpath = os.path.realpath(args.preprocess_dir)
    sys.stdout.write('preprocessing/intermediate file directory: '+ preprocess_dir_realpath +'\n')
    #
    #fofn_ref_realpath = os.path.realpath(args.fofn_ref)
    #sys.stdout.write('fofn_ref file: '+ fofn_ref_realpath +'\n')
    #
    dir_path = os.getcwd()
    sys.stdout.write('working directory: '+ dir_path +'\n')
    #
    input_file_name = args.ofa_inp
    sys.stdout.write('input file: '+ input_file_name +'\n')
    #
    fofn_bed_file = os.path.realpath(args.fofn_bed)
    sys.stdout.write('fofn_bed file: '+ fofn_bed_file +'\n')
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

        sys.stdout.write(filter_result+'\t'+te_loc+'\t'+line) #+'\t'+str(cnt_rd))

        left_clipped_rd, right_clipped_rd, left_discord_rd, right_discord_rd, test_class_score = 0, 0, 0, 0, 0

#--------------------------------------------------------------

def analysis_pat_check(    inp_fa_file, inp_fa_map_file ):
    #
    with open(inp_fa_file, 'r') as fa_file:
        fa_file_lines = fa_file.readlines()
    fa_file.close()
    #
    mapeed_seq_idx = []
    try:
        with open(inp_fa_map_file, 'r') as map_file:
            map_file_lines = map_file.readlines()
        map_file.close()
        #
        for line in  map_file_lines:
            mapeed_seq_idx.append( line.strip().split()[0] )
        #
    except FileNotFoundError:
        pass
    #
    fa_file_lines_itr = iter( fa_file_lines )
    #
    pat_test_out_lines = []
    #
    for line in fa_file_lines_itr:
        #
        info_line = line
        seq_line = next( fa_file_lines_itr )
        #
        seq_idx = (info_line.strip().split()[0]).strip('>')
        #
        if seq_idx in mapeed_seq_idx:
            continue
        pat_flag, pat_type = pat_check( seq_line, 9, 1 )
        if pat_flag == 1:
            pat_test_out_lines.append([seq_idx, pat_type ])
        #
    return ( pat_test_out_lines )

#---------------------------------------------------------------

def flt_discord( discord_mate_file ):
    #
    with open(discord_mate_file, 'r') as inp_file:
        inp_file_lines = inp_file.readlines()
    inp_file.close()
    #
    inp_file_lines_itr = iter( inp_file_lines )
    #
    dict_info = {}
    for line in inp_file_lines_itr:
        #
        info_line = line
        seq_id = line.strip().split()[1]
        seq_line = next( inp_file_lines_itr )
        seq = Seq( seq_line.strip() )
        #
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

#--------------------------------------------------------------
def find_clipped_ends( iterator_reads, insert_guess, args ):

    array_p_ps = []
    array_n_ps = []
    #
    insert_size = args.isz_inp
    clipped_length = args.cll_inp
    anchor_length = args.ahl_inp
    insert_range = args.rdl_inp
    min_mapq = args.mpq_inp
    #
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
    #    
    return(array_p_s, array_n_s)
#--------------------------------------------------------------
def exec_analyze(args):
    #
    preprocess_dir_realpath = os.path.realpath(args.preprocess_dir)
    sys.stdout.write('preprocessing/intermediate file directory: '+ preprocess_dir_realpath +'\n')
    #
    fofn_ref_realpath = os.path.realpath(args.fofn_ref)
    sys.stdout.write('preprocessing/intermediate file directory: '+ preprocess_dir_realpath +'\n')
    #
    dir_path = os.getcwd()
    sys.stdout.write('working directory: '+ dir_path +'\n')
    #up_dir_path = os.path.dirname(os.path.realpath(dir_path))
    #
    bam_full = args.bam_inp
    sys.stdout.write('Input bam file: '+bam_full+'\n')
    #
    bam_short_name = bam_full.split('/')[-1][:-4]
    sys.stdout.write('bam short name: '+ bam_short_name +'\n')
    #
    #discord_bam = up_dir_path+'/preprocessed_files/'+bam_full.split('/')[-1][:-4]+'_discord.bam'
    #discord_bam = dir_path+'/preprocessed_files/'+bam_full.split('/')[-1][:-4]+'_discord.bam'
    #discord_bam = args.preprocess_dir + '/'+bam_full.split('/')[-1][:-4]+'_discord.bam'
    discord_bam = preprocess_dir_realpath + '/'+bam_short_name+'_discord.bam'
    sys.stdout.write('Discordant bam file: '+discord_bam+'\n')
    if not check_file(discord_bam + ".bai"):
        sys.stderr.write('cannot find index for: '+discord_bam+', making now\n')
        pysam.index(discord_bam)
    #
    #te_type_file = up_dir_path+'/preprocessed_files/te_ref_type.fa'
    #te_type_file = dir_path+'/preprocessed_files/te_ref_type.fa'
    #preprocess_dir_realpath = os.path.realpath(args.preprocess_dir)
    te_type_file = preprocess_dir_realpath + '/te_ref_type.fa'
    #
    sys.stdout.write('\n----Run parameters----\n')
    #
    insert_range = args.rdl_inp # =read_length
    sys.stdout.write('Input read length: '+str(insert_range)+'\n')    
    #
    clipped_length = args.cll_inp
    sys.stdout.write('Minimum clipped read length: '+str(clipped_length)+'\n')
    #
    anchor_length = args.ahl_inp
    sys.stdout.write('Minimum anchor read length: '+str(anchor_length)+'\n')
    #
    cliped_end_range = args.cer_inp
    sys.stdout.write('Range to put clipped read in group: '+str(cliped_end_range)+'\n')
    #
    clipped_search_interval = args.csi_inp
    sys.stdout.write('Clipped search interval: '+str(clipped_search_interval)+'\n')    
    #
    min_reads_het = args.mrh_inp
    sys.stdout.write('Minimum reads to call hetrozygous: '+str(min_reads_het)+'\n')
    #
    insert_size = args.isz_inp
    sys.stdout.write('Insert size estimate: '+str(insert_size)+'\n')
    #
    min_reads_ends = args.mre_inp
    sys.stdout.write('Minimum read at each end: '+str(min_reads_ends)+'\n') # ??
    #
    qual_interval_inp = args.qii_inp
    sys.stdout.write('Mapping quality gap: '+str(qual_interval_inp)+'\n')
    #
    num_interval_inp = args.nii_inp
    sys.stdout.write('Number of quality gap searchs: '+str(num_interval_inp)+'\n')
    #
    min_mapq = args.mpq_inp
    sys.stdout.write('minimum mapping quality: '+str(min_mapq)+'\n')
    #
    min_mapq_uniq = args.mpqu_inp
    sys.stdout.write('minimum mapping quality for uniq test: '+str(min_mapq_uniq)+'\n\n')
    #
    list_file = args.list_inp
    #
    # Read initial prediction list
    with open(list_file, 'r') as inp_file:
        inp_file_lines = inp_file.readlines()
    inp_file.close()
    #
    #
    eprint("----------\n"  + "remap and find insertion points\n" + "----------\n")
    for lf_line in inp_file_lines:
        #
        if lf_line.startswith('#'):
            continue
        #
        chrom = lf_line.strip().split()[1]
        insert_guess = int(lf_line.strip().split()[2])
        #
        int_file_name = str(chrom)+'_'+str(insert_guess)
        samfile = pysam.AlignmentFile(bam_full, 'rb')
        #subprocess.run(['mkdir' , '-p' , int_file_name])
        #os.chdir(int_file_name)
        tmp_censor_dir = preprocess_dir_realpath + '/' + 'censor_results/' + int_file_name
        subprocess.run(['mkdir' , '-p' , tmp_censor_dir])
        os.chdir(tmp_censor_dir)
        output_file = open(str(chrom)+'_'+str(lf_line.split()[2])+'.out', 'w')
        #    
        #samfile = pysam.AlignmentFile("../" + bam_full, 'rb')
        insert_guess_start_range = insert_guess-(insert_size+insert_range)
        insert_guess_end_range = insert_guess+(insert_size+insert_range)
        if insert_guess_start_range < 1:
            insert_guess_start_range = 1
        iterator_reads = samfile.fetch(chrom, insert_guess_start_range, insert_guess_end_range)
        #
        # Convert iterator to list for multiple usages
        iterator_reads_list = list( iterator_reads )
        # Close samfile. Don't close it before copying iterator to list, will cause segfault
        samfile.close()
        #
        # discover clipped reads 
        array_p, array_n = find_clipped_ends( iterator_reads_list, insert_guess, args ) 
        #seperate funct for adding more functionality later
        #
        output_file.write('(+)end clipped at: ' + ', '.join(list(map(str, array_p))) + '\n')
        output_file.write('(-)end clipped at: ' + ', '.join(list(map(str, array_n))) + '\n')
        #
        clipped_ends_p, clipped_ends_n = [], []
        #
        if len(array_p) == 0:
            output_file.write("No clipped reads on (+) end near provided guess. \
                        You may wanna change anchor_length/clipped_length parameters!\n")
        else:
            clipped_ends_p = break_points(array_p, min_reads_ends, cliped_end_range)
        #
        if len(array_n) == 0:
            output_file.write("No clipped reads in (-) end near provided guess. \
                        You may wanna change anchor_length/clipped_length parameters!\n")
        else:
            clipped_ends_n = break_points(array_n, min_reads_ends, cliped_end_range)
        #----------------------------------------------------------------------------------------------
        calc_tsd_flag = 'y'
        #
        if len(clipped_ends_p) == 0:
            if len(clipped_ends_n) > 0: # <- chnage this order, bring array_p first
                clipped_ends_p = clipped_ends_n.copy()
            elif len(array_p) > 0:
                clipped_ends_p = array_p.copy()
        #
        if len(clipped_ends_n) == 0:
            if len(clipped_ends_p) > 0: # <- change this order, bring array_n first
                clipped_ends_n = clipped_ends_p.copy()
            elif len(array_n) > 0:
                clipped_ends_n = array_n.copy()
        #
        if len(clipped_ends_p) == 0:
            if len(clipped_ends_n) > 0:
                clipped_ends_p = clipped_ends_n.copy()
            else:
                 clipped_ends_p.append(insert_guess)
            output_file.write("No clipped cluster on (+) end. You may wanna change \
                        anchor_length/clipped_length parameters!\n")
            calc_tsd_flag = 'n'
        #
        if len(clipped_ends_n) == 0:
            if len(clipped_ends_p) > 0:
                clipped_ends_n = clipped_ends_p.copy()
            else:
                clipped_ends_n.append(insert_guess)
            output_file.write("No clipped cluster on (-) end. You may wanna change \
                        anchor_length/clipped_length parameters!\n")
            calc_tsd_flag = 'n'
        #
        # find insertion point *NEAREST* to provided guess
        insert_point = 0
        gap_ends = 0 # make first gap_end = max_tsd + 2*range of clipped read search (2*5)
        loop_counter = 0
        insert_point_set_flag = 0
        break_flag = 0
        while (insert_point == 0):
            if break_flag:
                break
            if gap_ends < 2*insert_range:
                gap_ends += clipped_search_interval
            pend = 0
            nstrt = 0
            gap = insert_range
            for i in clipped_ends_p:
                if break_flag:
                    break
                for j in clipped_ends_n:
                    loop_counter += 1
                    if loop_counter > 1000:
                        break_flag = 1
                        break
                    if abs(j-i) < gap_ends:
                        if abs(insert_guess-(abs(j+i)/2)) <= gap:
                            gap = abs(insert_guess-(abs(j-i)/2))
                            pend = i
                            nstrt = j
            if pend != 0:
                insert_point = round((pend+nstrt)/2)
                insert_point_out = pend
                insert_point_set_flag = 1
        if insert_point_set_flag == 0:
            insert_point_out = insert_guess
            insert_point = insert_guess
            warning = "WARNING: timed out estimating insertion point for " + str(chrom) + "_" + str(insert_guess)
            eprint(warning)
            output_file.write(warning)
            
        output_file.write('TE insertion point closest to provided guess is: ' + \
                    str(insert_point_out) + ', which is (+)cluster ' + str(pend)+'\n')
        #----------------------------------------------------------------------------------------
        if calc_tsd_flag == 'y':
            tsd_length = pend - nstrt
            if tsd_length <= 0:
                tsd_length = 0
        if calc_tsd_flag == 'n':
            tsd_length = 'NA'

        output_file.write('\n'+'TSD length (+-5) is: '+str(tsd_length)+'\n')

        het_tsd_length = abs(pend - nstrt)

        #-----------------------------------------------------------------------------------------
        # Writes clipped part of reads to a file
        # 
        output_file.write("\nDetails of clipped reads at insertion point:\n-----------------------------------------------\n")
        #
        insert_range_new = round(gap_ends/2) + 5

        file_reads_p = open(str(int_file_name)+'_reads_p.fa','w')
        file_reads_n = open(str(int_file_name)+'_reads_n.fa','w')

        cnt_rd_p, cnt_rd_n = 0, 0

        for read in iterator_reads_list:
            if read.reference_start in range(insert_point-insert_range, insert_point+insert_range) and \
                read.reference_end in range(insert_point-insert_size, insert_point+insert_size) and \
                (read.cigartuples != None) and (read.is_supplementary != True) and ( read.has_tag('XA') or \
                read.has_tag('SA') or read.mapping_quality >= min_mapq ) and (read.mapping_quality >= min_mapq_uniq):

                write_flag = check_uniq_mapping( read, args )

                # Clipped at p-end of TE
                if read.cigartuples[-1][0] == 4 and read.cigartuples[-1][1] > clipped_length and write_flag == 'y':
                    if read.reference_end-read.reference_start > anchor_length and read.reference_end in \
                        range(insert_point-insert_range_new, insert_point+insert_range_new):
                        #
                        cnt_rd_p += 1
                        file_reads_p.write('>' + str(cnt_rd_p) + ' ' + read.query_name  + '\n' \
                            + read.query_sequence[read.infer_query_length()-read.cigartuples[-1][1]-1:read.infer_query_length()] + '\n')
        #                if printlev == 2:
                        output_file.write('ref start = ' + str(read.reference_start) + '; ref end = ' \
                            + str(read.reference_end) + '; #mapped bp = ' + str(read.reference_end-read.reference_start) \
                            + '; #clipped bp = ' + str(read.cigartuples[-1][1]) + '; cigar = ' + str(read.cigarstring) \
                            + '; id = ' + str(read.query_name)+'\n')
                        #
                        output_file.write(read.query_sequence+'\n')
                        #
                        output_file.write(read.query_sequence[0:read.reference_end-read.reference_start-1] + ' * ' \
                            + read.query_sequence[read.infer_query_length()-read.cigartuples[-1][1]-1:read.infer_query_length()] + ' (+)\n')
                        # -1 in above line is for making it consistent with 0 based index

                # Clipped at n-end of TE
                if read.cigartuples[0][0] == 4 and read.cigartuples[0][1] > clipped_length and write_flag == 'y':
                    if read.reference_end-read.reference_start > anchor_length and read.reference_start in \
                        range(insert_point-insert_range_new, insert_point+insert_range_new):
                        #
                        cnt_rd_n += 1
                        file_reads_n.write('>' + str(cnt_rd_n) + ' ' + read.query_name  + '\n' \
                            + read.query_sequence[0:read.cigartuples[0][1]-1] + '\n')
        #                if printlev == 2:
                        output_file.write('ref start = ' + str(read.reference_start) + '; ref end = ' + str(read.reference_end) \
                            + '; #mapped bp = ' + str(read.reference_end-read.reference_start) + '; #clipped bp = ' \
                            + str(read.cigartuples[0][1]) + '; cigar = ' + str(read.cigarstring) + '; id = ' + str(read.query_name)+'\n')
                        #
                        output_file.write(read.query_sequence+'\n')
                        #
                        output_file.write(read.query_sequence[0:read.cigartuples[0][1]-1] + ' * ' \
                            + read.query_sequence[read.cigartuples[0][1]-1:read.cigartuples[0][1]+read.reference_end-read.reference_start] + ' (-)\n')

        file_reads_p.close()
        file_reads_n.close()

        output_file.write("\nTotal number of clipped reads on (+) side: "+str(cnt_rd_p))
        output_file.write("\nTotal number of clipped reads on (-) side: "+str(cnt_rd_n)+'\n')

        #--------------------------------------------------------------------------------------------------------------------------
        # Get TE type (L1, SINE, etc)
        #--------------------------------------------------------------------------------------------------------------------------
        #print(cnt_rd_p)
        #print(cnt_rd_n)

        if cnt_rd_p > 0: # Number of clipped reads at p-end of TE >1
        #    proc = Popen(['censor.ncbi', int_file_name+'_reads_p.fa' , '-lib' , te_type_file], shell=True, stdout=PIPE, stderr=PIPE)
        #    out, err = proc.communicate()
        #    print (err, proc.returncode)
            subprocess.run(['censor.ncbi', int_file_name+'_reads_p.fa' , '-lib' , te_type_file])
        if cnt_rd_n > 0: # Number of clipped reads at n-end of TE >1
#            proc = Popen(['censor.ncbi', int_file_name+'_reads_n.fa' , '-lib' , te_type_file], shell=True, stdout=PIPE, stderr=PIPE)
#            out, err = proc.communicate()
            subprocess.run(['censor.ncbi', int_file_name+'_reads_n.fa' , '-lib' , te_type_file])
        #---------------------------------------------------------------------------------------------------------------------------
        output_file.write("\nMapping information of clipped part of read\n------------------------------------------------\n")
        #
        clipped_read_p_flag = 'n'
        clipped_read_n_flag = 'n'
        #
        clipped_read_p_flag, type_clipped_p, output_write_lines = \
                read_type_info('clipped', 'p', int_file_name+'_reads_p.fa.map', args)
        #
        output_file.write(''.join(output_write_lines))
        del output_write_lines[:]
        #
        clipped_read_n_flag, type_clipped_n, output_write_lines = \
                read_type_info('clipped', 'n', int_file_name+'_reads_n.fa.map', args)
        #
        output_file.write(''.join(output_write_lines))
        del output_write_lines[:]
        #
        len_calc_flag = 'n'
        #fofn_ref_filename = dir_path + "/" + args.fofn_ref
        len_calc_flag, te_class_file, output_write_lines = te_type_setup( int_file_name+'_reads_p.fa.map', \
            int_file_name+'_reads_n.fa.map', type_clipped_p, type_clipped_n, cnt_rd_p, cnt_rd_n, fofn_ref_realpath)
        output_file.write(''.join(output_write_lines))
        del output_write_lines[:]
        #
        # Remove intermediate files
        subprocess.run(['rm' , '-rf' , '*.fa.*' , 'censor.ncbi.*' , 'error.log' , 'formatdb.log'])
        #call(['rm' , '-rf' , '*.fa.map' , 'censor.ncbi.*' , '*.log'])
        #os.chdir('..')
        #continue
        #---------------------------------------------------------------------------------------------------------------------------
        if len_calc_flag == 'y':
            subprocess.run(['censor.ncbi', int_file_name+'_reads_p.fa' , '-lib' , te_class_file])
            subprocess.run(['censor.ncbi', int_file_name+'_reads_n.fa' , '-lib' , te_class_file])
            if check_file("%s_reads_p.fa.map" % int_file_name) != True or \
                        check_file("%s_reads_n.fa.map" % int_file_name) !=True:
                len_calc_flag = 'n'

        if len_calc_flag == 'y':
            output_file.write("\nMapping information of clipped part of read to TE class.\
                        \n---------------------------------------------------------\n")        
            class_out_p, output_write_lines = read_class_info( 'clipped', 'p', int_file_name+'_reads_p.fa.map', args )
            output_file.write(''.join(output_write_lines))
            del output_write_lines[:]

            class_out_n, output_write_lines = read_class_info( 'clipped', 'n', int_file_name+'_reads_n.fa.map', args )
            output_file.write(''.join(output_write_lines))
            del output_write_lines[:]

            len_arr=[]
            clipped_type_both = None    

            if class_out_p[0][0] == class_out_n[0][0]:
                clipped_type_both = 'Y'
                output_file.write( str( class_out_p[0][0] )+' '+ str( calc_length( class_out_p, class_out_n ) )+'\n' )
                len_arr.append( [class_out_p[0][0], calc_length( class_out_p, class_out_n )] )
            
            if class_out_p[0][0] != class_out_n[0][0]:
                clipped_type_both = 'N'

                #subprocess.run(['rm', '-rf' , 'p.clipped.fa*' , 'n.clipped.fa*'])

                tmp = pysam.faidx( te_class_file, class_out_n[0][0] )
                f_tmp = open( class_out_n[0][0]+'.fa', 'w' )
                f_tmp.write(tmp)
                f_tmp.close()

                tmp = pysam.faidx( te_class_file, class_out_p[0][0] )
                f_tmp = open( class_out_p[0][0]+'.fa', 'w' )
                f_tmp.write(tmp)
                f_tmp.close()

                subprocess.run(['cp' , int_file_name+'_reads_n.fa' , 'n.clipped.fa'])
                subprocess.run(['censor.ncbi' , 'n.clipped.fa' , '-lib' , class_out_p[0][0]+'.fa'])
                #
                if check_file('n.clipped.fa.map') == True:
                    output_file.write(str(class_out_p[0][0]) + ' ' + str( get_class('n.clipped.fa.map', args)[0][0] ) + ' ' + \
                                    str( calc_length( class_out_p, get_class( 'n.clipped.fa.map', args ) ) ) + '\n') 
                    len_arr.append([class_out_p[0][0], calc_length(class_out_p, get_class('n.clipped.fa.map', args))])
            
                subprocess.run(['cp' , int_file_name+'_reads_p.fa' , 'p.clipped.fa'])
                subprocess.run(['censor.ncbi' , 'p.clipped.fa' , '-lib' , class_out_n[0][0]+'.fa'])
                if check_file('p.clipped.fa.map') == True:
                    output_file.write(str(class_out_n[0][0]) + ' ' + str(get_class('p.clipped.fa.map', args)[0][0]) + ' ' + \
                                    str(calc_length(class_out_n, get_class('p.clipped.fa.map', args))) + '\n') 
                    len_arr.append([class_out_n[0][0], calc_length(class_out_n, get_class('p.clipped.fa.map', args))])


            print_tup(len_arr, output_file, '-' , ' ')
            #subprocess.run(['rm', '-rf', '*.fa.*', 'censor.ncbi.*', 'error.log', 'formatdb.log'])
        
        if len_calc_flag == 'n':
            output_file.write("\nMapping information of clipped part of read to TE class.\
                        \n---------------------------------------------------------\n")
            output_file.write('Not attempting TE length calculation...\n')

    #---------------------------------------------------------------------------------------------------------------------------
    # Het vs Hom
    #---------------------------------------------------------------------------------------------------------------------------
        #eprint("----------\n"  + "determining if het or not\n" + "----------\n")
        cnt_het = 0
        for read in iterator_reads_list:
            if read.cigartuples != None:  # prevents accidental kill. remember, read.reference_start is 0 based.
                if read.reference_start < int(insert_point_out-het_tsd_length-5) \
                    and read.reference_end > int(insert_point_out+het_tsd_length+5):
                    cnt_het += 1

        output_file.write("\n\n------------------------------------------------\n")

        if cnt_het >= min_reads_het:
            output_file.write("This is a hetrozygous insertion")
            het_hom = 'Hetrozygous'
        else:
            output_file.write("This is a homozygous insertion")
            het_hom = 'Homozygous'

        output_file.write("\n------------------------------------------------\n\n")
    
    #----------------------------------------------------------------------------------------------------------------------------
    # Poly A/T
    #----------------------------------------------------------------------------------------------------------------------------
       # eprint("----------\n"  + "determining polyA\n" + "----------\n")
    # 1. p-end    
        pat_test_out_p = analysis_pat_check( int_file_name+'_reads_p.fa', int_file_name+'_reads_p.fa.map' )
        # num_pa_p = sum(x.count('A') for x in pat_test_out_p)        
        # num_pt_p = sum(x.count('T') for x in pat_test_out_p)
        num_pat_p = len(pat_test_out_p)    
        open(int_file_name+'_pat_p', 'w').write( '\n'.join('\t'.join(str(dat) for dat in itm) for itm in pat_test_out_p) )
        #
    # 2. n-end
        pat_test_out_n = analysis_pat_check( int_file_name+'_reads_n.fa', int_file_name+'_reads_n.fa.map' )        
        num_pat_n = len(pat_test_out_n)    
        open(int_file_name+'_pat_n', 'w').write( '\n'.join('\t'.join(str(dat) for dat in itm) for itm in pat_test_out_n) )

    #----------------------------------------------------------------------------------------------------------------------------

        output_file.write( "\n-----Final results ( Clipped )-----\n" )
        output_file.write( '@clipped\t'+ chrom +'\t'+ str(insert_guess) +'\t'+ str(insert_point_out)+'\t' )
        output_file.write( het_hom +'\t'+ str(cnt_het) +'\t'+ str(tsd_length) +'\t' )
        #
        if clipped_read_p_flag == 'y':
            if cnt_rd_p > 0:
                type_c_p_out = str(float("{0:.2f}".format((type_clipped_p[0][2]/cnt_rd_p)*100)))
            else:
                type_c_p_out = 'NA'
            #output_file.write(str(type_clipped_p[0][0])+'\t'+ str(type_clipped_p[0][1])+'\t'+ \
            #    str(type_clipped_p[0][2]) +'\t'+ str(float("{0:.2f}".format((type_clipped_p[0][2]/cnt_rd_p)*100)))+'\t')
            output_file.write(str(type_clipped_p[0][0])+'\t'+ str(type_clipped_p[0][1])+'\t'+  str(type_clipped_p[0][2]) +'\t'+ type_c_p_out +'\t')
        if clipped_read_p_flag == 'n':
            output_file.write("NA\tNA\t0\t0\t")
        #
        if clipped_read_n_flag == 'y':
            if cnt_rd_n > 0:
                type_c_n_out = str(float("{0:.2f}".format((type_clipped_n[0][2]/cnt_rd_n)*100)))
            else:
                type_c_n_out = 'NA'
            #output_file.write(str(type_clipped_n[0][0])+'\t'+ str(type_clipped_n[0][1])+'\t'+ \
            #    str(type_clipped_n[0][2]) +'\t'+ str(float("{0:.2f}".format((type_clipped_n[0][2]/cnt_rd_n)*100)))+'\t')
            output_file.write(str(type_clipped_n[0][0])+'\t'+ str(type_clipped_n[0][1])+'\t'+ str(type_clipped_n[0][2]) +'\t' + type_c_n_out + '\t')
        if clipped_read_n_flag == 'n':
            output_file.write("NA\tNA\t0\t0\t")
        #
        if len_calc_flag == 'y':
            if len(len_arr) == 2:
                if class_out_p[0][3] < class_out_n[0][3]:
                    output_file.write(class_out_n[0][0]+'\t'+ clipped_type_both+'\t'+ str(len_arr[1][1])+'\t')
                else:
                    output_file.write(class_out_p[0][0]+'\t'+ clipped_type_both+'\t'+ str(len_arr[0][1])+'\t')
            if len(len_arr) == 1:
                output_file.write(str(len_arr[0][0])+'\t'+ clipped_type_both+'\t'+ str(len_arr[0][1])+'\t')
            if len(len_arr) == 0:
                output_file.write('NA\tNA\t0\t')
            del len_arr
        if len_calc_flag == 'n':
            output_file.write('NA\tNA\t0\t')
        #
        output_file.write(str(gap_ends)+'\t')
        #
        output_file.write(str(num_pat_p) + '\t' + str(num_pat_n) + '\t' )
        #
        del array_p[:], array_n[:], clipped_ends_p[:], clipped_ends_n[:]
        del iterator_reads_list[:]

        output_file.close()
        #
        #os.chdir('..')
        os.chdir(dir_path)
    #*****************************************************************************************************************************
    #*****************************************************************************************************************************
    #Discordant mate pair analysis
    #
    eprint("----------\n"  + "discordant mate pair analysis\n" + "----------\n")
    samfile_idx = pysam.AlignmentFile(discord_bam, "rb")
    id_index = pysam.IndexedReads(samfile_idx)
    id_index.build()

    for lf_line in inp_file_lines:
        #
        if lf_line.startswith('#'):
            continue
        #    
        samfile = pysam.AlignmentFile(bam_full, "rb")
        chrom = lf_line.strip().split()[1]
        insert_guess = int(lf_line.strip().split()[2])
        #
        int_file_name = str(chrom)+'_'+str(insert_guess)
        #os.chdir(int_file_name)
        tmp_censor_dir = preprocess_dir_realpath + '/' + 'censor_results/' + int_file_name
        os.chdir(tmp_censor_dir)
        #
        del insert_point
        with open( int_file_name+'.out', 'r') as tmpfile:
            for line in tmpfile:
                if line.startswith('@clipped'):    
                    insert_point = int(line.split()[3])
        tmpfile.close()
        #
        output_file = open(str(chrom)+'_'+str(insert_guess)+'.out', 'a+')
        #
        #samfile = pysam.AlignmentFile("../" + bam_full, "rb")
        #iterator_reads = samfile.fetch(chrom, insert_point-insert_size, insert_point+insert_size)
        insert_point_start = insert_point-insert_size 
        insert_point_end = insert_point+insert_size 
        if insert_point_start < 1:
            insert_point_start = 1
        iterator_reads = samfile.fetch(chrom, insert_point_start, insert_point_end)
        iterator_reads_list = list( iterator_reads )
        samfile.close()
        #
        #----------------------------------------------------------------------------------------------------------------------------
        output_file.write("\n\n\nDiscordant mate pair analysis\n------------------------------------------------\n")
        #
        file_discord_p = open( int_file_name+'_discord_p', 'w')
        file_discord_n = open( int_file_name+'_discord_n', 'w')
        #
        for read in iterator_reads_list:
            if (read.is_paired == True) and (read.is_proper_pair != True) and (read.cigarstring != None) \
                and (read.is_supplementary != True) and ( read.has_tag('XA') or read.has_tag('SA') \
                or read.mapping_quality >= min_mapq ) and (read.mapping_quality >= min_mapq_uniq):    
                #
                write_flag = check_uniq_mapping( read, args )
                #    
                #if read.get_reference_positions()[0] in range(insert_point-insert_size, insert_point) and write_flag == 'y':
                if read.get_reference_positions()[0] in range(insert_point_start, insert_point) and write_flag == 'y':
                    file_discord_p.write(read.query_name + '\t' + str(read.flag) + '\t' + str(read.reference_start) + '\t' + str(read.reference_end) + '\n')        
                #
                #if read.get_reference_positions()[-1] in range(insert_point, insert_point+insert_size+1) and write_flag == 'y':
                if read.get_reference_positions()[-1] in range(insert_point, insert_point_end+1) and write_flag == 'y':
                    file_discord_n.write(read.query_name + '\t' + str(read.flag) + '\t' + str(read.reference_start) + '\t' + str(read.reference_end) + '\n')
        #    
        file_discord_p.close()
        file_discord_n.close()
        #----------------------------------------------------------------------------------------------------------------------------
        #
        with open(str(int_file_name)+'_discord_p', 'r') as f1, open(str(int_file_name)+'_discord_n', 'r') as f2:
            p_lines = f1.readlines()
            n_lines = f2.readlines()
        f1.close()
        f2.close()    
        #
        cfl_lines = []
        for i in p_lines:
            for j in n_lines:
        #        if ( i.strip().split()[0]==j.strip().split()[0] ) and ( i.strip().split()[1]==j.strip().split()[1] ) :
                if i.strip() == j.strip():
                    cfl_lines.append(i)
        #
        # Split common lines
        cfl_split_p = []
        cfl_split_n = []
        for cline in cfl_lines:
            if ( insert_point - int(cline.strip().split()[2]) ) >= ( int(cline.strip().split()[3]) - insert_point ):
                cfl_split_p.append(cline)
            else:
                cfl_split_n.append(cline)    
        #    
        with open(int_file_name+'_discord_pu', 'w') as tmp_pu:
            #
            for line in p_lines:
                if line not in cfl_lines: # not a duplicate
                    tmp_pu.write(line)
            #
            tmp_pu.write(''.join(cfl_split_p))  # <-- common read split 
            #
        tmp_pu.close()    
        #
        with open(int_file_name+'_discord_nu', 'w') as tmp_nu:
            #
            for line in n_lines:
                if line not in cfl_lines: # not a duplicate
                    tmp_nu.write(line)
            #
            tmp_nu.write(''.join(cfl_split_n))  # <-- common read split
            #
        tmp_nu.close()
        #
        # write common line to a file
        open(str(int_file_name)+'_discord_c', 'w').write(''.join(cfl_lines))
        #
        output_file.write( 'Number of discordant reads at (+) end: ' + str( num_of_lines( int_file_name+'_discord_pu' ) ) + '\n' )
        output_file.write( 'Number of discordant reads at (-) end: ' + str( num_of_lines( int_file_name+'_discord_nu' ) ) + '\n' )
        # Removing duplicate lines done
        #----------------------------------------------------------------------------------------------------------------------------
        # Call mate of a pair

        cnt_dm_p = 0
        if check_file( int_file_name+'_discord_pu' ) == True:
            file_discord_mate_p = open(str(int_file_name)+'_discord_mate_p.fa','w')
            for text_p in open(str(int_file_name)+'_discord_pu','r'):
                iterator = id_index.find(text_p.strip().split()[0])
                for read in iterator:
                    if read.query_name == text_p.split()[0] and (int(read.flag) & 0x40) != (int(text_p.split()[1]) & 0x40): # 1 is first and 2nd is second
                        cnt_dm_p += 1
                        file_discord_mate_p.write('>' + str(cnt_dm_p) + ' ' + read.query_name + '\n' + read.query_sequence + '\n')

            file_discord_mate_p.close()

        cnt_dm_n = 0
        if check_file( int_file_name+'_discord_nu' ) == True:
            file_discord_mate_n = open(str(int_file_name)+'_discord_mate_n.fa','w')
            for text_n in open(str(int_file_name)+'_discord_nu','r'):
                iterator = id_index.find(text_n.strip().split()[0])
                for read in iterator:
                    if read.query_name == text_n.split()[0] and (int(read.flag) & 0x40) != (int(text_n.split()[1]) & 0x40):
                        cnt_dm_n += 1
                        file_discord_mate_n.write('>' + str(cnt_dm_n) + ' ' + read.query_name + '\n' + read.query_sequence + '\n')

            file_discord_mate_n.close()

        output_file.write("Number of discordant mates for (+) end: " + str(cnt_dm_p)+'\n')
        output_file.write("Number of discordant mates for (-) end: " + str(cnt_dm_n)+'\n')
        #-------------------------------------------------------------------------------------------------------------
        if args.flt_inp:
            #
            if check_file( int_file_name+'_discord_mate_p.fa' ):
                flt_discord( int_file_name+'_discord_mate_p.fa' )
                #
            if check_file( int_file_name+'_discord_mate_n.fa' ):
                flt_discord( int_file_name+'_discord_mate_n.fa' )
                #
        #---------------------------------------------------------------------------------------------------------------
        if cnt_dm_p > 0:
            subprocess.run(['censor.ncbi' , int_file_name+'_discord_mate_p.fa' , '-lib' , te_type_file])
        if cnt_dm_n > 0:
            subprocess.run(['censor.ncbi' , int_file_name+'_discord_mate_n.fa' , '-lib' , te_type_file])
        #--------------------------------------------------------------------------------------------------------------
        output_file.write("\nMapping information of mates of discordant reads\n------------------------------------------------\n")
        #
        discord_read_p_flag = 'n'
        discord_read_n_flag = 'n'
        #
        discord_read_p_flag, type_discord_p, output_write_lines = read_type_info('discord', 'p', int_file_name+'_discord_mate_p.fa.map', args )
        output_file.write(''.join(output_write_lines))
        del output_write_lines[:]
        #
        discord_read_n_flag, type_discord_n, output_write_lines = read_type_info('discord', 'n', int_file_name+'_discord_mate_n.fa.map', args )
        output_file.write(''.join(output_write_lines))
        del output_write_lines[:]
        #
        discord_class_flag = 'n'
        #fofn_ref_filename = dir_path + "/" + args.fofn_ref
        discord_class_flag, te_class_file, output_write_lines = te_type_setup( int_file_name+'_discord_mate_p.fa.map', \
                int_file_name+'_discord_mate_n.fa.map', type_discord_p, type_discord_n, cnt_dm_p, cnt_dm_n, fofn_ref_realpath )
        output_file.write(''.join(output_write_lines))
        del output_write_lines[:]
        #
        subprocess.run(['rm', '-rf', '*.fa.*', 'censor.ncbi.*', 'error.log', 'formatdb.log'])
        #call(['rm', '-rf', '*.fa.map', 'censor.ncbi.*', '*.log' ])    
        #--------------------------------------------------------------------------------------------------------------
        #Type estimation for discordant mates
        if discord_class_flag == 'y':
            subprocess.run(['censor.ncbi' , int_file_name+'_discord_mate_p.fa' , '-lib' , te_class_file])
            subprocess.run(['censor.ncbi' , int_file_name+'_discord_mate_n.fa' , '-lib' , te_class_file])
            #
            if check_file( int_file_name+'_discord_mate_p.fa.map' ) !=True or check_file( int_file_name+'_discord_mate_n.fa.map' ) !=True:
                discord_class_flag = 'n'
            #
        if discord_class_flag == 'y':
            output_file.write("\nMapping information of discordant mates to TE class.\n---------------------------------------------------------\n")

            class_out_discord_p, output_write_lines = read_class_info('discord', 'p', int_file_name+'_discord_mate_p.fa.map', args )
            output_file.write(''.join(output_write_lines))
            del output_write_lines[:]
        
            class_out_discord_n, output_write_lines = read_class_info('discord', 'n', int_file_name+'_discord_mate_n.fa.map', args )
            output_file.write(''.join(output_write_lines))
            del output_write_lines[:]

            discord_class_both = None
            if class_out_discord_p[0][0] == class_out_discord_n[0][0]:
                discord_class_both = 'Y'
            else:
                discord_class_both = 'N'

        subprocess.run(['rm', '-rf', '*.fa.*', 'censor.ncbi.*', 'error.log', 'formatdb.log'])
        
        output_file.write("\n-----Final results (Discordant reads)-----\n")
        output_file.write( '@discord\t' + chrom + '\t' + str(insert_guess) + '\t' + str(insert_point) + '\t' )

        if discord_read_p_flag == 'y':
            output_file.write(str(type_discord_p[0][0])+'\t'+str(type_discord_p[0][1])+'\t'+str(type_discord_p[0][2]) +'\t' \
                                +str(float("{0:.2f}".format((type_discord_p[0][2]/cnt_dm_p)*100)))+'\t' )
        if discord_read_p_flag == 'n':
            output_file.write("NA\tNA\t0\t0\t")

        if discord_read_n_flag == 'y':
            output_file.write(str(type_discord_n[0][0])+'\t'+str(type_discord_n[0][1])+'\t'+str(type_discord_n[0][2]) +'\t' \
                                + str(float("{0:.2f}".format((type_discord_n[0][2]/cnt_dm_n)*100)))+'\t' )
        if discord_read_n_flag == 'n':
            output_file.write("NA\tNA\t0\t0\t")

        if discord_class_flag == 'y':
            if class_out_discord_p[0][3] < class_out_discord_n[0][3]:
                output_file.write(class_out_discord_n[0][0]+'\t'+ discord_class_both+'\t' )
            else:
                output_file.write(class_out_discord_p[0][0]+'\t'+discord_class_both+'\t')
        if discord_class_flag == 'n':
            output_file.write("NA\tNA\t")
        #
        output_file.write('\n')
        #os.chdir('..')
        os.chdir(dir_path)

        output_file.close()
        del iterator_reads_list[:]

    #sys.exit()
    samfile_idx.close()
    #--------------------------------------------------------------------------------------------------------------
    # Write combined output file.
    eprint("----------\n"  + "writing output files\n" + "----------\n")
    #final_output = open('final_results.tsv', 'w')
    with open(args.output_file, 'w') as final_output:
        final_output.write("Type\t# Chromosome\tInitial_Guess\tActual_insertion_point\tHetrozygous/Homozygous\t"
            + "#reads_forHet\tTSD_length\t(+)clipped_type\t(+)clipped_type_quality\t"
            + "#reads_supporting_(+)clipped_type\t%reads_supporting_(+)clipped_type\t(-)clipped_type\t"
            + "(-)clipped_type_quality\t#reads_supporting_(-)clipped_type\t%reads_supporting_(-)clipped_type\t"
            + "Estimated_TE_class\tBoth_clipped_end_support\tEstimated_TE_Length(bp)\tgap_bw_ends\t"
            + "num_pat_p\tnum_pat_n\t"
            + "(+)discord_mate_type\t(+)discord_mate_type_quality\t#reads_supporting_(+)discord_mate_typet\t"
            + "%reads_supporting_(+)discord_mate_type\t(-)discord_mate_type\t(-)discord_mate_type_quality\t"
            + "#reads_supporting_(-)discord_mate_type\t%reads_supporting_(-)discord_mate_type\t"
            + "Estimated_TE_class-discord\tBoth_end_discord_mate_support\tspecial_comment\n")

        for lf_line in inp_file_lines:
            #
            if lf_line.startswith('#'):
                continue
            #
            chrom = lf_line.split()[1]
            insert_guess = int(lf_line.split()[2])
            int_file_name = str(chrom)+'_'+str(insert_guess)
            tmp_censor_dir = preprocess_dir_realpath + '/' + 'censor_results/' + int_file_name
            #output_file = open('./'+int_file_name+'/'+str(chrom)+'_'+str(insert_guess)+'.out', 'r')
            output_file = open(tmp_censor_dir+'/'+str(chrom)+'_'+str(insert_guess)+'.out', 'r')
            for line in output_file:
                if line.startswith('@clipped'):
                    final_output.write(lf_line.split()[0]+'\t'+line.split('\t', 1)[1].rstrip('\n'))
                if line.startswith('@discord'):
                    final_output.write(line.split('\t', 4)[4])
            #
            output_file.close()
        final_output.close()    

def main():
    FUNCTION_MAP = {
            'preprocess' : exec_preprocess, 
            'discover' : exec_discover,
            'nadiscover' : exec_nadiscover, 
            'analyze' : exec_analyze,
            'cluster2d' : exec_cluster2d,
            'filter' : exec_filter,
            'filter_p' : exec_filter_p
            }

    parser = argparse.ArgumentParser()
    #subparsers = parser.add_subparsers(dest='command', required=True)
    subparsers = parser.add_subparsers(dest='command')
    subparsers.required = True

    sp_preprocess = subparsers.add_parser('preprocess', help="preprocess argument")
    sp_preprocess.add_argument('-i', '--input_bam', action='store', dest='bam_inp', required=True,
        help='input Bam(.bam) file of aligned reads')
    sp_preprocess.add_argument('-r', '--ref', action='store', dest='fofn_ref', required=True,
        help='File with reference sequence paths, see README.md for more info')
    sp_preprocess.add_argument('-p', '--preprocess_dir', action='store',
        dest='preprocess_dir', default='preprocessed_files',
        help='directory to store preprocessing output files (default: preprocessed_files)')
    sp_preprocess.add_argument('--cll', action='store', dest='cll_inp', type=int, default=25, help='Minimum clipped length(bp)')

    sp_discover = subparsers.add_parser('discover', help="discover argument")
    sp_discover.add_argument('-i', '--input_bam', action='store', dest='bam_inp', required=True, 
        help='input Bam(.bam) file of aligned reads')
    #sp_discover.add_argument('-ref', action='store', dest='fofn_ref', required=True, help='FoFn for reference sequence')
    sp_discover.add_argument('-r', '--ref', action='store', dest='fofn_ref', required=True,
        help='File with reference sequence paths, see README.md for more info')
    sp_discover.add_argument('--isz', action='store', dest='isz_inp', type=int, default=340, help='insert Size estimate')
    sp_discover.add_argument('--rdl', action='store', dest='rdl_inp', type=int, default=150, help='Average read length')
    sp_discover.add_argument('--drd', action='store', dest='drd_inp', type=int, default=10, help='discord read clust denst')
    sp_discover.add_argument('--cct', action='store', dest='cct_inp', type=int, default=200, help='Coverage cutoff input')
    sp_discover.add_argument('--cll', action='store', dest='cll_inp', type=int, default=25, help='Minimum clipped length(bp)')
    sp_discover.add_argument('--mpq', action='store', dest='mpq_inp', type=int, default=30, help='Minimum mapping quality')
    sp_discover.add_argument('--mpqu', action='store', dest='mpqu_inp', type=int, default=1, help='Minimum mapping quality')
    sp_discover.add_argument('-p', '--preprocess_dir', action='store',
        dest='preprocess_dir', default='preprocessed_files',
        help='directory used to store preprocessing output files (default: preprocessed_files)')

    sp_analyze = subparsers.add_parser('analyze', help="analyze argument")
    sp_analyze.add_argument('-i','--input_bam', action='store', dest='bam_inp', required=True, 
        help='input Bam(.bam) file of aligned reads')
    #sp_analyze.add_argument('-ref', action='store', dest='fofn_ref', required=True, help='FoFn for reference sequence')
    sp_analyze.add_argument('-r', '--ref', action='store', dest='fofn_ref', required=True,
        help='File with reference sequence paths, see README.md for more info')
    sp_analyze.add_argument('-p', '--preprocess_dir', action='store',
        dest='preprocess_dir', default='preprocessed_files',
        help='directory used to store preprocessing output files (default: preprocessed_files)')
    sp_analyze.add_argument('--inp', action='store', dest='list_inp', required=True, help='Input list of insertions')
    sp_analyze.add_argument('-o', '--output_file', action='store',
        dest='output_file', default='final_results.tsv',
        help='Tab-delimited output file of potential TE insertions(default: final_resutls.tsv)')
    sp_analyze.add_argument('--rdl', action='store', dest='rdl_inp', type=int, default=150, help='Average read length')
    sp_analyze.add_argument('--cll', action='store', dest='cll_inp', type=int, default=25, help='Minimum clipped length(bp)')
    sp_analyze.add_argument('--ahl', action='store', dest='ahl_inp', type=int, default=30, help='Minimum anchor length(bp)')
    sp_analyze.add_argument('--cer', action='store', dest='cer_inp', type=int, default=5, help='Range of clipped reads at a end to put in a group')
    sp_analyze.add_argument('--csi', action='store', dest='csi_inp', type=int, default=20, help='Clipped read search interval')
    sp_analyze.add_argument('--mre', action='store', dest='mre_inp', type=int, default=4, help='min read for breakpoint')
    sp_analyze.add_argument('--mrh', action='store', dest='mrh_inp', type=int, default=3, help='Minimum reads to call hetrozygous insertion')
    sp_analyze.add_argument('--isz', action='store', dest='isz_inp', type=int, default=340, help='insert Size estimate')
    sp_analyze.add_argument('--qii', action='store', dest='qii_inp', type=float, default=0.05, help='Interval for mapping quality')
    sp_analyze.add_argument('--nii', action='store', dest='nii_inp', type=int, default=6, help='Number of intervals')
    sp_analyze.add_argument('--mpq', action='store', dest='mpq_inp', type=int, default=30, help='Minimum mapping quality')
    sp_analyze.add_argument('--mpqu', action='store', dest='mpqu_inp', type=int, default=1, help='Minimum mapping quality uniq test')
    sp_analyze.add_argument('--flt', action='store_true', dest='flt_inp', default=False, help='Filter discord mate files')

    sp_nadiscover = subparsers.add_parser('nadiscover', help="nadiscover argument")
    sp_nadiscover.add_argument('-i','--input_bam', action='store', dest='bam_inp', required=True, 
        help='input Bam(.bam) file of aligned reads')
    #sp_nadiscover.add_argument('-ref', action='store', dest='fofn_ref', required=True, help='FoFn for reference sequence')
    sp_nadiscover.add_argument('-r', '--ref', action='store', dest='fofn_ref', required=True,
        help='File with reference sequence paths, see README.md for more info')
    sp_nadiscover.add_argument('-o', '--output_file', action='store',
        dest='output_file', default='initial_predictions_noalign.txt',
        help='Tab-delimited output file of initial set of TE insertions (default: initial_predictions_noalign.txt)')
    sp_nadiscover.add_argument('-p', '--preprocess_dir', action='store',
        dest='preprocess_dir', default='preprocessed_files',
        help='directory used to store preprocessing output files (default: preprocessed_files)')
    sp_nadiscover.add_argument('--cll', action='store', dest='cll_inp', type=int, default=25, help='Minimum clipped length(bp)')
    sp_nadiscover.add_argument('--isz', action='store', dest='isz_inp', type=int, default=340, help='insert Size estimate')
    sp_nadiscover.add_argument('--rdl', action='store', dest='rdl_inp', type=int, default=150, help='Average read length')
    sp_nadiscover.add_argument('--drd', action='store', dest='drd_inp', type=int, default=5, help='discord read clust denst')
    sp_nadiscover.add_argument('--cct', action='store', dest='cct_inp', type=int, default=200, help='Coverage cutoff input')
    sp_nadiscover.add_argument('--all', action='store_true', dest='flg_all', default=False, help='Only clipped or all?')
    sp_nadiscover.add_argument('--mrg', action='store_true', dest='merged', default=False, help='Merge aligned predictions?')
    sp_nadiscover.add_argument('--pat', action='store_true', dest='pat_inp', default=False, help='Poly A/T search')
    sp_nadiscover.add_argument('--nas', action='store_true', dest='nas_inp', default=False, help='Non-alignment ref bam search')
    sp_nadiscover.add_argument('--mpq', action='store', dest='mpq_inp', type=int, default=30, help='Minimum mapping quality')
    sp_nadiscover.add_argument('--pql', action='store', dest='pql_inp', type=int, default=9, help='poly A/T Length')
    sp_nadiscover.add_argument('--pmm', action='store', dest='pmm_inp', type=int, default=1, help='poly A/T mismatch')
    sp_nadiscover.add_argument('--mpqu', action='store', dest='mpqu_inp', type=int, default=1, help='Minimum mapping quality uniq test')
    sp_nadiscover.add_argument('--bed', action='store', dest='rmsk_bed', help='FoFn for existing repeat elements')

    sp_cluster2d = subparsers.add_parser('cluster2d', help="cluster2d argument")
    sp_cluster2d.add_argument('-i','--input_bam', action='store', dest='bam_inp', required=True, 
        help='input Bam(.bam) file of aligned reads')
    #sp_cluster2d.add_argument('-ref', action='store', dest='fofn_ref', required=True, help='FoFn for reference sequence')
    sp_cluster2d.add_argument('-r', '--ref', action='store', dest='fofn_ref', required=True,
        help='File with reference sequence paths, see README.md for more info')
    sp_cluster2d.add_argument('-p', '--preprocess_dir', action='store',
        dest='preprocess_dir', default='preprocessed_files',
        help='directory used to store preprocessing output files (default: preprocessed_files)')
    sp_cluster2d.add_argument('--isz', action='store', dest='isz_inp', type=int, default=340, help='insert Size estimate')
    sp_cluster2d.add_argument('--rdl', action='store', dest='rdl_inp', type=int, default=150, help='Average read length')
    sp_cluster2d.add_argument('--drd', action='store', dest='drd_inp', type=int, default=5, help='discord read clust denst')
    sp_cluster2d.add_argument('--cct', action='store', dest='cct_inp', type=int, default=200, help='Coverage cutoff input')
    sp_cluster2d.add_argument('--all', action='store_true', dest='flg_all', default=False, help='Only clipped or all?')

    sp_filter = subparsers.add_parser('filter', help="Filter argument")
    sp_filter.add_argument('--ofa', action='store', dest='ofa_inp', required=True, help='output file from analyze section')
    sp_filter.add_argument('--bed', action='store', dest='fofn_bed', required=True, help='FoFn for existing repeat elements')
    sp_filter.add_argument('-p', '--preprocess_dir', action='store',
        dest='preprocess_dir', default='preprocessed_files',
        help='directory used to store preprocessing output files (default: preprocessed_files)')
    sp_filter.add_argument('--qlm', action='store', dest='qlm_inp', type=float, default=0.85, help='Lowest limit for alignment quality')
    sp_filter.add_argument('--tcr', action='store', dest='tcr_inp', type=int, default=5, help='Minimum number of clipped reads')
    sp_filter.add_argument('--trd', action='store', dest='trd_inp', type=int, default=10, help='Minimum total [clipped+discordant] reads')
    #sp_filter.add_argument('-ref', action='store', dest='fofn_ref', help='FoFn for reference sequence')
    sp_filter.add_argument('--rp', action='store', dest='rp_inp', type=float, default=10.0, help='read percent value')
    sp_filter.add_argument('--rdl', action='store', dest='rdl_inp', type=int, default=150, help='Average read length')
    sp_filter.add_argument('--isz', action='store', dest='isz_inp', type=int, default=340, help='insert Size estimate')

    sp_filter_p = subparsers.add_parser('filter_p', help="Filter argument")
    sp_filter_p.add_argument('--ofa', action='store', dest='ofa_inp', required=True, help='output file from analyze section')
    sp_filter_p.add_argument('--bed', action='store', dest='fofn_bed', required=True, help='FoFn for existing repeat elements')
    sp_filter_p.add_argument('-p', '--preprocess_dir', action='store',
        dest='preprocess_dir', default='preprocessed_files',
        help='directory used to store preprocessing output files (default: preprocessed_files)')
    sp_filter_p.add_argument('--qlm', action='store', dest='qlm_inp', type=float, default=0.75, help='Lowest limit for alignment quality')
    sp_filter_p.add_argument('--tcr', action='store', dest='tcr_inp', type=int, default=2, help='Minimum number of clipped reads')
    sp_filter_p.add_argument('--trd', action='store', dest='trd_inp', type=int, default=5, help='Minimum total [clipped+discordant] reads')
    #sp_filter_p.add_argument('-ref', action='store', dest='fofn_ref', help='FoFn for reference sequence')
    sp_filter_p.add_argument('--rp', action='store', dest='rp_inp', type=float, default=10.0, help='read percent value')
    sp_filter_p.add_argument('--rdl', action='store', dest='rdl_inp', type=int, default=100, help='Average read length')
    sp_filter_p.add_argument('--isz', action='store', dest='isz_inp', type=int, default=369, help='insert Size estimate')

    args = parser.parse_args()
    funct = FUNCTION_MAP[args.command]
    funct(args)

if __name__ == '__main__':
    import sys
    import os
    import copy
    import pysam
    import numpy as np
    import argparse
    import re
    import subprocess
    from subprocess import Popen, PIPE
    from Bio.Seq import Seq
    from Bio.Sequencing.Applications import BwaIndexCommandline
    from Bio.Sequencing.Applications import BwaAlignCommandline
    from Bio.Sequencing.Applications import BwaSamseCommandline
    main() 

