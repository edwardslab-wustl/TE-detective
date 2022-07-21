import sys
import os
import subprocess

import numpy as np
import pysam
from Bio.Sequencing.Applications import BwaIndexCommandline
from Bio.Sequencing.Applications import BwaAlignCommandline
from Bio.Sequencing.Applications import BwaSamseCommandline

from TEdetective.io_functions import eprint
from TEdetective.general_functions import check_uniq_mapping, break_points_2d

def exec_preprocess(args):
    log_FH=open(args.log_file, 'w')
    dir_path = os.getcwd()
    log_FH.write('working directory: '+ dir_path +'\n')
    
    preprocess_dir_realpath = os.path.realpath(args.preprocess_dir)
    log_FH.write('preprocessing/intermediate file directory: '+ preprocess_dir_realpath +'\n')
    
    fofn_ref_realpath = os.path.realpath(args.fofn_ref)
    log_FH.write('fofn_ref file: '+ fofn_ref_realpath +'\n')
    
    bam_full = os.path.realpath(args.bam_inp)
    log_FH.write('Input bam file: '+str(bam_full)+'\n')
    
    bam_short_name = bam_full.split('/')[-1][:-4]
    log_FH.write('bam short name: '+ bam_short_name +'\n')
    
    clipped_length = args.cll_inp
    log_FH.write('Minimum clipped length: '+str(clipped_length)+'\n')
    
    log_FH.write('working directory: '+os.getcwd()+'\n')

    discord_bam = bam_short_name+'_discord.bam'
    clipped_bam = bam_short_name+'_clipped.bam'
    #clipped_bam_cmpl = bam_short_name+'_clipped_cmpl.bam'
    subprocess.run(['mkdir' , '-p' , preprocess_dir_realpath ])

    samfile = pysam.AlignmentFile(bam_full, "rb")
    newsam_d = pysam.AlignmentFile(preprocess_dir_realpath + '/' + discord_bam, 'wb', template=samfile)
    newsam_c = pysam.AlignmentFile(preprocess_dir_realpath + '/' + clipped_bam, 'wb', template=samfile)
    #newsam_c_cmpl = pysam.AlignmentFile(preprocess_dir_realpath + '/' + clipped_bam_cmpl, 'wb', template=samfile)

    for read in samfile.fetch():
        if read.is_paired == True and read.is_proper_pair != True:
            newsam_d.write(read)
        if read.cigartuples != None:
            if ( (read.cigartuples[-1][0] == 4)
                  and ( read.cigartuples[-1][1] > clipped_length ) 
                  and ((read.cigartuples[0][0] == 4 and read.cigartuples[0][1] > 5) != True)):
                #clipped side = R
                a = pysam.AlignedSegment()
                a.query_name = read.query_name
                #a.query_sequence=read.query_sequence[read.infer_query_length() \
                #            - read.cigartuples[-1][1]-1:read.infer_query_length()]
                start = read.infer_query_length() - read.cigartuples[-1][1] #double checked this
                end = read.infer_query_length()
                a.query_sequence=read.query_sequence[ start : end ]
                a.flag = read.flag
                a.reference_id = read.reference_id
                a.reference_start = read.reference_start
                a.mapping_quality = read.mapping_quality
                #a.mapping_quality = 0 # read.mapping_quality
                #a.cigarstring = []
                #a.cigarstring.append((read.cigartuples[-1][0], read.cigartuples[-1][1]))
                cigar = []
                cigar.append((read.cigartuples[-1][0], read.cigartuples[-1][1])) 
                a.cigartuples = cigar
                #a.cigar = []
                #a.cigar.append((read.cigartuples[-1][0], read.cigartuples[-1][1]))
                a.next_reference_id = read.next_reference_id
                a.next_reference_start=read.next_reference_start
                a.template_length=read.template_length
                #a.query_qualities = read.query_qualities[read.infer_query_length() \
                #            - read.cigartuples[-1][1]-1:read.infer_query_length()]
                a.query_qualities = read.query_qualities[start:end]
                a.set_tags(read.get_tags())
                a.set_tag("ZS",'R')
                write_flag = check_uniq_mapping( read, args )
                a.set_tag("ZU",write_flag)
                a.is_supplementary = read.is_supplementary
                #eprint("top", a.cigar, a.cigarstring, str(len(a.query_sequence)), a.query_name)
                newsam_c.write(a)
                #newsam_c_cmpl.write(read)
            elif ( (read.cigartuples[0][0] == 4 ) 
                      and ( read.cigartuples[0][1] > clipped_length ) 
                      and ((read.cigartuples[-1][0] == 4 and read.cigartuples[-1][1] > 5) != True) ):
                #clipped side = L
                a = pysam.AlignedSegment()
                a.query_name = read.query_name
                #a.query_sequence=read.query_sequence[0:read.cigartuples[0][1]-1]
                a.query_sequence=read.query_sequence[0:read.cigartuples[0][1]]
                a.flag = read.flag
                a.reference_id = read.reference_id
                a.reference_start = read.reference_start
                #a.mapping_quality = 0 #read.mapping_quality
                a.mapping_quality = read.mapping_quality
                cigar = []
                cigar.append((read.cigartuples[0][0], read.cigartuples[0][1])) 
                a.cigartuples = cigar
                #a.cigar = []
                #a.cigar.append((read.cigartuples[0][0], read.cigartuples[0][1])) 
                a.next_reference_id = read.next_reference_id
                a.next_reference_start=read.next_reference_start
                a.template_length=read.template_length
                #a.query_qualities=read.query_qualities[0:read.cigartuples[0][1]-1]
                #a.query_qualities=read.query_qualities[0:read.cigartuples[0][1]]
                a.set_tags(read.get_tags())
                a.set_tag("ZS",'L')
                write_flag = check_uniq_mapping( read, args )
                a.set_tag("ZU",write_flag)
                a.is_supplementary = read.is_supplementary
                #eprint("bot", a.cigarstring, str(len(a.query_sequence)), a.query_name)
                newsam_c.write(a)
                #newsam_c_cmpl.write(read)
    newsam_d.close()
    newsam_c.close()
    #newsam_c_cmpl.close()
    samfile.close()
    pysam.index(preprocess_dir_realpath + '/' + discord_bam)
    pysam.index(preprocess_dir_realpath + '/' + clipped_bam)
    #pysam.index(preprocess_dir_realpath + '/' + clipped_bam_cmpl)
    log_FH.write('Created %s/%s\n' % (preprocess_dir_realpath,discord_bam))
    log_FH.write('Created %s/%s\n' % (preprocess_dir_realpath,clipped_bam))
    #log_FH.write('Created %s/%s\n' % (preprocess_dir_realpath,clipped_bam_cmpl))

    ref_type_file_name = []
    with open(fofn_ref_realpath, 'r') as ref_type_file_file:
        for line in ref_type_file_file:
            ref_type_file_name.append([line.split()[0], line.split()[1]])

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

    log_FH.write('Created ' + preprocess_dir_realpath + '/te_ref_type.fa\n')
    log_FH.write('Created ' + preprocess_dir_realpath + '/te_ref_type_bwa.fa\n')
    log_FH.close()

 
def preprocess_setup_arg_parser(parser):
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument('-i', '--input_bam', action='store', dest='bam_inp', required=True,
        help='input Bam(.bam) file of aligned reads')
    parser_required.add_argument('-r', '--ref', action='store', dest='fofn_ref', required=True,
        help='File with reference sequence paths, see README.md for more info')
    parser.add_argument('-p', '--preprocess_dir', action='store',
        dest='preprocess_dir', default='preprocessed_files',
        help='directory to store preprocessing output files (default: preprocessed_files)')
    parser.add_argument('--min_clipped_len', action='store', dest='cll_inp', type=int, default=25,
        help='Minimum clipped length(bp) (default: 25)')
    parser.add_argument('--min_map_qual', action='store', dest='mpq_inp', type=int, default=30,
        help='Minimum mapping quality (default: 30)')
    parser.add_argument('--log_file', action='store',
        dest='log_file', default='preprocess.log',
        help='run log file (default: preprocess.log)')
    parser._action_groups.reverse()
    return parser   