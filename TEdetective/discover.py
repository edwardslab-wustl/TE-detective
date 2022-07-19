import sys
import os
import subprocess

import pysam
import numpy as np
from Bio.Sequencing.Applications import BwaIndexCommandline
from Bio.Sequencing.Applications import BwaAlignCommandline
from Bio.Sequencing.Applications import BwaSamseCommandline

from TEdetective.io_functions import eprint
from TEdetective.general_functions import check_uniq_mapping, break_points_2d

def exec_discover(args):
    
    log_FH=open(args.log_file, 'w')
    dir_path = os.getcwd() #os.path.dirname(os.path.realpath(__file__))
    log_FH.write('Working directory: '+str(dir_path)+'\n')
    
    preprocess_dir_realpath = os.path.realpath(args.preprocess_dir)
    log_FH.write('preprocessing/intermediate file directory: '+ preprocess_dir_realpath +'\n')
    
    fofn_ref_realpath = os.path.realpath(args.fofn_ref)
    log_FH.write('fofn_ref file: '+ fofn_ref_realpath +'\n')
    
    bam_full = os.path.realpath(args.bam_inp)
    log_FH.write('bam file name: '+str(bam_full)+'\n')
    
    bam_short_name = bam_full.split('/')[-1][:-4]
    log_FH.write('bam short name: '+ bam_short_name +'\n')
    
    insert_size = args.isz_inp
    log_FH.write('Insert size estimate: '+str(insert_size)+'\n')
    
    discord_rd_clust_denst = args.drd_inp # <--change this variable name, confusing.
    log_FH.write('Number of reads in a cluster to call it insertion: '+str(discord_rd_clust_denst)+'\n')
    
    read_length = args.rdl_inp
    log_FH.write('Average read length: '+str(read_length)+'\n')
    
    coverage_cutoff = args.cct_inp
    log_FH.write('Coverage cutoff to skip a region: '+str(coverage_cutoff)+'\n')
    
    min_mapq = args.mpq_inp
    log_FH.write('minimum mapping quality: '+str(min_mapq)+'\n')
    
    min_mapq_uniq = args.mpqu_inp
    log_FH.write('minimum mapping quality: '+str(min_mapq_uniq)+'\n')
        
    clipped_length = args.cll_inp
    log_FH.write('Minimum clipped length: '+str(clipped_length)+'\n')
    
    log_FH.write('Writing initial predictions to: '+ args.output_file +'\n')
    
    subprocess.run(['mkdir' , '-p' , preprocess_dir_realpath])
    
    reference_genome = preprocess_dir_realpath+'/te_ref_type_bwa.fa'
    index_cmd = BwaIndexCommandline(infile=reference_genome, algorithm='bwtsw')
    index_cmd()
    
    read_bam = preprocess_dir_realpath + '/'+bam_short_name+'_discord.bam'
    output_sai_file = preprocess_dir_realpath + '/aligned_reads.sai'
    read_group='@RG ID:foo  SM:bar'
    align_cmd = BwaAlignCommandline(reference=reference_genome, b='b', read_file=read_bam)
    align_cmd(stdout=output_sai_file)
    
    output_sam_file = preprocess_dir_realpath + '/aligned_reads.sam'
    samse_cmd = BwaSamseCommandline(reference=reference_genome, read_file=read_bam, sai_file=output_sai_file)
    samse_cmd(stdout=output_sam_file)
    
    read_bam_clipped = preprocess_dir_realpath+'/'+bam_short_name+'_clipped.bam'
    output_sai_file_clipped = preprocess_dir_realpath+ '/aligned_reads_clipped.sai'
    read_group='@RG ID:foo  SM:bar'
    align_cmd_clipped = BwaAlignCommandline(reference=reference_genome, b='b', read_file=read_bam_clipped)
    align_cmd_clipped(stdout=output_sai_file_clipped)
    
    output_sam_file_clipped = preprocess_dir_realpath + '/aligned_reads_clipped.sam'
    samse_cmd_clipped = BwaSamseCommandline(reference=reference_genome, read_file=read_bam_clipped, sai_file=output_sai_file_clipped)
    samse_cmd_clipped(stdout=output_sam_file_clipped)    
    
    ref_type_file_name = []
    with open(fofn_ref_realpath, 'r') as ref_type_file_file:
        for line in ref_type_file_file:
            ref_type_file_name.append([line.split()[0], line.split()[1]])
    #
    # n is 0 based *********
    dict_aof = {}
    for cnt_1 in range( len(ref_type_file_name) ):
        dict_aof[ref_type_file_name[cnt_1][0]] = []
    
    with open(output_sam_file,'r') as output_sam_file_fl:
        output_sam_file_fl_lines = output_sam_file_fl.readlines()
    output_sam_file_fl.close
    
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
    
    for cnt_1 in range( len(ref_type_file_name) ):
        with open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_ref-aligned_line-num_id.dat', 'w') as temp_fl:
            temp_fl.write('\n'.join(dict_aof[ref_type_file_name[cnt_1][0]]))    
        temp_fl.close
    dict_aof.clear()
    
    ######
    dict_aofc = {}
    for cnt_1 in range( len(ref_type_file_name) ):
        dict_aofc[ref_type_file_name[cnt_1][0]] = []

    with open(output_sam_file_clipped,'r') as output_sam_file_clipped_fl:
        output_sam_file_clipped_fl_lines = output_sam_file_clipped_fl.readlines()
    output_sam_file_clipped_fl.close()
    
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
    
    for cnt_1 in range( len(ref_type_file_name) ):
        with open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_clipped_ref-aligned_line-num_id.dat', 'w') as temp_fl:
            temp_fl.write('\n'.join(dict_aofc[ref_type_file_name[cnt_1][0]]))
        temp_fl.close
    dict_aofc.clear()
    
    # Extract ids of mapped discordant reads.
    for cnt_1 in range( len(ref_type_file_name) ):
        line_number = 0
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
        raw_out_file = open(preprocess_dir_realpath +'/'+ref_type_file_name[cnt_1][0]+'_read-bam_id_flag.dat','w')
        raw_out_file.write('\n'.join(raw_out_file_line))
        del raw_out_file_line[:]
        raw_out_file.close()
        dict_rali.clear()

#    Extract ids of uniquly mapped clipped reads
    for cnt_1 in range( len(ref_type_file_name) ):
        line_number = 0
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
            if line_number in dict_crali:
                #eprint(line_number, read.query_name,
                #       dict_crali[line_number], read.is_supplementary,
                #       read.cigartuples, read.mapping_quality, read.cigarstring)
                if ( (dict_crali[line_number] == read.query_name)
                    and (read.is_supplementary != True)
                    and (read.cigartuples != None)
                    and (read.has_tag('XA') or read.has_tag('SA') or read.mapping_quality >= min_mapq)
                    and (read.mapping_quality >= min_mapq_uniq) ):

                    write_flag = check_uniq_mapping( read, args )
                    #eprint(line_number, read.query_name, write_flag)
                    if write_flag == 'y':
                        clipped_side = 'X'
                        if  ( read.cigartuples[-1][0] == 4 ) and ( read.cigartuples[-1][1] > clipped_length ):
                            clipped_side = 'R'
                        elif ( read.cigartuples[0][0] == 4 ) and ( read.cigartuples[0][1] > clipped_length ):
                            clipped_side = 'L'
                        raw_out_file_line.append(read.query_name +' '+ str(read.flag) + ' ' \
                            + str(read.reference_name) + ' ' + str(read.reference_start) + ' ' + str(read.reference_end) + ' ' + clipped_side)
        samfile.close()
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
        mate_out_file = open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_read-bam_mate_id_pos.dat', 'w')
        mate_out_file.write('\n'.join(mate_out_file_line))
        del mate_out_file_line[:]
        mate_out_file.close()
    samfile_idx.close()
    
    #read_positions_clusters_file = open('initial_predictions.txt', 'w')
    read_positions_clusters_file = open(args.output_file, 'w')
    for cnt_1 in range( len(ref_type_file_name) ):
        flag_read_position = 'y'
        read_positions = []
        with open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_read-bam_mate_id_pos.dat', 'r') as mate_id_dat:
            mate_id_dat_lines = mate_id_dat.readlines()
        mate_id_dat.close()

        for line in mate_id_dat_lines:
            read_positions.append((line.split()[2], int(line.split()[3])))
        del mate_id_dat_lines[:]

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
                #print(np.mean(coverage_values_np[read_length:-1*read_length]))
                mean_coverage = np.mean(coverage_values_np[read_length:-1*read_length])
                log_FH.write( bam_short_name + ".bam mean coverage: " + str(mean_coverage) + "\n")
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
    log_FH.close()
    
def discover_setup_arg_parser(parser):
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument('-i', '--input_bam', action='store', dest='bam_inp', required=True, 
        help='Input Bam(.bam) file of aligned reads')
    parser_required.add_argument('-r', '--ref', action='store', dest='fofn_ref', required=True,
        help='File with reference sequence paths, see README.md for more info')
    parser.add_argument('-o', '--output_file', action='store',
        dest='output_file', default='initial_predictions.txt',
        help='Tab-delimited file of initial set of TE insertions (default: initial_predictions.txt)')
    parser.add_argument('-p', '--preprocess_dir', action='store',
        dest='preprocess_dir', default='preprocessed_files',
        help='directory used to store preprocessing output files (default: preprocessed_files)')
    parser.add_argument('--insert_size_est', action='store', dest='isz_inp', type=int, default=340, 
        help='Insert size estimate (default: 340)')
    parser.add_argument('--read_length', action='store', dest='rdl_inp', type=int, default=150, 
        help='Average read length (default: 150)')
    parser.add_argument('--discord_cluster_dens', action='store', dest='drd_inp', type=int, default=10, 
        help='Discord read cluster density (default: 10)')
    parser.add_argument('--coverage_cutoff', action='store', dest='cct_inp', type=int, default=200, 
        help='Coverage cutoff input (default: 200)')
    parser.add_argument('--min_clipped_len', action='store', dest='cll_inp', type=int, default=25,
        help='Minimum clipped length(bp) (default: 25)')
    parser.add_argument('--min_map_qual', action='store', dest='mpq_inp', type=int, default=30, 
        help='Minimum mapping quality (default: 30)')
    parser.add_argument('--map_qual_uniq', action='store', dest='mpqu_inp', type=int, default=1, 
        help='Minimum mapping quality (default: 1)')
    parser.add_argument('--log_file', action='store',
        dest='log_file', default='discover.log',
        help='run log file (default: discover.log)')
    parser._action_groups.reverse()
    return parser