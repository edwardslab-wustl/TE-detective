import sys, os
import pysam
import numpy as np
from TEdetective.general_functions import break_points_2d

def exec_cluster2d(args):    
    #
    log_FH=open(args.log_file, 'w')
    preprocess_dir_realpath = os.path.realpath(args.preprocess_dir)
    log_FH.write('preprocessing/intermediate file directory: '+ preprocess_dir_realpath +'\n')
    #
    fofn_ref_realpath = os.path.realpath(args.fofn_ref)
    log_FH.write('fofn_ref file: '+ fofn_ref_realpath +'\n')
    #
    dir_path = os.getcwd()
    log_FH.write('working directory: '+ dir_path +'\n')
    #
    log_FH.write('Writing recluster initial predictions to: '+ args.output_file +'\n')
    #
    bam_full = args.bam_inp
    insert_size = args.isz_inp
    discord_rd_clust_denst = args.drd_inp
    read_length = args.rdl_inp
    coverage_cutoff    = args.cct_inp
    ref_type_file_name = []
    with open(fofn_ref_realpath, 'r') as ref_type_file_file:
        for line in ref_type_file_file:
            ref_type_file_name.append([line.split()[0], line.split()[1]])

    #read_positions_clusters_file = open('recluster_initial_predictions.txt', 'w')
    read_positions_clusters_file = open(args.output_file, 'w')
    for cnt_1 in range(0, len(ref_type_file_name)):
        flag_read_position = 'y'
        read_positions = []
        #
        if args.flg_all:
            with open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_read-bam_mate_id_pos.dat', 'r') as mate_id_dat:
                mate_id_dat_lines = mate_id_dat.readlines()
            mate_id_dat.close()

            for line in mate_id_dat_lines:
                read_positions.append((line.split()[2], int(line.split()[3])))
            del mate_id_dat_lines[:]
        #
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
    log_FH.close()

def cluster2d_setup_arg_parser(parser):
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument('-i','--input_bam', action='store', dest='bam_inp', required=True, 
        help='input Bam(.bam) file of aligned reads')
    parser_required.add_argument('-r', '--ref', action='store', dest='fofn_ref', required=True,
        help='File with reference sequence paths, see README.md for more info')
    parser.add_argument('-o', '--output_file', action='store',
        dest='output_file', default='recluster_initial_predictions.txt',
        help='Tab-delimited file of initial set of TE insertions (default: recluster_initial_predictions.txt)')
    parser.add_argument('-p', '--preprocess_dir', action='store',
        dest='preprocess_dir', default='preprocessed_files',
        help='directory used to store preprocessing output files (default: preprocessed_files)')
    parser.add_argument('--insert_size_est', action='store', dest='isz_inp', type=int, default=340,
        help='insert Size estimate (default: 340)')
    parser.add_argument('--read_length', action='store', dest='rdl_inp', type=int, default=150, 
        help='Average read length (default: 150)')
    parser.add_argument('--discord_cluster_dens', action='store', dest='drd_inp', type=int, default=5, 
        help='Discord read cluster density (default: 5)')
    parser.add_argument('--coverage_cutoff', action='store', dest='cct_inp', type=int, default=200, 
        help='Coverage cutoff input (default: 200)')
    parser.add_argument('--all', action='store_true', dest='flg_all', default=False, 
        help='use all reads instead of only clipped (default: False)')
    parser.add_argument('--log_file', action='store',
        dest='log_file', default='cluster2d.log',
        help='run log file (default: cluster2d.log)')
    parser._action_groups.reverse()
    return parser