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
    bam_full = args.bam_inp
    insert_size = args.isz_inp
    discord_rd_clust_denst = args.drd_inp
    read_length = args.rdl_inp
    coverage_cutoff    = args.cct_inp
    ref_type_file_name = []
    with open(fofn_ref_realpath, 'r') as ref_type_file_file:
        for line in ref_type_file_file:
            ref_type_file_name.append([line.split()[0], line.split()[1]])

    read_positions_clusters_file = open('recluster_initial_predictions.txt', 'w')
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