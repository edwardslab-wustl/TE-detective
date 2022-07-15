import os

import pysam
import numpy as np

from TEdetective.general_functions import check_uniq_mapping, break_points_2d, pat_check, break_points
from TEdetective.nadiscover_functions import alt_mapped_pos, write_pat_clipped_dat, write_pat_clipped_dat_new
from TEdetective.io_functions import eprint

def exec_nadiscover(args):
    log_FH=open(args.log_file, 'w')
    dir_path = os.getcwd() #os.path.dirname(os.path.realpath(__file__))
    log_FH.write('Working directory: '+str(dir_path)+'\n')
    
    preprocess_dir_realpath = os.path.realpath(args.preprocess_dir)
    log_FH.write('preprocessing/intermediate file directory: '+ preprocess_dir_realpath +'\n')
    
    fofn_ref_realpath = os.path.realpath(args.fofn_ref)
    log_FH.write('fofn_ref file: '+ fofn_ref_realpath +'\n')
    
    bam_full = args.bam_inp
    log_FH.write('bam file name: '+str(bam_full)+'\n')
    
    bam_short_name = bam_full.split('/')[-1][:-4]
    log_FH.write('bam short name: '+ bam_short_name +'\n')
    
    insert_size = args.isz_inp
    log_FH.write('Insert size estimate: '+str(insert_size)+'\n')
    
    discord_rd_clust_denst = args.drd_inp
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
    
    log_FH.write('Writing initial predictions noalign to: '+ args.output_file +'\n')
    
    pat_query_len = args.pql_inp
    pat_mis_match = args.pmm_inp
    rmsk_bed = args.rmsk_bed
    read_bam_clipped = preprocess_dir_realpath +'/'+bam_short_name+'_clipped.bam'
    
    if args.pat_inp:
        pat_out_file = open(preprocess_dir_realpath+'/pat_clipped_read-bam_id_flag.dat','w')
        #pat_out_file_lines = write_pat_clipped_dat(read_bam_clipped, min_mapq_uniq, clipped_length, pat_out_file, pat_query_len, pat_mis_match, args)
        pat_out_file_lines = write_pat_clipped_dat_new(bam_full, read_bam_clipped, min_mapq_uniq, clipped_length, pat_out_file, pat_query_len, pat_mis_match, args)
        
    ref_type_file_name = []
    with open(fofn_ref_realpath, 'r') as ref_type_file_file:
        for line in ref_type_file_file:
            ref_type_file_name.append([line.split()[0], line.split()[1]])
    
    if args.nas_inp:
        # Process RMSK file.
        dict_ref_bed = {}
        for items in ref_type_file_name:
            dict_ref_bed[items[0]] = []
        #
        with open(rmsk_bed, 'r') as rmsk_bed_file:
            rmsk_bed_file_lines = rmsk_bed_file.readlines()
        rmsk_bed_file.close()
        
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
        open(preprocess_dir_realpath+'/rmsk_cons.bed', 'w').write(''.join(rmsk_bed_items))
        
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
        discord_bam = preprocess_dir_realpath + '/'+bam_short_name+'_discord.bam'
        log_FH.write('Discordant bam file: '+discord_bam+'\n')
        
        samfile_discord = pysam.AlignmentFile(discord_bam, 'rb')
        for read in samfile_discord.fetch():
            discord_pos_list.append( read.query_name +'\t'+ str(read.flag) +'\t'+ read.reference_name +'\t'+ str(read.reference_start) +'\t'+ str(read.reference_end)+'\n')
        samfile_discord.close()
        open(preprocess_dir_realpath+'/discord_pos_list.txt', 'w').write(''.join(discord_pos_list))
        
        if args.flg_all:
            clipped_pos_list = []
            clipped_bam = preprocess_dir_realpath + '/'+bam_short_name+'_clipped.bam'
            log_FH.write('Clipped bam file: '+clipped_bam+'\n')
            
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

            open(preprocess_dir_realpath+'/clipped_pos_list.txt', 'w').write(''.join(clipped_pos_list))
            samfile_clipped.close()
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
        
        for cnt_1 in range(0, len(ref_type_file_name)):
            try:
                open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_noalign_read-bam_id_flag.dat', 'w').write(''.join(dict_discord[ref_type_file_name[cnt_1][0]]))
            except KeyError:
                open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_noalign_read-bam_id_flag.dat', 'w').write('')
        
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
            for cnt_1 in range(0, len(ref_type_file_name)):
                try:
                    open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_noalign_clipped_read-bam_id_flag.dat', 'w').write(''.join(dict_clipped[ref_type_file_name[cnt_1][0]]))
                except KeyError:
                    open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_noalign_clipped_read-bam_id_flag.dat', 'w').write('')
        
        # Search for mate and their position
        samfile_idx = pysam.AlignmentFile(discord_bam, 'rb')
        id_index = pysam.IndexedReads(samfile_idx)
        id_index.build()
        
        for cnt_1 in range( len(ref_type_file_name) ):
            with open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_noalign_read-bam_id_flag.dat' ,'r') as read_bam_dat:
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
            mate_out_file = open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_noalign_read-bam_mate_id_pos.dat', 'w')
            mate_out_file.write('\n'.join(mate_out_file_line))
            del mate_out_file_line[:]
            mate_out_file.close()
        samfile_idx.close()
        # nas condition done
        
    read_positions_clusters_file = open(args.output_file, 'w')
    for cnt_1 in range( len(ref_type_file_name) ):
        read_positions = []
        dict_pat_test = {}
        
        # nas condition
        if args.nas_inp:
            with open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_noalign_read-bam_mate_id_pos.dat', 'r') as noalign_mate_id_dat:
                noalign_mate_id_dat_lines = noalign_mate_id_dat.readlines()
            noalign_mate_id_dat.close()
            #
            if args.merged:
                with open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_read-bam_mate_id_pos.dat', 'r') as align_mate_id_dat:
                    align_mate_id_dat_lines = align_mate_id_dat.readlines()
                align_mate_id_dat.close()
                dict_align = {}
                align_mate_read_positions = []
                
                for line in align_mate_id_dat_lines:
                    align_mate_read_positions.append( ( line.split()[2], int(line.split()[3]) ) )
                    item = line.strip().split()[0]+'_'+line.strip().split()[1]+'_'+line.strip().split()[2]+'_'+line.strip().split()[3]
                    dict_align[item] = 1
                del align_mate_id_dat_lines[:]
                noalign_mate_read_positions = []
                
                for line in noalign_mate_id_dat_lines:
                    test_item = line.strip().split()[0]+'_'+line.strip().split()[1]+'_'+line.strip().split()[2]+'_'+line.strip().split()[3]
                    try:
                        if dict_align[test_item] == 1:
                            continue
                    except KeyError:
                        noalign_mate_read_positions.append( ( line.split()[2], int(line.split()[3]) ) )
                    
                # del noalign_mate_id_dat_lines[:] need this if no merged
                mate_read_positions = align_mate_read_positions + noalign_mate_read_positions
                #print(mate_read_positions)
                #print(align_mate_read_positions)
                #print(noalign_mate_read_positions)
                del align_mate_read_positions[:]
                del noalign_mate_read_positions[:]
                dict_align.clear()
                
                # Clipped
                with open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_clipped_read-bam_id_flag.dat', 'r') as align_id_dat:
                    align_id_dat_lines = align_id_dat.readlines()
                align_id_dat.close()
                dict_clp_align = {}
                align_read_positions = []
                for line in align_id_dat_lines:
                    align_read_positions.append( ( line.split()[2], int(line.split()[3]) ) )
                    item = line.strip().split()[0]+'_'+line.strip().split()[1]+'_'+line.strip().split()[2]+'_'+line.strip().split()[3]
                    dict_clp_align[item] = 1
                    dict_pat_test[line.strip().split()[0]] = line.strip().split()[1]
                del align_id_dat_lines[:]

                noalign_read_positions = []
                if args.flg_all:
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
                del align_read_positions[:]
                del noalign_read_positions[:]
                dict_clp_align.clear()
                read_positions = mate_read_positions + clpd_read_positions
                del mate_read_positions[:]
                del clpd_read_positions[:]
            if not args.merged:
                mate_read_positions = []
                for line in noalign_mate_id_dat_lines:
                    mate_read_positions.append( ( line.split()[2], int(line.split()[3]) ) )

                clpd_read_positions = []
                if args.flg_all:
                    with open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_noalign_clipped_read-bam_id_flag.dat', 'r') as noalign_id_dat:
                        noalign_id_dat_lines = noalign_id_dat.readlines()
                    noalign_id_dat.close()
                    for line in noalign_id_dat_lines:
                        clpd_read_positions.append( ( line.strip().split()[2], int(line.strip().split()[3]) ) )
                        dict_pat_test[line.strip().split()[0]] = line.strip().split()[1]
                    del noalign_id_dat_lines[:]
                    
                read_positions = mate_read_positions + clpd_read_positions
                del mate_read_positions[:]
                del clpd_read_positions[:]
                
            del noalign_mate_id_dat_lines[:]
            # nas consition ends
            
        if not args.nas_inp:
            #print(args.nas_inp)
            log_FH.write('--nonaligned_search flag is ' + str(args.nas_inp))
                
            with open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_read-bam_mate_id_pos.dat', 'r') as mate_id_dat:
                mate_id_dat_lines = mate_id_dat.readlines()
            mate_id_dat.close()

            for line in mate_id_dat_lines:
                read_positions.append((line.strip().split()[2], int(line.strip().split()[3])))
            del mate_id_dat_lines[:]
            
            with open(preprocess_dir_realpath+'/'+ref_type_file_name[cnt_1][0]+'_clipped_read-bam_id_flag.dat', 'r') as clpd_read_id_dat:
                clpd_read_id_dat_lines = clpd_read_id_dat.readlines()
            clpd_read_id_dat.close()
            
            for line in clpd_read_id_dat_lines:
                read_positions.append((line.strip().split()[2], int(line.strip().split()[3])))
                dict_pat_test[line.strip().split()[0]] = line.strip().split()[1]            
    
            del clpd_read_id_dat_lines[:]
        read_positions_sorted = sorted(read_positions, key=lambda x: (x[0], x[1]))
        
        if args.pat_inp: # polyA/T search
            min_non_pat_rd_denst = 3 #min_non_pat_rd_denst
            pat_dat_lines = []
            for line in pat_out_file_lines:
                try:
                    if dict_pat_test[ line.strip().split()[0] ] ==  line.strip().split()[1]:
                        continue
                except (KeyError, ValueError):
                    pat_dat_lines.append( ( line.strip().split()[2], int(line.strip().split()[3]) ) )
            dict_pat_test.clear()
            pat_dat_lines_sorted = sorted(pat_dat_lines, key=lambda x: (x[0], x[1]))
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
                
                line_num += 1
                if pos_track != int( int(items[1]) /10000 ):
                    for i in range(pos_track, int( int(items[1])/10000)):
                        if pos_track + 1 == int( int(items[1])/10000):
                            dict_pat[chr_track+'_idx'][ int( int(items[1])/10000) ] = line_num
                            pos_track = int( int(items[1])/10000 )
                        if pos_track + 1 < int( int(items[1])/10000):
                            dict_pat[chr_track+'_idx'][pos_track+1] = dict_pat[chr_track+'_idx'][pos_track]
                            pos_track = pos_track + 1
            
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
            read_positions_clusters = read_positions_clusters_pat.copy()
            del read_positions_clusters_pat[:]
            dict_pat.clear()
            
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
    log_FH.close()
    
def nadiscover_setup_arg_parser (parser):
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument('-i','--input_bam', action='store', dest='bam_inp', required=True, 
        help='input Bam(.bam) file of aligned reads')
    parser_required.add_argument('-r', '--ref', action='store', dest='fofn_ref', required=True,
        help='File with reference sequence paths, see README.md for more info')
    parser.add_argument('--bed', action='store', dest='rmsk_bed', help='FoFn for existing repeat elements')
    parser.add_argument('-o', '--output_file', action='store',
        dest='output_file', default='initial_predictions_noalign.txt',
        help='Tab-delimited output file of initial set of TE insertions (default: initial_predictions_noalign.txt)')
    parser.add_argument('-p', '--preprocess_dir', action='store',
        dest='preprocess_dir', default='preprocessed_files',
        help='directory used to store preprocessing output files (default: preprocessed_files)')
    parser.add_argument('--min_clipped_len', action='store', dest='cll_inp', type=int, default=25,
        help='Minimum clipped length(bp) (default: 25)')
    parser.add_argument('--insert_size_est', action='store', dest='isz_inp', type=int, default=340,
        help='insert Size estimate (default: 340)')
    parser.add_argument('--read_length', action='store', dest='rdl_inp', type=int, default=150,
        help='Average read length (default: 150)')
    parser.add_argument('--discord_cluster_dens', action='store', dest='drd_inp', type=int, default=5,
        help='discord read cluster density (default: 5)')
    parser.add_argument('--coverage_cutoff', action='store', dest='cct_inp', type=int, default=200,
        help='Coverage cutoff input (default: 200)')
    parser.add_argument('--all', action='store_true', dest='flg_all', default=False,
        help='use all reads instead of only clipped (default: False)')
    parser.add_argument('--merge_aligned', action='store_true', dest='merged', default=False,
        help='Merge aligned predictions (default: False)')
    parser.add_argument('--nonaligned_search', action='store_true', dest='nas_inp', default=False,
        help='Perform non-alignment ref bam search (default: False)')
    parser.add_argument('--min_map_qual', action='store', dest='mpq_inp', type=int, default=30,
        help='Minimum mapping quality (default: 30)')
    parser.add_argument('--map_qual_uniq', action='store', dest='mpqu_inp', type=int, default=1,
        help='Minimum mapping quality unique test (default: 1)')
    parser.add_argument('--polyA', action='store_true', dest='pat_inp', default=False,
        help='Perform poly A/T search (default: False)')
    parser.add_argument('--polyA_len', action='store', dest='pql_inp', type=int, default=9,
        help='poly A/T Length (default: 9)')
    parser.add_argument('--polyA_mismatch', action='store', dest='pmm_inp', type=int, default=1,
        help='poly A/T mismatch (default: 1)')
    parser.add_argument('--log_file', action='store',
        dest='log_file', default='nadiscover.log',
        help='run log file (default: nadiscover.log)')
    parser._action_groups.reverse()
    return parser
    