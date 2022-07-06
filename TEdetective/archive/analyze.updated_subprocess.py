import os, subprocess
import pysam
from TEdetective.io_functions import check_file, read_type_info, te_type_setup, read_class_info, get_class, print_tup, num_of_lines
from TEdetective.general_functions import break_points, check_uniq_mapping, calc_length
from TEdetective.analyze_functions import find_clipped_ends,flt_discord, analysis_pat_check


def exec_analyze(args):
    #
    log_FH=open(args.log_file, 'w')
    preprocess_dir_realpath = os.path.realpath(args.preprocess_dir)
    log_FH.write('preprocessing/intermediate file directory: '+ preprocess_dir_realpath +'\n')
    #
    fofn_ref_realpath = os.path.realpath(args.fofn_ref)
    log_FH.write('preprocessing/intermediate file directory: '+ preprocess_dir_realpath +'\n')
    #
    dir_path = os.getcwd()
    log_FH.write('working directory: '+ dir_path +'\n')
    #up_dir_path = os.path.dirname(os.path.realpath(dir_path))
    #
    bam_full = args.bam_inp
    log_FH.write('Input bam file: '+bam_full+'\n')
    #
    bam_short_name = bam_full.split('/')[-1][:-4]
    log_FH.write('bam short name: '+ bam_short_name +'\n')
    #
    discord_bam = preprocess_dir_realpath + '/'+bam_short_name+'_discord.bam'
    log_FH.write('Discordant bam file: '+discord_bam+'\n')
    if not check_file(discord_bam + ".bai"):
        log_FH.write('cannot find index for: '+discord_bam+', making now\n')
        pysam.index(discord_bam)
    #
    te_type_file = preprocess_dir_realpath + '/te_ref_type.fa'
    #
    log_FH.write('\n----Run parameters----\n')
    #
    insert_range = args.rdl_inp # =read_length
    log_FH.write('Input read length: '+str(insert_range)+'\n')    
    #
    clipped_length = args.cll_inp
    log_FH.write('Minimum clipped read length: '+str(clipped_length)+'\n')
    #
    anchor_length = args.ahl_inp
    log_FH.write('Minimum anchor read length: '+str(anchor_length)+'\n')
    #
    cliped_end_range = args.cer_inp
    log_FH.write('Range to put clipped read in group: '+str(cliped_end_range)+'\n')
    #
    clipped_search_interval = args.csi_inp
    log_FH.write('Clipped search interval: '+str(clipped_search_interval)+'\n')    
    #
    min_reads_het = args.mrh_inp
    log_FH.write('Minimum reads to call hetrozygous: '+str(min_reads_het)+'\n')
    #
    insert_size = args.isz_inp
    log_FH.write('Insert size estimate: '+str(insert_size)+'\n')
    #
    min_reads_ends = args.mre_inp
    log_FH.write('Minimum read at each end: '+str(min_reads_ends)+'\n') # ??
    #
    qual_interval_inp = args.qii_inp
    log_FH.write('Mapping quality gap: '+str(qual_interval_inp)+'\n')
    #
    num_interval_inp = args.nii_inp
    log_FH.write('Number of quality gap searchs: '+str(num_interval_inp)+'\n')
    #
    min_mapq = args.mpq_inp
    log_FH.write('minimum mapping quality: '+str(min_mapq)+'\n')
    #
    min_mapq_uniq = args.mpqu_inp
    log_FH.write('minimum mapping quality for uniq test: '+str(min_mapq_uniq)+'\n\n')
    #
    list_file = args.list_inp
    #
    # Read initial prediction list
    with open(list_file, 'r') as inp_file:
        inp_file_lines = inp_file.readlines()
    inp_file.close()
    #
    #
    #eprint("----------\n"  + "remap and find insertion points\n" + "----------\n")
    log_FH.write("----------\n"  + "remap and find insertion points\n" + "----------\n")
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
            #eprint(warning)
            log_FH.write(warning + "\n")
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
            log_FH.write("censor p-end: " + int_file_name + "\n")
            #subprocess.run(['censor.ncbi', int_file_name+'_reads_p.fa' , '-lib' , te_type_file])
            result = subprocess.run(['censor.ncbi', int_file_name+'_reads_p.fa' , '-lib' , te_type_file], capture_output=True, text=True)
            log_FH.write(result.stdout)
            log_FH.write(result.stderr)
            log_FH.write("\n\n")
        if cnt_rd_n > 0: # Number of clipped reads at n-end of TE >1
#            proc = Popen(['censor.ncbi', int_file_name+'_reads_n.fa' , '-lib' , te_type_file], shell=True, stdout=PIPE, stderr=PIPE)
#            out, err = proc.communicate()
            log_FH.write("censor n-end: " + int_file_name + "\n")
            #subprocess.run(['censor.ncbi', int_file_name+'_reads_n.fa' , '-lib' , te_type_file])
            result = subprocess.run(['censor.ncbi', int_file_name+'_reads_n.fa' , '-lib' , te_type_file], capture_output=True, text=True)
            log_FH.write(result.stdout)
            log_FH.write(result.stderr)
            log_FH.write("\n\n")
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
            #subprocess.run(['censor.ncbi', int_file_name+'_reads_p.fa' , '-lib' , te_class_file])
            #subprocess.run(['censor.ncbi', int_file_name+'_reads_n.fa' , '-lib' , te_class_file])
            log_FH.write("censor p reads: " + int_file_name + "\n")
            result = subprocess.run(['censor.ncbi', int_file_name+'_reads_p.fa' , '-lib' , te_class_file], capture_output=True, text=True)
            log_FH.write(result.stdout)
            log_FH.write(result.stderr)
            log_FH.write("\n\n")
            log_FH.write("censor n reads: " + int_file_name + "\n")
            result = subprocess.run(['censor.ncbi', int_file_name+'_reads_n.fa' , '-lib' , te_class_file], capture_output=True, text=True)
            log_FH.write(result.stdout)
            log_FH.write(result.stderr)
            log_FH.write("\n\n")
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
                log_FH.write("censor n clipped: " + int_file_name + "\n")
                result = subprocess.run(['censor.ncbi' , 'n.clipped.fa' , '-lib' , class_out_p[0][0]+'.fa'], capture_output=True, text=True)
                log_FH.write(result.stdout)
                log_FH.write(result.stderr)
                log_FH.write("\n\n")
                #
                if check_file('n.clipped.fa.map') == True:
                    output_file.write(str(class_out_p[0][0]) + ' ' + str( get_class('n.clipped.fa.map', args)[0][0] ) + ' ' + \
                                    str( calc_length( class_out_p, get_class( 'n.clipped.fa.map', args ) ) ) + '\n') 
                    len_arr.append([class_out_p[0][0], calc_length(class_out_p, get_class('n.clipped.fa.map', args))])
            
                subprocess.run(['cp' , int_file_name+'_reads_p.fa' , 'p.clipped.fa'])
                #subprocess.run(['censor.ncbi' , 'p.clipped.fa' , '-lib' , class_out_n[0][0]+'.fa'])
                log_FH.write("censor n clipped: " + int_file_name + "\n")
                result = subprocess.run(['censor.ncbi' , 'p.clipped.fa' , '-lib' , class_out_n[0][0]+'.fa'], capture_output=True, text=True)
                log_FH.write(result.stdout)
                log_FH.write(result.stderr)
                log_FH.write("\n\n")
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
            output_file.write(str(type_clipped_p[0][0])+'\t'+ str(type_clipped_p[0][1])+'\t'+  str(type_clipped_p[0][2]) +'\t'+ type_c_p_out +'\t')
        if clipped_read_p_flag == 'n':
            output_file.write("NA\tNA\t0\t0\t")
        #
        if clipped_read_n_flag == 'y':
            if cnt_rd_n > 0:
                type_c_n_out = str(float("{0:.2f}".format((type_clipped_n[0][2]/cnt_rd_n)*100)))
            else:
                type_c_n_out = 'NA'
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
    #eprint("----------\n"  + "discordant mate pair analysis\n" + "----------\n")
    log_FH.write("----------\n"  + "discordant mate pair analysis\n" + "----------\n")
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
            #subprocess.run(['censor.ncbi' , int_file_name+'_discord_mate_p.fa' , '-lib' , te_type_file])
            log_FH.write("censor p discord: " + int_file_name + "\n")
            result = subprocess.run(['censor.ncbi' , int_file_name+'_discord_mate_p.fa' , '-lib' , te_type_file], capture_output=True, text=True)
            log_FH.write(result.stdout)
            log_FH.write(result.stderr)
            log_FH.write("\n\n")
        if cnt_dm_n > 0:
            #subprocess.run(['censor.ncbi' , int_file_name+'_discord_mate_n.fa' , '-lib' , te_type_file])
            log_FH.write("censor n discord: " + int_file_name + "\n")
            result = subprocess.run(['censor.ncbi' , int_file_name+'_discord_mate_n.fa' , '-lib' , te_type_file], capture_output=True, text=True)
            log_FH.write(result.stdout)
            log_FH.write(result.stderr)
            log_FH.write("\n\n")
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
            #subprocess.run(['censor.ncbi' , int_file_name+'_discord_mate_p.fa' , '-lib' , te_class_file])
            #subprocess.run(['censor.ncbi' , int_file_name+'_discord_mate_n.fa' , '-lib' , te_class_file])
            log_FH.write("censor p discord mates: " + int_file_name + "\n")
            result = subprocess.run(['censor.ncbi' , int_file_name+'_discord_mate_p.fa' , '-lib' , te_class_file], capture_output=True, text=True)
            log_FH.write(result.stdout)
            log_FH.write(result.stderr)
            log_FH.write("\n\n")
            log_FH.write("censor n discord mates: " + int_file_name + "\n")
            result = subprocess.run(['censor.ncbi' , int_file_name+'_discord_mate_n.fa' , '-lib' , te_class_file], capture_output=True, text=True)
            log_FH.write(result.stdout)
            log_FH.write(result.stderr)
            log_FH.write("\n\n")
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
    #eprint("----------\n"  + "writing output files\n" + "----------\n")
    log_FH.write("----------\n"  + "writing output files\n" + "----------\n")
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
    log_FH.close()