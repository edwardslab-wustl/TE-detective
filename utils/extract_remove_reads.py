#!/home/mksingh/anaconda3/bin/python

def exec_extract_reads():

	# Step 1. Extract discordant reads
	# Step 2. Extract region reads
	# Step 3. Find mates of extracted region reads.
	# Step 4. Remove duplicate reads 

	read_range = args.rng_inp
	# ./tmp dir for sorting
	if not os.path.isdir('tmp'):
		os.makedirs('tmp')
	#
	bam_full = args.bam_inp	
	reg_file = args.reg_inp
	#
	with open(reg_file, 'r') as reg_file:
		reg_file_line = reg_file.readlines()
	reg_file.close()
	#	
	# Readin region info	
	reg_info = []
	for item in reg_file_line:
		reg_info.append( (item.strip().split()[1], item.strip().split()[2]) )

	samfile = pysam.AlignmentFile(bam_full, "rb")
	#
	# extract discordant reads
	discord_bam = bam_full.split('/')[-1][:-4]+'_discord.bam'
	newsam_d = pysam.AlignmentFile(discord_bam, "wb", template=samfile)
	for read in samfile.fetch():
		if read.is_paired == True and read.is_proper_pair != True:
			newsam_d.write(read)
	newsam_d.close()
	#
	# extract region reads
	extract_bam = bam_full.split('/')[-1][:-4]+'_extract.bam'
	newsam = pysam.AlignmentFile(extract_bam, 'wb', template=samfile) 
	#	
	for item in reg_info:
		chrom = item[0]
		start_pos = int(item[1]) - read_range
		end_pos = int(item[1]) + read_range
		for read in samfile.fetch(chrom, start_pos, end_pos):
			newsam.write(read)
	newsam.close()
	#call(['samtools', 'index' , extract_bam])
#	pysam.index(extract_bam)
	#
	# Sort extracted bam file
	sorted_extracted_bam = extract_bam.split('/')[-1][:-4]+'_sorted.bam'
	pysam.sort('-T', './tmp/aln.exted', '-o', sorted_extracted_bam, extract_bam)
	call(['samtools', 'index' , sorted_extracted_bam])
	#pysam.index(sorted_extracted_bam)	
	#
	# find mate in discordant bam file
	extracted_discord_mate = extract_bam.split('/')[-1][:-4]+'_discord_mate.bam'
	newsam_tm = pysam.AlignmentFile(extracted_discord_mate, 'wb', template=samfile)

	samfile_d = pysam.AlignmentFile(discord_bam, 'rb')
	id_index_d = pysam.IndexedReads(samfile_d)
	id_index_d.build()

	samfile_t = pysam.AlignmentFile(sorted_extracted_bam, 'rb')		
	for read in samfile_t.fetch():
		try:
			iterator = id_index_d.find(read.query_name)
			for tread in iterator:
				if (read.flag & 0x40) != (tread.flag & 0x40):
					newsam_tm.write(tread)
		except KeyError:
			continue

	newsam_tm.close()
	samfile_t.close()
	samfile_d.close()

	extracted_discord_mate_sorted = extracted_discord_mate.split('/')[-1][:-4]+'_sorted.bam'	
	#Sort
	pysam.sort('-T', './tmp/aln.dm', '-o', extracted_discord_mate_sorted, extracted_discord_mate)	
	# Close samfile
	samfile.close()

	# Merge Extracted and discordant mates
	extractd_rds_discord_mate_merged = bam_full.split('/')[-1][:-4]+'_extracted_merged_wdup.bam'
	pysam.merge('-f', extractd_rds_discord_mate_merged, sorted_extracted_bam, extracted_discord_mate_sorted)
	#
	#	print(extractd_rds_discord_mate_merged)
	call(['samtools', 'index' , extractd_rds_discord_mate_merged])
	#	pysam.index(extractd_rds_discord_mate_merged)

	# Remove Duplicates Avoid using picard for this purpose.
	final_extracted_file = bam_full.split('/')[-1][:-4]+'_extracted_final.bam'
	samfile_emwd= pysam.AlignmentFile(extractd_rds_discord_mate_merged, 'rb')
	#
	newsam_emwd = pysam.AlignmentFile(final_extracted_file, 'wb', template=samfile)
	#
	#def read_prop(read):
	#	uniq_read_id = str(read.query_name)+str(read.flag)+str(read.reference_name)+str(read.reference_start)
	#	return uniq_read_id

	read_track = []
	for read in samfile_emwd.fetch():
		uniq_read_id = str(read.query_name)+str(read.flag)+str(read.reference_name)+str(read.reference_start)
		#uniq_read_id = read_prop(read)
		if uniq_read_id not in read_track:
			newsam_emwd.write(read)
			read_track.append(uniq_read_id)

	samfile_emwd.close()
	newsam_emwd.close()

	call(['samtools', 'index' , final_extracted_file])

def exec_remove_reads():
	bam_full = args.bam_inp
	reg_file = args.reg_inp

	with open(reg_file, 'r') as reg_file:
		reg_file_line = reg_file.readlines()
	reg_file.close()

	reg_info = []
	for item in reg_file_line:
		reg_info.append((item.strip().split()[1], item.strip().split()[2]))

	reg_info_sorted = sorted(reg_info, key=lambda x: (x[0], x[1]))
	#print(reg_info_sorted)
	#
	remove_bam = bam_full.split('/')[-1][:-4]+'_remove.bam'
	samfile = pysam.AlignmentFile(bam_full, "rb")
	newsam = pysam.AlignmentFile(remove_bam, "wb", template=samfile)
	#
	chr_track = 'chr00'
	dict_reg = {}
	for item in reg_info_sorted:	
		if item[0] != chr_track:
			chr_track = item[0]
			dict_reg[chr_track] = {}
			for pos in range( (int(item[1])-500), (int(item[1])+501) ):
				dict_reg[chr_track][pos] = 1
			continue
		for pos in range( (int(item[1])-500), (int(item[1])+501) ):	
			dict_reg[chr_track][pos] = 1

	#read.reference_end
	#	aligned reference position of the read on the reference genome.
	#	reference_end points to one past the last aligned residue. 
	#	Returns None if not available (read is unmapped or no cigar alignment present).	

	for read in samfile.fetch():
		try:
			if ( dict_reg[read.reference_name][read.reference_start] == 1 ) or \
				( dict_reg[read.reference_name][read.reference_end] == 1 ):
				continue
		except KeyError:
			newsam.write(read)
	newsam.close()
	samfile.close()

if __name__ == '__main__':
	import sys
	import os
	import pysam
	import argparse
	from subprocess import call

	FUNCTION_MAP = {
			'extract_reads' : exec_extract_reads,
			'remove_reads' : exec_remove_reads
			}

	parser = argparse.ArgumentParser()
	parser.add_argument('module', choices=FUNCTION_MAP.keys())
	subparsers = parser.add_subparsers()

	sp_extract_reads = subparsers.add_parser('extract_reads_opt', help="extract reads argument")
	sp_extract_reads.add_argument('-bam', action='store', dest='bam_inp', help='Bam(.bam) file with full path')
	sp_extract_reads.add_argument('-reg', action='store', dest='reg_inp', help='region file with full path')
	sp_extract_reads.add_argument('-rng', action='store', dest='rng_inp', type=int, default=500, help='Read range')

	sp_remove_reads = subparsers.add_parser('remove_reads_opt', help="remove reads argument")
	sp_remove_reads.add_argument('-bam', action='store', dest='bam_inp', help='Bam(.bam) file with full path')
	sp_remove_reads.add_argument('-reg', action='store', dest='reg_inp', help='region file with full path')

	args = parser.parse_args()
	funct = FUNCTION_MAP[args.module]
	funct()
