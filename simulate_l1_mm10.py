#Type, Length and position of insertion
class te_seq:

	def __init__(self, te_class, cnt):
		self.te_class = te_class
		self.cnt = cnt
		self.utr3_len = 600
		self.genome_seq = '/scratch/mksingh/ecat11/l1_base/Mus_musculus.GRCm38.dna.toplevel.fa'
		with open('/scratch/mksingh/ecat11/l1_base/fli-l1_fltrd.csv', 'r') as l1base_fl:
		        self.l1base_fl_lines = l1base_fl.readlines()
		l1base_fl.close()

	def var_len(self):
		with open('/scratch/mksingh/ecat11/rmsk_ucsc/rmsk_ucsc_mm10.bed', 'r') as bed_fl:
			bed_fl_lines = bed_fl.readlines()
		bed_fl.close

		len_arr = []
		ref_chr = 'chr00'
		prev_start = 0
		prev_end = 0
		for line in bed_fl_lines:
			if line.split()[3] == 'LINE':
				current_start = int( line.split()[1] )
				current_end = int( line.split()[2] )
				if current_start > prev_end+1 or line.split()[0] != ref_chr:
					l1_length = prev_end - prev_start
					len_arr.append( l1_length )
					prev_start = current_start
					prev_end = current_end
					ref_chr = line.split()[0]
					continue
				prev_end = current_end	
		len_value = np.random.choice(len_arr)
		del bed_fl_lines[:]
		return(len_value)

	def find_seq(self):

		self.te_len = self.var_len()
		l1_info = self.l1base_fl_lines[ np.random.random_integers( 0, high=( len( self.l1base_fl_lines ) )-1 ) ]

		chrom = l1_info.split(',')[2][1:-1]
		loc_start = int(l1_info.split(',')[3][1:-1])
		loc_end = int(l1_info.split(',')[4][1:-1])
		strand = int(l1_info.split(',')[7][1:-1])
		orf2_end = int(l1_info.split(',')[18][1:-1])

		if strand == 1: #
			utr3_end = orf2_end + self.utr3_len
			start_idx = loc_start + utr3_end - self.te_len
			end_idx = loc_start + utr3_end
		if strand == -1: #
			utr3_end = orf2_end + self.utr3_len
			end_idx = loc_start - utr3_end + self.te_len # compared to start
			start_idx = loc_start - utr3_end

		self.reg = chrom+':'+str(start_idx)+'-'+str(end_idx)
		te_fasta = pysam.faidx(self.genome_seq, self.reg)
		seq_str = ''.join( te_fasta.split('\n')[1:] )
		self.seq = textwrap.fill( seq_str, width=70 )

		self.output_file_name = self.te_class+'_seq_for_insert_'+self.cnt+'.txt'
		with open(self.output_file_name , 'w') as output_file:
			output_file.write(self.seq)
		output_file.close()

		sys.stdout.write(self.reg+'\t'+self.output_file_name+'\t'+str(self.te_len)+'\n')
		return(self.reg, self.output_file_name)

class find_ref_seq:

	def __init__(self, ref_seq_name, te_class):
		self.ref_seq_name = ref_seq_name
		self.te_class = te_class
		self.chrom_size = {
				'chr1'	:	195471971,
				'chr2'	:	182113224,
				'chr3'	:	160039680,
				'chr4'	:	156508116,
				'chr5'	:	151834684,
				'chr6'	:	149736546,
				'chr7'	:	145441459,
				'chr8'	:	129401213,
				'chr9'	:	124595110,
				'chr10'	:	130694993,
				'chr11'	:	122082543,
				'chr12'	:	120129022,
				'chr13'	:	120421639,
				'chr14'	:	124902244,
				'chr15'	:	104043685,
				'chr16'	:	98207768,
				'chr17'	:	94987271,
				'chr18'	:	90702639,
				'chr19'	:	61431566,
				'chrX'	:	171031299,
				'chrY'	:	91744698,
			}
		
#	def list_chr(self):
		self.chrm_list = []
		with open(self.ref_seq_name, 'r') as seq_file:
			for line in seq_file:
				if line.startswith('>'):
					self.chrm_list.append(line[1:-1])
		find_chr = True
		while (find_chr==True):
			self.chrom = self.chrm_list[np.random.random_integers(0, high=(len(self.chrm_list)-1))]
			if self.chrom not in ['chrM', 'chrX', 'chrY']:
				find_chr = False
		sys.stdout.write(self.chrom+'\n')
		self.seqinfo = []
		self.non_te_reg = []
#		return(self.chrom)

	def prop_bed(self):
		find_ip = True
		while (find_ip == True):
			init_ip = np.random.random_integers(1, self.chrom_size[self.chrom])
			init_reg = self.chrom+':'+str(init_ip-1000)+'-'+str(init_ip+1000)
			init_reg_seq = pysam.faidx(self.ref_seq_name, init_reg)
			if init_reg_seq.lower().count('n') == 0:
				find_ip = False

		self.ip_te = init_ip
		return(self.ip_te)
	
	def seq_info(self):
		del self.seqinfo[:]
		self.chrseq = pysam.faidx(self.ref_seq_name, self.chrom)
	#	print(len(self.chrseq.split('\n')))
		poscnt = 0
		for i, line in enumerate(self.chrseq.split('\n')):
			if line.startswith('>') != True:
	#		if (len(line)-line.lower().count('n')) != 0 and line.startswith('>') != True:
				poscnt = poscnt + (len(line)) #-line.lower().count('n'))
				self.seqinfo.append([i+1, poscnt])
		return(self.seqinfo)

	def seq_line(self, chrom_line_inp):
		self.chrom_line = chrom_line_inp
		self.count_chrom_line = 'n'
		self.chrom_line_cnt = -1
		for i, line in enumerate(open(self.ref_seq_name, 'r')):
			if line.startswith('>') and line[1:-1] == self.chrom:
				self.count_chrom_line = 'y'
				self.chrom_line_cnt = 0
			if self.count_chrom_line == 'y':
				self.chrom_line_cnt += 1
			if self.chrom_line_cnt == self.chrom_line:
				break
		self.seq_line_chrom = i+1
		return(self.seq_line_chrom)	
	
	def find_ip(self):
		self.initial_insert_point = self.prop_bed() 
		sys.stdout.write(str(self.initial_insert_point)+'\n')

		for info in self.seq_info():
			if info[1] > self.initial_insert_point:
				sys.stdout.write(str(info[0])+' '+str(info[1])+'\n')
				self.chrom_line_insert = info[0]
				self.chrom_pos_insert = info[1]
				break
		self.tsd_len = np.random.choice([0, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20])
		tsd_reg = self.chrom+':'+str(self.chrom_pos_insert-self.tsd_len+1)+'-'+str(self.chrom_pos_insert)
		tsd_seq_init = pysam.faidx(self.ref_seq_name, tsd_reg)
		tsd_seq = ''.join( tsd_seq_init.split('\n')[1:] )

		self.tsd_out_fl = self.chrom+'_'+str(self.chrom_line_insert)+'.txt'
		with open(self.tsd_out_fl , 'w') as output_file:
			output_file.write(tsd_seq)
		output_file.close()

		self.seq_line_insert = self.seq_line(self.chrom_line_insert)
		return(self.seq_line_insert, self.te_class, self.chrom, self.chrom_line_insert, self.chrom_pos_insert, self.tsd_out_fl, self.tsd_len)	

class ContinueI(Exception):
	pass

def count_bases(input_file_name):
	total_base_cnt = 0
	with open(input_file_name, 'r') as inp_seq_file:
		for line in inp_seq_file:
			total_base_cnt = total_base_cnt + (len(line)-line.count('\n'))
	inp_seq_file.close()
	return(total_base_cnt) 

def setup_tes_for_writing(te_class_inp, num_of_te, min_sepetation_inp, name_inserts_seq, te_track_inp, inserts_pos_inp, ref_seq_file_inp): 

	for i in range(0, num_of_te):
		i_str = str(i)
		te_type, te_flnm = te_seq(te_class_inp, i_str).find_seq()
		name_inserts_seq.append([te_type, te_flnm])
	i = 0
	while(i < num_of_te):
		insert_estimate = find_ref_seq(ref_seq_file_inp, te_class_inp).find_ip()
		try:
			for itm in te_track_inp:
				if insert_estimate[1] == itm[1] and insert_estimate[4] in range(itm[4]-min_sepetation_inp, itm[4]+min_sepetation_inp):
					raise ContinueI()
			inserts_pos_inp.append(insert_estimate)
			te_track_inp.append(insert_estimate)
			i += 1
		except ContinueI:
#			sys.stdout.write(str(insert_estimate[4])+'\texception raised\n')
			continue


if __name__ == '__main__':
	import sys
	import os
	import pysam
	import numpy as np
	import textwrap

	iden = sys.argv[1]
	ref_seq_file = '../mm10.fa'
	ref_seq_te_file = 'mm10_withTE_'+iden+'.fa'
	insertion_summary_file = 'summary_of_insertions_'+iden+'.txt'

	num_line = 500
	min_sepetation = 1000
	te_track = []
	inserts_pos = []

	if num_line > 0:
		line_inserts_seq = []
		setup_tes_for_writing('LINE', num_line, min_sepetation, line_inserts_seq, te_track, inserts_pos, ref_seq_file)

	inserts_pos_sorted = sorted(inserts_pos, key=lambda x : x[0])	

	te_inserted_output_file = open(ref_seq_te_file, 'w') 
	summary_file = open(insertion_summary_file, 'w')

	line_start = 0
	line_seq_idx = 0
	sine_seq_idx = 0
	ltr_seq_idx = 0

	for pos in inserts_pos_sorted:
		seq_file_tmp = open(ref_seq_file, 'r')
		for i, line_sq in enumerate(seq_file_tmp):	
			if i in range(line_start, pos[0]):
				te_inserted_output_file.write(line_sq)
		seq_file_tmp.close()
		line_start = pos[0]
		
		summary_file.write(str(pos[1])+' '+str(pos[2])+' '+str(pos[4])+' '+str(pos[6])+' ')

		with open(line_inserts_seq[line_seq_idx][1], 'r') as line_fl:
			line_fl_ln = line_fl.readlines()
		line_fl.close()
		l1_seq_str = ''.join( itm.strip() for itm in line_fl_ln )
		with open(pos[5], 'r') as tsd_fl:				
			tsd_ln = tsd_fl.readlines()
		tsd_fl.close()
		tsd_ln_str = ''.join( itm.strip() for itm in tsd_ln )
#		print(tsd_ln_str)

		l1_tsd_seq = l1_seq_str+tsd_ln_str
#		print(l1_tsd_seq)

		l1_tsd_seq_fa = textwrap.fill( l1_tsd_seq, width=60 )
#		l1_tsd_seq_fa = l1_tsd_seq_fa+'\n'
#		print(l1_tsd_seq_fa)

		te_inserted_output_file.write(l1_tsd_seq_fa)
		summary_file.write(line_inserts_seq[line_seq_idx][0]+' '+str(count_bases(line_inserts_seq[line_seq_idx][1]))+'\n' )
		line_seq_idx += 1

	seq_file_tmp = open(ref_seq_file, 'r')
	for i, line_sq in enumerate(seq_file_tmp):
		if i >= line_start:
			te_inserted_output_file.write(line_sq)
	seq_file_tmp.close()

	summary_file.close()
	te_inserted_output_file.close()

