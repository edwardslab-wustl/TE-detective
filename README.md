# TEdetective

![Scheme_v1](images/Scheme_v1.png)

*Figure 1: Schematic representation of transposable element insertions.*

```````
```````

![Scheme_v2](images/Scheme_v2.png)

*Figure 2: TE-Detective approach and modules.*

## Requirements
 * [NCBI blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
 * [Censor](https://www.girinst.org/downloads/software/censor/)
 * [BWA](https://github.com/lh3/bwa)
 * Python Packages
     *  numpy
     *  pysam
     *  biopython

## Installation
Clone and install with pip:

````
    git clone https://github.com/edwardslab-wustl/TE-detective.git
    cd TE-detective
    pip install .
````
## Docker Image

   [TE-Detective Docker image](https://hub.docker.com/r/edwardslab/te-detective)
````
	Be sure the following are add to your PATH:
		/opt/conda/bin
		/usr/local/bwa/
		/usr/local/blast-2.2.26/bin
		/usr/local/utils/

	Please set the environmental variable:
		BLASTDIR=/usr/local/blast-2.2.26/bin
````



## Usage

### Input files

````
    1. BAM file
	preferably prepared using following alignment command:
	bwa mem -M -Y -R $RG_LINE ref.fa test_1.fq test_2.fq | samtools view -b -S - > test_ref.bam

    2. ref_fofn file
	Reference file specifying repeat to be examined (e.g. LINE) and location of Reference sequence of repeat elements:
	file can be space- or tab-delimited
	first field is the repeat name
	second field is the full path to a fasta file containing the reference sequences for the repeat element
	Reference sequences of repeat elements can be obtained from Repbase(https://www.girinst.org/server/RepBase/index.php) or other resources.
	See example file ref_fofn in the example_data folder.

	3. bed formatted file of known repeat locations corresponding to the genome version used for alignment.
	You can download the repeatmasker track data from the [UCSC Genome Browser](https://hgdownload.soe.ucsc.edu/downloads.html) and filter with something like:
	zcat rmsk.txt.gz | awk '{print $6"\t"$7"\t"$8"\t"$12;}' > rmsk_hg19.bed

````

### Test example

````
    cd example_data
	TE_detective preprocess -i test_sim.bam -r ref_fofn 
    TE_detective discover -i test_sim.bam -r ref_fofn 
    TE_detective nadiscover -i test_sim.bam -r ref_fofn --polyA --discord_cluster_dens 5  
    TE_detective analyze -i test_sim.bam -r ref_fofn --inp initial_predictions.txt 
    TE_detective cluster2d -i test_sim.bam -r ref_fofn 
    TE_detective filter -i final_results.tsv --bed rmsk_ucsc_mm10.bed

	or 

	cd example_data
	sh run_example.sh
````

### CEU-Trio use case example.

	Setup and assumptions:

		In your working directory ($DIR) you need four subdirectories: one for each sample, and then one for polymorphic analysis (NA12878, NA12891, NA12892, polymorph)
		You have downloaded, aligned and sorted the bam file for each of the three individuals using the instructions above.
		Create a file called ref_fofn in $DIR that has a single line with the repeat element and location (complete path) of a fasta file for the repeat element of interest. See above and/or example ref_fofn file in example_data.

	1. Insertion prediction in child (NA12878).

		cd $DIR/NA12878/
		TE_detective preprocess -i NA12878_hg19_sorted.bam -r ../ref_fofn
		TE_detective discover -i NA12878_hg19_sorted.bam -r ../ref_fofn --read_length 100 --insert_size_est 383
		TE_detective nadiscover -i NA12878_hg19_sorted.bam -r ../ref_fofn -o initial_predictions_NA12878.txt --polyA --read_length 100 --insert_size_est 383
		TE_detective analyze -i NA12878_hg19_sorted.bam  -r ../ref_fofn -o final_results_NA12878.txt --inp initial_predictions_NA12878.txt --read_length 100 --insert_size_est 383

	2. Insertion prediction in one parent (NA12891).

		cd $DIR/NA12891/
		TE_detective preprocess -i NA12891_hg19_sorted.bam -r ../ref_fofn
		TE_detective discover -i NA12891_hg19_sorted.bam -r ../ref_fofn --read_length 100 --insert_size_est 383
		TE_detective nadiscover -i NA12891_hg19_sorted.bam -r ../ref_fofn -o initial_predictions_NA12891.txt --polyA --read_length 100 --insert_size_est 383
		TE_detective analyze -i NA12891_hg19_sorted.bam  -r ../ref_fofn -o final_results_NA12891.txt --inp initial_predictions_NA12891.txt --read_length 100 --insert_size_est 383

	3. Insertion prediction in other parent (NA12892).

		cd $DIR/NA12892/
		TE_detective preprocess -i NA12892_hg19_sorted.bam -r ../ref_fofn
		TE_detective discover -i NA12892_hg19_sorted.bam -r ../ref_fofn --read_length 100 --insert_size_est 383
		TE_detective nadiscover -i NA12892_hg19_sorted.bam -r ../ref_fofn -o initial_predictions_NA12892.txt --polyA --read_length 100 --insert_size_est 383
		TE_detective analyze -i NA12892_hg19_sorted.bam  -r ../ref_fofn -o final_results_NA12892.txt --inp initial_predictions_NA12892.txt --read_length 100 --insert_size_est 383

	4. Polymorphic insertion prediciton in child:
	
		cd $DIR/polymorph	
		TE_detective analyze -o final_results_NA12878_NA12891.txt -i $DIR/NA12891/NA12891_hg19_sorted.bam -r ../ref_fofn --inp $DIR/NA12878/initial_predictions_NA12878.txt -p $DIR/NA12891/preprocessed_files --read_length 100 --insert_size_est 439
		TE_detective analyze -o final_results_NA12878_NA12892.txt -i $DIR/NA12892/NA12892_hg19_sorted.bam -r ../ref_fofn --inp $DIR/NA12878/initial_predictions_NA12878.txt -p $DIR/NA12892/preprocessed_files --read_length 100 --insert_size_est 439

	5. Apply a filter on final_results_NA12878.txt, final_results_NA12878_NA12891.txt and final_results_NA12878_NA12892.txt. After applying filter, a new insertion in child (NA12878) would be those which are found in final_result_NA12878.txt but not in final_result_NA12878_NA12891.txt or final_result_NA12878_NA12892.txt.


## License information
[Censor](http://www.girinst.org/censor/index.php) is distributed under the GPL license.  See details in [Kohany et. al. Bioinformatics 2006](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-7-474).

[NCBI Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) is freely available to the public for use as a "United States Government Work".  See details [here](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/scripts/projects/blast/LICENSE).

### Detailed usage for each command/module

````

1. Preprocess

usage: TE_detective preprocess [-h] -i BAM_INP -r FOFN_REF [-p PREPROCESS_DIR]
                               [--min_clipped_len CLL_INP]
                               [--log_file LOG_FILE]

Processes the input files (indexed BAM file and indexed fasta file), extracts
discordant and clipped reads, as well as creates other files needed in
subsequent steps. Outputs all files to directory specified by --preprocess_dir

required arguments:
  -i BAM_INP, --input_bam BAM_INP
                        input Bam(.bam) file of aligned reads
  -r FOFN_REF, --ref FOFN_REF
                        File with reference sequence paths, see README.md for
                        more info

optional arguments:
  -h, --help            show this help message and exit
  -p PREPROCESS_DIR, --preprocess_dir PREPROCESS_DIR
                        directory to store preprocessing output files
                        (default: preprocessed_files)
  --min_clipped_len CLL_INP
                        Minimum clipped length(bp) (default: 25)
  --log_file LOG_FILE   run log file (default: preprocess.log)


2. Discover

usage: TE_detective discover [-h] -i BAM_INP -r FOFN_REF [-o OUTPUT_FILE]
                             [-p PREPROCESS_DIR] [--insert_size_est ISZ_INP]
                             [--read_length RDL_INP]
                             [--discord_cluster_dens DRD_INP]
                             [--coverage_cutoff CCT_INP]
                             [--min_clipped_len CLL_INP]
                             [--min_map_qual MPQ_INP]
                             [--map_qual_uniq MPQU_INP] [--log_file LOG_FILE]

Uses output from preprocessing step and makes an initial list of candidate
insertions

required arguments:
  -i BAM_INP, --input_bam BAM_INP
                        Input Bam(.bam) file of aligned reads
  -r FOFN_REF, --ref FOFN_REF
                        File with reference sequence paths, see README.md for
                        more info

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Tab-delimited file of initial set of TE insertions
                        (default: initial_predictions.txt)
  -p PREPROCESS_DIR, --preprocess_dir PREPROCESS_DIR
                        directory used to store preprocessing output files
                        (default: preprocessed_files)
  --insert_size_est ISZ_INP
                        Insert size estimate (default: 340)
  --read_length RDL_INP
                        Average read length (default: 150)
  --discord_cluster_dens DRD_INP
                        Discord read cluster density (default: 10)
  --coverage_cutoff CCT_INP
                        Coverage cutoff input (default: 200)
  --min_clipped_len CLL_INP
                        Minimum clipped length(bp) (default: 25)
  --min_map_qual MPQ_INP
                        Minimum mapping quality (default: 30)
  --map_qual_uniq MPQU_INP
                        Minimum mapping quality (default: 1)
  --log_file LOG_FILE   run log file (default: discover.log)


3. Nadiscover

usage: TE_detective nadiscover [-h] -i BAM_INP -r FOFN_REF [--bed RMSK_BED]
                               [-o OUTPUT_FILE] [-p PREPROCESS_DIR]
                               [--min_clipped_len CLL_INP]
                               [--insert_size_est ISZ_INP]
                               [--read_length RDL_INP]
                               [--discord_cluster_dens DRD_INP]
                               [--coverage_cutoff CCT_INP] [--all]
                               [--merge_aligned] [--nonaligned_search]
                               [--min_map_qual MPQ_INP]
                               [--map_qual_uniq MPQU_INP] [--polyA]
                               [--polyA_len PQL_INP]
                               [--polyA_mismatch PMM_INP]
                               [--log_file LOG_FILE]

Performs nonalignment part of the discovery step. Module adds poly A/T
information into predictions made by discovery step. This module performs
initial searches as well, but without using BWA aligner for clipped and
discordant read alignment to TE reference sequence. Instead, a bed file of
masked regions is provided as input, and alignment information from input BAM
file is used.

required arguments:
  -i BAM_INP, --input_bam BAM_INP
                        input Bam(.bam) file of aligned reads
  -r FOFN_REF, --ref FOFN_REF
                        File with reference sequence paths, see README.md for
                        more info

optional arguments:
  -h, --help            show this help message and exit
  --bed RMSK_BED        FoFn for existing repeat elements
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Tab-delimited output file of initial set of TE
                        insertions (default: initial_predictions_noalign.txt)
  -p PREPROCESS_DIR, --preprocess_dir PREPROCESS_DIR
                        directory used to store preprocessing output files
                        (default: preprocessed_files)
  --min_clipped_len CLL_INP
                        Minimum clipped length(bp) (default: 25)
  --insert_size_est ISZ_INP
                        insert Size estimate (default: 340)
  --read_length RDL_INP
                        Average read length (default: 150)
  --discord_cluster_dens DRD_INP
                        discord read cluster density (default: 5)
  --coverage_cutoff CCT_INP
                        Coverage cutoff input (default: 200)
  --all                 use all reads instead of only clipped (default: False)
  --merge_aligned       Merge aligned predictions (default: False)
  --nonaligned_search   Perform non-alignment ref bam search (default: False)
  --min_map_qual MPQ_INP
                        Minimum mapping quality (default: 30)
  --map_qual_uniq MPQU_INP
                        Minimum mapping quality unique test (default: 1)
  --polyA               Perform poly A/T search (default: False)
  --polyA_len PQL_INP   poly A/T Length (default: 9)
  --polyA_mismatch PMM_INP
                        poly A/T mismatch (default: 1)
  --log_file LOG_FILE   run log file (default: nadiscover.log)


4. Analyze (Realignment step from figure 2)

usage: TE_detective analyze [-h] -i BAM_INP -r FOFN_REF --inp LIST_INP
                            [-p PREPROCESS_DIR] [-o OUTPUT_FILE]
                            [--read_length RDL_INP]
                            [--min_clipped_len CLL_INP]
                            [--min_anchor_len AHL_INP]
                            [--clipped_read_range CER_INP]
                            [--clipped_search_interval CSI_INP]
                            [--min_breakpt_reads MRE_INP]
                            [--min_het_reads MRH_INP]
                            [--insert_size_est ISZ_INP]
                            [--mapping_qual_interval QII_INP]
                            [--intervals NII_INP] [--min_map_qual MPQ_INP]
                            [--map_qual_uniq MPQU_INP]
                            [--filter_discord_mates] [--log_file LOG_FILE]

Realigns reads around a predicted insertion point. Can be used to refine
initial predictions from the discover step, or to find evidence of potential
insertions in a different sample (e.g. for polymorphic subtraction). Filter
output from detailed analysis section. User can filter results using the
filter module or manually by importing them into Excel or any other tool

required arguments:
  -i BAM_INP, --input_bam BAM_INP
                        input Bam(.bam) file of aligned reads
  -r FOFN_REF, --ref FOFN_REF
                        File with reference sequence paths, see README.md for
                        more info
  --inp LIST_INP        Input list of insertions

optional arguments:
  -h, --help            show this help message and exit
  -p PREPROCESS_DIR, --preprocess_dir PREPROCESS_DIR
                        directory used to store preprocessing output files
                        (default: preprocessed_files)
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Tab-delimited output file of potential TE
                        insertions(default: final_resutls.tsv)
  --read_length RDL_INP
                        Average read length (default: 150)
  --min_clipped_len CLL_INP
                        Minimum clipped length(bp) (default: 25)
  --min_anchor_len AHL_INP
                        Minimum anchor length(bp) (defualt: 30)
  --clipped_read_range CER_INP
                        Range of clipped reads at a end to put in a group
                        (default: 5)
  --clipped_search_interval CSI_INP
                        Clipped read search interval (default: 20)
  --min_breakpt_reads MRE_INP
                        Min read for breakpoint (default: 4)
  --min_het_reads MRH_INP
                        Minimum reads to call hetrozygous insertion (default:
                        3)
  --insert_size_est ISZ_INP
                        Insert size estimate (default: 340)
  --mapping_qual_interval QII_INP
                        Interval for mapping quality (default: 0.05)
  --intervals NII_INP   Number of intervals (default: 6)
  --min_map_qual MPQ_INP
                        Minimum mapping quality (default: 30)
  --map_qual_uniq MPQU_INP
                        Minimum mapping quality uniq test (default: 1)
  --filter_discord_mates
                        Filter discord mate files (default: False)
  --log_file LOG_FILE   run log file (default: analyze.log)


5. Cluster2D

usage: TE_detective cluster2d [-h] -i BAM_INP -r FOFN_REF [-o OUTPUT_FILE]
                              [-p PREPROCESS_DIR] [--insert_size_est ISZ_INP]
                              [--read_length RDL_INP]
                              [--discord_cluster_dens DRD_INP]
                              [--coverage_cutoff CCT_INP] [--all]
                              [--log_file LOG_FILE]

Optional module to change the discordant read clustering density for initial
prediction without realigning everything. For example, if
--discord_cluster_dens was set to 10 for initial discovery step and user want
to see predictions with --discord_cluster_dens = 5. Uses intermediate files
from discover section and generates new prediction file.

required arguments:
  -i BAM_INP, --input_bam BAM_INP
                        input Bam(.bam) file of aligned reads
  -r FOFN_REF, --ref FOFN_REF
                        File with reference sequence paths, see README.md for
                        more info

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Tab-delimited file of initial set of TE insertions
                        (default: recluster_initial_predictions.txt)
  -p PREPROCESS_DIR, --preprocess_dir PREPROCESS_DIR
                        directory used to store preprocessing output files
                        (default: preprocessed_files)
  --insert_size_est ISZ_INP
                        insert Size estimate (default: 340)
  --read_length RDL_INP
                        Average read length (default: 150)
  --discord_cluster_dens DRD_INP
                        Discord read cluster density (default: 5)
  --coverage_cutoff CCT_INP
                        Coverage cutoff input (default: 200)
  --all                 use all reads instead of only clipped (default: False)
  --log_file LOG_FILE   run log file (default: cluster2d.log)


6. Filter
	
usage: TE_detective filter [-h] -i OFA_INP -b FOFN_BED [-p PREPROCESS_DIR]
                           [--align_qual_lim QLM_INP]
                           [--min_clipped_reads TCR_INP]
                           [--min_clipped_and_dischord_reads TRD_INP]
                           [--read_percent RP_INP] [--read_length RDL_INP]
                           [--insert_size_est ISZ_INP] [--log_file LOG_FILE]

Filter output from analyze step.# Filteration step code looks like this if
total_clipped_rd >= tcr or ( (total_clipped_rd >= mtcr ) and (
(total_clipped_rd_wpat+total_discord_rd) >= trd ) ): filter_result = 'PASS'
elif total_discord_rd >= odrd: filter_result = 'PASS_D' # This flag says
passed based on only discordant reads.

required arguments:
  -i OFA_INP, --input_file OFA_INP
                        use the output file from analyze section as input
  -b FOFN_BED, --bed FOFN_BED
                        File containg a list of files to existing repeat
                        elements. List the full path for each file. See
                        example in example_data

optional arguments:
  -h, --help            show this help message and exit
  -p PREPROCESS_DIR, --preprocess_dir PREPROCESS_DIR
                        directory used to store preprocessing output files
                        (default: preprocessed_files)
  --align_qual_lim QLM_INP
                        Lowest limit for alignment quality (default: 0.85)
  --min_clipped_reads TCR_INP
                        Minimum number of clipped reads (default: 5)
  --min_clipped_and_dischord_reads TRD_INP
                        Minimum total [clipped+discordant] reads (default: 10)
  --read_percent RP_INP
                        read percent value (default: 10.0)
  --read_length RDL_INP
                        Average read length (default: 150)
  --insert_size_est ISZ_INP
                        insert Size estimate (default: 340)
  --log_file LOG_FILE   run log file (default: filter.log)

````
