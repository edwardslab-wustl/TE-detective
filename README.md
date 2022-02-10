# TEdetective

![Scheme_v1](images/Scheme_v1.png)

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

## Usage

### Input files

### Example

````
    cd example_data
    TE_detective preprocess -bam test_sim.bam -ref ref_fofn
    TE_detective discover -bam test_sim.bam -ref ref_fofn
    TE_detective analyze -bam test_sim.bam -ref ref_fofn -inp initial_predictions.txt
    TE_detective nadiscover -bam test_sim.bam -ref ref_fofn
    TE_detective cluster2d -bam test_sim.bam -ref ref_fofn
    TE_detective filter -ofa final_results -bed rmsk_ucsc_mm10.bed
````

## License information
[Censor](http://www.girinst.org/censor/index.php) is distributed under the GPL license.  See details in [Kohany et. al. Bioinformatics 2006](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-7-474).

[NCBI Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) is freely available to the public for use as a "United States Government Work".  See details [here](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/scripts/projects/blast/LICENSE).

### Parameters

````

1. Preprocess:

	This module processes the input files and create input files for next steps. 

	-bam  : Input indexed bam file (aligned with bwa -mem).
	-ref  : File of file name of TE reference fasta file. Please provide file name with absolute path.
	-cll  : Length of clipped reads to be extracted. 

2. Discover:
	
	This module makes initial list of candidate insertions. 	

	-bam  : Input indexed bam file (aligned with bwa -mem).
	-ref  : File of file name of TE reference fasta file. Please provide file name with absolute path.
	-isz  : Insert size estimate. ( = mean_insert_size + 2 * insert_size_standard_deviation – read_length)
	-rdl  : Average read length of bam file (can be estimated using picard).
	-drd  : Number of supporting reads (discordant + clipped) for calling an insertion.
	-cct  : A region with coverage more than this will be ignored from initial prediction. 
	-ccl  : Minimum length of clipped read to be analyzed.
	-mpq  : Minimum mapping quality of a read.
	-mpqu : Value of a mapping quality which is used by uniqness testing algorithm (used for clipped reads).

3. Analyze:

	Detailed analysis of initial list of candidate insertions. 	

	-bam  : Input indexed bam file (aligned with bwa -mem).
	-ref  : File of file name of TE reference fasta file. Please provide file name with absolute path.
	-inp  : Initial prediction file from discover step.
	-rdl  : Average read length of bam file.
	-cll  : Minimum length of clipped read to be analyzed.
	-ahl  : Minimum anchor length of clipped reads.
	-cer  : Range of clipped reads ends to be put in a group.
	-csi  : Clipped read search interval ( from initial prediction )
	-mre  : Minimum reads to support a breakpoint (while determining exact insertion point.)
	-mrh  : Minimum supporting reads to call hetrozygous insertion.
	-isz  : Insert size estimate. ( = mean_insert_size + 2 * insert_size_standard_deviation – read_length)
	-qii  : Interval for mapping quality of reads from Censor output.
	-nii  : Number of intervals for mapping quality of reads from Censor output to be searched and printed.
	-mpq  : Minimum mapping quality of a read.
	-mpqu : Value of a mapping quality which is used by uniqness testing algorithm (used for clipped reads).

4. Nadiscover:

	This module makes initial list of candidate insertions without use of bwa aligner. Instead, a bed file of masked regions are provided as input. 

	-bam  : Input indexed bam file (aligned with bwa -mem).
	-ref  : File of file name of TE reference fasta file. Please provide file name with absolute path.
	-cll  : Minimum length of clipped read to be analyzed.
	-isz  : Insert size estimate. ( = mean_insert_size + 2 * insert_size_standard_deviation – read_length)
	-rdl  : Average read length of bam file.
	-drd  : Number of supporting reads (discordant + clipped) for calling an insertion. 
	-cct  : A region with coverage more than this will be ignored from prediction.
	-all  : If set false, only clipped reads will be taken into consideration. 
	-mrg  : Set true if you want to merge this part of analysis with alignment module of initial prediction.
	-pat  : Include P/T analysis in prediction.  
	-nas  : Non-alignment ref bam search
	-mpq  : Minimum mapping quality of a read.
	-pql  : Poly A/T Length
	-pmm  : Maximum poly A/T mismatch
	-mpqu : Value of a mapping quality which is used by uniqness testing algorithm (used for clipped reads).
	-bed  : BED file of existing repeat elements ( START	END	TE_CLASS ) 

5. Cluster2D:

	Use this module if you want to change the clustering desnity for initial prediction. This module uses intermediate files from discover section and generates new prediction file.

	-bam  : Input indexed bam file (aligned with bwa -mem).
	-ref  : File of file name of TE reference fasta file. Please provide file name with absolute path.
	-isz  : Insert size estimate. ( = mean_insert_size + 2 * insert_size_standard_deviation – read_length)
	-rdl  : Average read length of bam file.
	-drd  : Number of supporting reads (discordant + clipped) for calling an insertion.
	-cct  : A region with coverage more than this will be ignored from prediction.
	-all  : If set false, only clipped reads will be taken into consideration.

6. Filter:
	
	Filter output from detailed analysis section.

	-ofa  : Output file from analyze section
	-bed  : BED file of existing repeat elements ( START    END     TE_CLASS )
	-qlm  : Lowest limit for Censor alignment quality
	-tcr  : Minimum number of clipped reads
	-trd  : Minimum total [clipped+discordant] reads
	-ref  : File of file name of TE reference fasta file. Please provide file name with absolute path.
	-rp   : Supporting read percernt as a fraction of total read.
	-rdl  : Average read length of bam file.
	-isz  : Insert size estimate. ( = mean_insert_size + 2 * insert_size_standard_deviation – read_length)

````
