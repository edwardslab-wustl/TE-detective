# TEdetective

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

### Example

''''
    cd example_data
    TE_detective preprocess -bam test_sim.bam -ref ref_fofn
    TE_detective discover -bam test_sim.bam -ref ref_fofn
    TE_detective analyze -bam test_sim.bam -ref ref_fofn -inp initial_predictions.txt
    TE_detective nadiscover -bam test_sim.bam -ref ref_fofn
    TE_detective cluster2d -bam test_sim.bam -ref ref_fofn
''''

## License information
[Censor](http://www.girinst.org/censor/index.php) is distributed under the GPL license.  See details in [Kohany et. al. Bioinformatics 2006](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-7-474).
