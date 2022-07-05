################################################################
# Python based tool to identify novel transposable element insertions from WGS data.
# Contributors: Manoj Kumar Singh (manoj@wustl.edu),
#        and John Edwards (jredwards@wustl.edu)
#
# External dependencies: censor, NCBI blast (provided with package)
#
################################################################

import argparse
from TEdetective.preprocess import exec_preprocess
from TEdetective.discover import exec_discover
from TEdetective.nadiscover import exec_nadiscover
from TEdetective.cluster2d import exec_cluster2d
from TEdetective.filters import exec_filter, exec_filter_p, exec_filter_p_ceu
from TEdetective.analyze import exec_analyze
from TEdetective.polymorph_screen import exec_polymorph, polymorph_setup_arg_parser

def main():
    FUNCTION_MAP = {
            'preprocess' : exec_preprocess, 
            'discover' : exec_discover,
            'nadiscover' : exec_nadiscover, 
            'analyze' : exec_analyze,
            'cluster2d' : exec_cluster2d,
            'filter' : exec_filter,
            'filter_p' : exec_filter_p,
            'filter_p_ceu' : exec_filter_p_ceu,
            'polymorph' : exec_polymorph
            }

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='command')
    subparsers.required = True
    
    polymorph_desc=str("polymorphic subtraction module")
    sp_polymorph = subparsers.add_parser('polymorph',
        help=polymorph_desc, description=polymorph_desc)
    sp_polymorph = polymorph_setup_arg_parser(sp_polymorph)
    
    preprocess_desc=str("Processes the input files (indexed BAM file and indexed " +
        "fasta file), extracts discordant and clipped reads, as well as  " +
        "creates other files needed in subsequent steps. Outputs all files " +
        "to directory specified by --preprocess_dir")
    sp_preprocess = subparsers.add_parser('preprocess',
        description=preprocess_desc, help=preprocess_desc)
    sp_preprocess_required = sp_preprocess.add_argument_group('required arguments')
    sp_preprocess_required.add_argument('-i', '--input_bam', action='store', dest='bam_inp', required=True,
        help='input Bam(.bam) file of aligned reads')
    sp_preprocess_required.add_argument('-r', '--ref', action='store', dest='fofn_ref', required=True,
        help='File with reference sequence paths, see README.md for more info')
    sp_preprocess.add_argument('-p', '--preprocess_dir', action='store',
        dest='preprocess_dir', default='preprocessed_files',
        help='directory to store preprocessing output files (default: preprocessed_files)')
    sp_preprocess.add_argument('--min_clipped_len', action='store', dest='cll_inp', type=int, default=25,
        help='Minimum clipped length(bp) (default: 25)')
    sp_preprocess.add_argument('--log_file', action='store',
        dest='log_file', default='preprocess.log',
        help='run log file (default: preprocess.log)')
    sp_preprocess._action_groups.reverse()

    discover_desc=str("Uses output from preprocessing step and makes an initial list of candidate insertions")
    sp_discover = subparsers.add_parser('discover',
        help=discover_desc, description=discover_desc)
    sp_discover_required = sp_discover.add_argument_group('required arguments')
    sp_discover_required.add_argument('-i', '--input_bam', action='store', dest='bam_inp', required=True, 
        help='Input Bam(.bam) file of aligned reads')
    sp_discover_required.add_argument('-r', '--ref', action='store', dest='fofn_ref', required=True,
        help='File with reference sequence paths, see README.md for more info')
    sp_discover.add_argument('-o', '--output_file', action='store',
        dest='output_file', default='initial_predictions.txt',
        help='Tab-delimited file of initial set of TE insertions (default: initial_predictions.txt)')
    sp_discover.add_argument('-p', '--preprocess_dir', action='store',
        dest='preprocess_dir', default='preprocessed_files',
        help='directory used to store preprocessing output files (default: preprocessed_files)')
    sp_discover.add_argument('--insert_size_est', action='store', dest='isz_inp', type=int, default=340, 
        help='Insert size estimate (default: 340)')
    sp_discover.add_argument('--read_length', action='store', dest='rdl_inp', type=int, default=150, 
        help='Average read length (default: 150)')
    sp_discover.add_argument('--discord_cluster_dens', action='store', dest='drd_inp', type=int, default=10, 
        help='Discord read cluster density (default: 10)')
    sp_discover.add_argument('--coverage_cutoff', action='store', dest='cct_inp', type=int, default=200, 
        help='Coverage cutoff input (default: 200)')
    sp_discover.add_argument('--min_clipped_len', action='store', dest='cll_inp', type=int, default=25,
        help='Minimum clipped length(bp) (default: 25)')
    sp_discover.add_argument('--min_map_qual', action='store', dest='mpq_inp', type=int, default=30, 
        help='Minimum mapping quality (default: 30)')
    sp_discover.add_argument('--map_qual_uniq', action='store', dest='mpqu_inp', type=int, default=1, 
        help='Minimum mapping quality (default: 1)')
    sp_discover.add_argument('--log_file', action='store',
        dest='log_file', default='discover.log',
        help='run log file (default: discover.log)')
    sp_discover._action_groups.reverse()

    analyze_desc=str("Realigns reads around a predicted insertion point. " +
          "Can be used to refine initial predictions from the " +
          "discover step, or to find evidence of potential " +
          "insertions in a different sample (e.g. for polymorphic " +
          "subtraction). Filter output from detailed analysis " +
          "section. User can filter results using the filter " +
          "module or manually by importing them into Excel or " +
          "any other tool")
    sp_analyze = subparsers.add_parser('analyze', 
        description=analyze_desc, help=analyze_desc)
    sp_analyze_required = sp_analyze.add_argument_group('required arguments')
    sp_analyze_required.add_argument('-i','--input_bam', action='store', dest='bam_inp', required=True, 
        help='input Bam(.bam) file of aligned reads')
    sp_analyze_required.add_argument('-r', '--ref', action='store', dest='fofn_ref', required=True,
        help='File with reference sequence paths, see README.md for more info')
    sp_analyze_required.add_argument('--inp', action='store', dest='list_inp', required=True, 
        help='Input list of insertions')
    sp_analyze.add_argument('-p', '--preprocess_dir', action='store',
        dest='preprocess_dir', default='preprocessed_files',
        help='directory used to store preprocessing output files (default: preprocessed_files)')
    sp_analyze.add_argument('-o', '--output_file', action='store',
        dest='output_file', default='final_results.tsv',
        help='Tab-delimited output file of potential TE insertions(default: final_resutls.tsv)')
    sp_analyze.add_argument('--read_length', action='store', dest='rdl_inp', type=int, default=150, 
        help='Average read length (default: 150)')
    sp_analyze.add_argument('--min_clipped_len', action='store', dest='cll_inp', type=int, default=25, 
        help='Minimum clipped length(bp) (default: 25)')
    sp_analyze.add_argument('--min_anchor_len', action='store', dest='ahl_inp', type=int, default=30, 
        help='Minimum anchor length(bp) (defualt: 30)')
    sp_analyze.add_argument('--clipped_read_range', action='store', dest='cer_inp', type=int, default=5,
        help='Range of clipped reads at a end to put in a group (default: 5)')
    sp_analyze.add_argument('--clipped_search_interval', action='store', dest='csi_inp', type=int, default=20, 
        help='Clipped read search interval (default: 20)')
    sp_analyze.add_argument('--min_breakpt_reads', action='store', dest='mre_inp', type=int, default=4, 
        help='Min read for breakpoint (default: 4)')
    sp_analyze.add_argument('--min_het_reads', action='store', dest='mrh_inp', type=int, default=3, 
        help='Minimum reads to call hetrozygous insertion (default: 3)')
    sp_analyze.add_argument('--insert_size_est', action='store', dest='isz_inp', type=int, default=340,
        help='Insert size estimate (default: 340)')
    sp_analyze.add_argument('--mapping_qual_interval', action='store', dest='qii_inp', type=float, default=0.05, 
        help='Interval for mapping quality (default: 0.05)')
    sp_analyze.add_argument('--intervals', action='store', dest='nii_inp', type=int, default=6,
        help='Number of intervals (default: 6)')
    sp_analyze.add_argument('--min_map_qual', action='store', dest='mpq_inp', type=int, default=30, 
        help='Minimum mapping quality (default: 30)')
    sp_analyze.add_argument('--map_qual_uniq', action='store', dest='mpqu_inp', type=int, default=1,
        help='Minimum mapping quality uniq test (default: 1)')
    sp_analyze.add_argument('--filter_discord_mates', action='store_true', dest='flt_inp', default=False,
        help='Filter discord mate files (default: False)')
    sp_analyze.add_argument('--log_file', action='store',
        dest='log_file', default='analyze.log',
        help='run log file (default: analyze.log)')
    sp_analyze._action_groups.reverse()

    nadiscover_desc=str("Performs nonalignment part of the discovery step. " +
          "Module adds poly A/T information into predictions made " +
          "by discovery step. This module performs initial searches " +
          "as well, but without using BWA aligner for clipped and " +
          "discordant read alignment to TE reference sequence. " +
          "Instead, a bed file of masked regions is provided as " +
          "input, and alignment information from input BAM file is used.")
    sp_nadiscover = subparsers.add_parser('nadiscover', 
        description=nadiscover_desc,help=nadiscover_desc)
    sp_nadiscover_required = sp_nadiscover.add_argument_group('required arguments')
    sp_nadiscover_required.add_argument('-i','--input_bam', action='store', dest='bam_inp', required=True, 
        help='input Bam(.bam) file of aligned reads')
    sp_nadiscover_required.add_argument('-r', '--ref', action='store', dest='fofn_ref', required=True,
        help='File with reference sequence paths, see README.md for more info')
    sp_nadiscover.add_argument('--bed', action='store', dest='rmsk_bed', help='FoFn for existing repeat elements')
    sp_nadiscover.add_argument('-o', '--output_file', action='store',
        dest='output_file', default='initial_predictions_noalign.txt',
        help='Tab-delimited output file of initial set of TE insertions (default: initial_predictions_noalign.txt)')
    sp_nadiscover.add_argument('-p', '--preprocess_dir', action='store',
        dest='preprocess_dir', default='preprocessed_files',
        help='directory used to store preprocessing output files (default: preprocessed_files)')
    sp_nadiscover.add_argument('--min_clipped_len', action='store', dest='cll_inp', type=int, default=25,
        help='Minimum clipped length(bp) (default: 25)')
    sp_nadiscover.add_argument('--insert_size_est', action='store', dest='isz_inp', type=int, default=340,
        help='insert Size estimate (default: 340)')
    sp_nadiscover.add_argument('--read_length', action='store', dest='rdl_inp', type=int, default=150,
        help='Average read length (default: 150)')
    sp_nadiscover.add_argument('--discord_cluster_dens', action='store', dest='drd_inp', type=int, default=5,
        help='discord read cluster density (default: 5)')
    sp_nadiscover.add_argument('--coverage_cutoff', action='store', dest='cct_inp', type=int, default=200,
        help='Coverage cutoff input (default: 200)')
    sp_nadiscover.add_argument('--all', action='store_true', dest='flg_all', default=False,
        help='use all reads instead of only clipped (default: False)')
    sp_nadiscover.add_argument('--merge_aligned', action='store_true', dest='merged', default=False,
        help='Merge aligned predictions (default: False)')
    sp_nadiscover.add_argument('--nonaligned_search', action='store_true', dest='nas_inp', default=False,
        help='Perform non-alignment ref bam search (default: False)')
    sp_nadiscover.add_argument('--min_map_qual', action='store', dest='mpq_inp', type=int, default=30,
        help='Minimum mapping quality (default: 30)')
    sp_nadiscover.add_argument('--map_qual_uniq', action='store', dest='mpqu_inp', type=int, default=1,
        help='Minimum mapping quality unique test (default: 1)')
    sp_nadiscover.add_argument('--polyA', action='store_true', dest='pat_inp', default=False,
        help='Perform poly A/T search (default: False)')
    sp_nadiscover.add_argument('--polyA_len', action='store', dest='pql_inp', type=int, default=9,
        help='poly A/T Length (default: 9)')
    sp_nadiscover.add_argument('--polyA_mismatch', action='store', dest='pmm_inp', type=int, default=1,
        help='poly A/T mismatch (default: 1)')
    sp_nadiscover.add_argument('--log_file', action='store',
        dest='log_file', default='nadiscover.log',
        help='run log file (default: nadiscover.log)')
    sp_nadiscover._action_groups.reverse()

    cluster2d_desc=str("Optional module to change the discordant read clustering density for " +
          "initial prediction without realigning everything. For example, if " +
          "--discord_cluster_dens was set to 10 for initial discovery step " +
          "and user want to see predictions with --discord_cluster_dens = 5. " +
          "Uses intermediate files from discover section and generates new prediction file.")
    sp_cluster2d = subparsers.add_parser('cluster2d', 
        help=cluster2d_desc, description=cluster2d_desc)
    sp_cluster2d_required = sp_cluster2d.add_argument_group('required arguments')
    sp_cluster2d_required.add_argument('-i','--input_bam', action='store', dest='bam_inp', required=True, 
        help='input Bam(.bam) file of aligned reads')
    sp_cluster2d_required.add_argument('-r', '--ref', action='store', dest='fofn_ref', required=True,
        help='File with reference sequence paths, see README.md for more info')
    sp_cluster2d.add_argument('-o', '--output_file', action='store',
        dest='output_file', default='recluster_initial_predictions.txt',
        help='Tab-delimited file of initial set of TE insertions (default: recluster_initial_predictions.txt)')
    sp_cluster2d.add_argument('-p', '--preprocess_dir', action='store',
        dest='preprocess_dir', default='preprocessed_files',
        help='directory used to store preprocessing output files (default: preprocessed_files)')
    sp_cluster2d.add_argument('--insert_size_est', action='store', dest='isz_inp', type=int, default=340,
        help='insert Size estimate (default: 340)')
    sp_cluster2d.add_argument('--read_length', action='store', dest='rdl_inp', type=int, default=150, 
        help='Average read length (default: 150)')
    sp_cluster2d.add_argument('--discord_cluster_dens', action='store', dest='drd_inp', type=int, default=5, 
        help='Discord read cluster density (default: 5)')
    sp_cluster2d.add_argument('--coverage_cutoff', action='store', dest='cct_inp', type=int, default=200, 
        help='Coverage cutoff input (default: 200)')
    sp_cluster2d.add_argument('--all', action='store_true', dest='flg_all', default=False, 
        help='use all reads instead of only clipped (default: False)')
    sp_cluster2d.add_argument('--log_file', action='store',
        dest='log_file', default='cluster2d.log',
        help='run log file (default: cluster2d.log)')
    sp_cluster2d._action_groups.reverse()

    filter_desc=str("Filter output from analyze step." +
            "# Filteration step code looks like this\n" +
            "if total_clipped_rd >= tcr or ( (total_clipped_rd >= mtcr ) " +
            "and ( (total_clipped_rd_wpat+total_discord_rd) >= trd ) ):\n" +
            "filter_result = 'PASS'\n" +
            "elif total_discord_rd >= odrd: \n" +
            "filter_result = 'PASS_D' # This flag says passed based on only discordant reads.")
    sp_filter = subparsers.add_parser('filter', 
        help=filter_desc, description=filter_desc)
    sp_filter_required = sp_filter.add_argument_group('required arguments')
    sp_filter_required.add_argument('-i', '--input_file', action='store', dest='ofa_inp', required=True,
        help='use the output file from analyze section as input')
    sp_filter_required.add_argument('-b', '--bed', action='store', dest='fofn_bed', required=True, 
        help='File containg a list of files to existing repeat elements. ' +
        'List the full path for each file. See example in example_data')
    #sp_filter.add_argument('-p', '--preprocess_dir', action='store',
    #    dest='preprocess_dir', default='preprocessed_files',
    #    help='directory used to store preprocessing output files (default: preprocessed_files)')
    sp_filter.add_argument('--align_qual_lim', action='store', dest='qlm_inp', type=float, default=0.85, 
        help='Lowest limit for alignment quality (default: 0.85)')
    sp_filter.add_argument('--min_clipped_reads', action='store', dest='tcr_inp', type=int, default=5, 
        help='Minimum number of clipped reads (default: 5)')
    sp_filter.add_argument('--min_clipped_and_dischord_reads', action='store', dest='trd_inp', type=int, default=10, 
        help='Minimum total [clipped+discordant] reads (default: 10)')
    sp_filter.add_argument('--read_percent', action='store', dest='rp_inp', type=float, default=10.0, 
        help='read percent value (default: 10.0)')
#    sp_filter.add_argument('--read_length', action='store', dest='rdl_inp', type=int, default=150, 
#        help='Average read length (default: 150)')
#    sp_filter.add_argument('--insert_size_est', action='store', dest='isz_inp', type=int, default=340, 
#        help='insert Size estimate (default: 340)')
    sp_filter.add_argument('--log_file', action='store',
        dest='log_file', default='filter.log',
        help='run log file (default: filter.log)')
    sp_filter._action_groups.reverse()

    filter_p_desc=str("alternative filter for output from analyze step")
    sp_filter_p = subparsers.add_parser('filter_p', 
        help=filter_p_desc, description=filter_p_desc)
    sp_filter_p_required = sp_filter_p.add_argument_group('required arguments')
    sp_filter_p_required.add_argument('-i', '--input_file', action='store', dest='ofa_inp', required=True,
        help='use the output file from analyze section as input')
    sp_filter_p_required.add_argument('-b', '--bed', action='store', dest='fofn_bed', required=True, 
        help='File containg a list of files to existing repeat elements. List the full path for each file. See example in example_dir')
    #sp_filter_p.add_argument('-p', '--preprocess_dir', action='store',
    #    dest='preprocess_dir', default='preprocessed_files',
    #    help='directory used to store preprocessing output files (default: preprocessed_files)')
    sp_filter_p.add_argument('--align_qual_lim', action='store', dest='qlm_inp', type=float, default=0.75,
        help='Lowest limit for alignment quality (default: 0.75)')
    sp_filter_p.add_argument('--min_clipped_reads', action='store', dest='tcr_inp', type=int, default=2, 
        help='Minimum number of clipped reads (default: 2)')
    sp_filter_p.add_argument('--min_clipped_and_dischord_reads', action='store', dest='trd_inp', type=int, default=5, 
        help='Minimum total [clipped+discordant] reads (default: 5)')
    sp_filter_p.add_argument('--read_percent', action='store', dest='rp_inp', type=float, default=10.0, 
        help='read percent value (default: 10.0)')
#    sp_filter_p.add_argument('--read_length', action='store', dest='rdl_inp', type=int, default=100, 
#        help='Average read length (default: 100)')
#    sp_filter_p.add_argument('--insert_size_est', action='store', dest='isz_inp', type=int, default=369, 
#        help='insert Size estimate (default: 369)')
    sp_filter_p.add_argument('--log_file', action='store',
        dest='log_file', default='filter_p.log',
        help='run log file (default: filter_p.log)')
    sp_filter_p._action_groups.reverse()

    filter_p_ceu_desc=str("alternative filter for output from analyze step")
    sp_filter_p_ceu = subparsers.add_parser('filter_p_ceu', 
        help=filter_p_ceu_desc, description=filter_p_ceu_desc)
    sp_filter_p_ceu_required = sp_filter_p_ceu.add_argument_group('required arguments')
    sp_filter_p_ceu_required.add_argument('-i', '--input_file', action='store', dest='ofa_inp', required=True,
        help='use the output file from analyze section as input')
    sp_filter_p_ceu_required.add_argument('-b', '--bed', action='store', dest='fofn_bed', required=True, 
        help='File containg a list of files to existing repeat elements. List the full path for each file. See example in example_dir')
    #sp_filter_p_ceu.add_argument('-p', '--preprocess_dir', action='store',
    #    dest='preprocess_dir', default='preprocessed_files',
    #    help='directory used to store preprocessing output files (default: preprocessed_files)')
    sp_filter_p_ceu.add_argument('--align_qual_lim', action='store', dest='qlm_inp', type=float, default=0.75,
        help='Lowest limit for alignment quality (default: 0.75)')
    sp_filter_p_ceu.add_argument('--min_clipped_reads', action='store', dest='tcr_inp', type=int, default=2, 
        help='Minimum number of clipped reads (default: 2)')
    sp_filter_p_ceu.add_argument('--min_clipped_and_dischord_reads', action='store', dest='trd_inp', type=int, default=5, 
        help='Minimum total [clipped+discordant] reads (default: 5)')
    sp_filter_p_ceu.add_argument('--read_percent', action='store', dest='rp_inp', type=float, default=10.0, 
        help='read percent value (default: 10.0)')
#    sp_filter_p_ceu.add_argument('--read_length', action='store', dest='rdl_inp', type=int, default=100, 
#        help='Average read length (default: 100)')
#    sp_filter_p_ceu.add_argument('--insert_size_est', action='store', dest='isz_inp', type=int, default=369, 
#        help='insert Size estimate (default: 369)')
    sp_filter_p_ceu.add_argument('--log_file', action='store',
        dest='log_file', default='filter_p_ceu.log',
        help='run log file (default: filter_p_ceu.log)')
    sp_filter_p_ceu._action_groups.reverse()

    args = parser.parse_args()
    funct = FUNCTION_MAP[args.command]
    funct(args)

if __name__ == '__main__':
    import argparse
    from TEdetective.preprocess import exec_preprocess
    from TEdetective.discover import exec_discover
    from TEdetective.nadiscover import exec_nadiscover
    from TEdetective.cluster2d import exec_cluster2d
    from TEdetective.filters import exec_filter, exec_filter_p, exec_filter_p_ceu
    from TEdetective.analyze import exec_analyze
    from TEdetective.polymorph_screen import exec_polymorph, polymorph_setup_arg_parser
    main() 

