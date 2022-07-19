################################################################
# Python based tool to identify novel transposable element insertions from WGS data.
# Contributors: Manoj Kumar Singh (manoj@wustl.edu),
#        and John Edwards (jredwards@wustl.edu)
#
# External dependencies: censor, NCBI blast (provided with package)
#
################################################################

import argparse

from TEdetective.preprocess import exec_preprocess, preprocess_setup_arg_parser
from TEdetective.discover import exec_discover, discover_setup_arg_parser
from TEdetective.nadiscover import exec_nadiscover, nadiscover_setup_arg_parser
from TEdetective.cluster2d import exec_cluster2d, cluster2d_setup_arg_parser
from TEdetective.old_filters import exec_old_filter, exec_old_filter_p, exec_old_filter_p_ceu, old_filter_setup_arg_parser
from TEdetective.analyze import exec_analyze, analyze_setup_arg_parser
from TEdetective.filter import exec_filter, filter_setup_arg_parser

def main():
    FUNCTION_MAP = {
            'preprocess' : exec_preprocess, 
            'discover' : exec_discover,
            'nadiscover' : exec_nadiscover, 
            'analyze' : exec_analyze,
            'cluster2d' : exec_cluster2d,
            'old_filter' : exec_old_filter,
            'old_filter_p' : exec_old_filter_p,
            'old_filter_p_ceu' : exec_old_filter_p_ceu,
            'filter' : exec_filter
            }

    parser = argparse.ArgumentParser()
    parser = setup_parsers(parser)   
    args = parser.parse_args()
    funct = FUNCTION_MAP[args.command]
    funct(args)

def setup_parsers(parser):
    subparsers = parser.add_subparsers(dest='command')
    subparsers.required = True
    
    preprocess_desc=str("Processes the input files (indexed BAM file and indexed " +
        "fasta file), extracts discordant and clipped reads, as well as  " +
        "creates other files needed in subsequent steps. Outputs all files " +
        "to directory specified by --preprocess_dir")
    sp_preprocess = subparsers.add_parser('preprocess',
        description=preprocess_desc, help=preprocess_desc)
    sp_preprocess = preprocess_setup_arg_parser(sp_preprocess)
    
    discover_desc=str("Uses output from preprocessing step and makes an initial list of candidate insertions")
    sp_discover = subparsers.add_parser('discover',
        help=discover_desc, description=discover_desc)
    sp_discover = discover_setup_arg_parser(sp_discover)
    
    nadiscover_desc=str("Performs nonalignment part of the discovery step. " +
          "Module adds poly A/T information into predictions made " +
          "by discovery step. This module performs initial searches " +
          "as well, but without using BWA aligner for clipped and " +
          "discordant read alignment to TE reference sequence. " +
          "Instead, a bed file of masked regions is provided as " +
          "input, and alignment information from input BAM file is used.")
    sp_nadiscover = subparsers.add_parser('nadiscover', 
        description=nadiscover_desc,help=nadiscover_desc)
    sp_nadiscover = nadiscover_setup_arg_parser(sp_nadiscover)

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
    sp_analyze = analyze_setup_arg_parser(sp_analyze)

    filter_desc=str("Filter module to filter predictions from nadiscover " +
        "based on read support, existing TEs, and polymorphic subtraction")
    sp_filter = subparsers.add_parser('filter',
        help=filter_desc, description=filter_desc)
    sp_filter = filter_setup_arg_parser(sp_filter)
    
    cluster2d_desc=str("Optional module to change the discordant read clustering density for " +
          "initial prediction without realigning everything. For example, if " +
          "--discord_cluster_dens was set to 10 for initial discovery step " +
          "and user want to see predictions with --discord_cluster_dens = 5. " +
          "Uses intermediate files from discover section and generates new prediction file.")
    sp_cluster2d = subparsers.add_parser('cluster2d', 
        help=cluster2d_desc, description=cluster2d_desc)
    sp_cluster2d = cluster2d_setup_arg_parser(sp_cluster2d)

    old_filter_desc=str("OLD. Filter output from analyze step." +
            "# Filteration step code looks like this\n" +
            "if total_clipped_rd >= tcr or ( (total_clipped_rd >= mtcr ) " +
            "and ( (total_clipped_rd_wpat+total_discord_rd) >= trd ) ):\n" +
            "filter_result = 'PASS'\n" +
            "elif total_discord_rd >= odrd: \n" +
            "filter_result = 'PASS_D' # This flag says passed based on only discordant reads.")
    sp_old_filter = subparsers.add_parser('old_filter', 
        help=old_filter_desc, description=old_filter_desc)
    sp_old_filter = old_filter_setup_arg_parser(sp_old_filter)
    
    old_filter_p_desc=str("OLD. alternative filter for output from analyze step")
    sp_old_filter_p = subparsers.add_parser('old_filter_p', 
        help=old_filter_p_desc, description=old_filter_p_desc)
    sp_old_filter_p = old_filter_setup_arg_parser(sp_old_filter_p)

    filter_p_ceu_desc=str("OLD. alternative filter for output from analyze step")
    sp_old_filter_p_ceu = subparsers.add_parser('old_filter_p_ceu', 
        help=filter_p_ceu_desc, description=filter_p_ceu_desc)
    sp_old_filter_p_ceu = old_filter_setup_arg_parser(sp_old_filter_p_ceu)
    
    return parser

if __name__ == '__main__':
    import argparse
    from TEdetective.preprocess import exec_preprocess, preprocess_setup_arg_parser
    from TEdetective.discover import exec_discover, discover_setup_arg_parser
    from TEdetective.nadiscover import exec_nadiscover, nadiscover_setup_arg_parser
    from TEdetective.cluster2d import exec_cluster2d, cluster2d_setup_arg_parser
    from TEdetective.old_filters import exec_old_filter, exec_old_filter_p, exec_old_filter_p_ceu, old_filter_setup_arg_parser
    from TEdetective.analyze import exec_analyze,analyze_setup_arg_parser
    from TEdetective.filter import exec_filter, filter_setup_arg_parser
    main() 

