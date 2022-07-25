import argparse 

from TEdetective.io_functions import eprint
import TEdetective.filter_functions as flt_fun
    
def main():
    parser = argparse.ArgumentParser()
    arg_parser = filter_setup_arg_parser(parser)
    args = arg_parser.parse_args()
    exec_filter(args)
    
def exec_filter(args):
    if args.verbose > 0:
        eprint("reading input file: " + args.input_file) 
    if args.verbose > 0:
        eprint("filtering input file.") 
    filter_results = flt_fun.filter_input_file(args.input_file, args.ini_filter, args.qual_thresh, args.te_type,args.no_polyA_info)
    filter_cnt = 0   
    filter_file_names = dict()
    if args.screen_file_list != 'None':
        for file_name in args.screen_file_list.split(','):
            filter_cnt += 1
            filter_file_names[filter_cnt] = file_name
            if args.verbose > 0:
                eprint("reading filter file: " + file_name) 
            filter_results = flt_fun.add_filter_data(filter_results, file_name, filter_cnt, args.pm_qual_thresh, args.filter, args.te_type,args.no_polyA_info)
    if args.results_screen_file_list != 'None':
        for file_name in args.results_screen_file_list.split(','):
            filter_cnt += 1
            filter_file_names[filter_cnt] = file_name
            if args.verbose > 0:
                eprint("reading filter file: " + file_name) 
            filter_results = flt_fun.add_filter_other_results(filter_results, file_name, args.insert_size, args.read_length)
    if args.bed_screen:
        if args.verbose > 0:
            eprint("filtering all " + args.te_type + " with bed file: " + args.bed_screen) 
        filter_results = flt_fun.add_filter_existing_data(filter_results, args.bed_screen, args.input_file, args.te_type)
        filter_cnt += 1
        filter_file_names[filter_cnt] = "in_existing_" + args.te_type
    if args.filter_alt_chrom:
        if args.verbose > 0:
            eprint("filtering alt chromosomes") 
        filter_results,total_alt_chrom = flt_fun.add_filter_alt_chrom(filter_results)
        filter_cnt += 1
        filter_file_names[filter_cnt] = "alt_chrom"
    if args.verbose > 0:
        eprint("writing results")
    total_initial_pass,total_pass,total_not_found,total_initial_predictions,total_filtered,results = flt_fun.calc_filter_results(args.input_file, filter_cnt, filter_results)
    flt_fun.write_stats(total_initial_pass,total_pass,total_not_found,total_initial_predictions,total_filtered,filter_file_names,args.stats_file)
    flt_fun.write_results_mask(results,filter_cnt,filter_file_names,args.output_file + ".mask")
    total_pass_all_filters = flt_fun.write_results(filter_results,filter_cnt,args.input_file,args.output_file)
    if args.verbose > 0:
        eprint("total_passing_all_filters: " + str(total_pass_all_filters))
    return


def filter_setup_arg_parser(parser):
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument('-i', '--input_file', action='store', 
        dest='input_file', required=True, 
        help="prediction list to be filtered (output from analyze module)")
    parser.add_argument('-s', '--screen_files', action='store', 
        dest='screen_file_list', default='None', 
        help="comma separated list of prediction file(s) for filtering, set to None to skip filter (default: None)")
    parser.add_argument('--results_screen_files', action='store', 
        dest='results_screen_file_list', default='None', 
        help="comma separated list of results file(s) for filtering, filters TEs within insert_size - read_length. Set to None to skip filter (default: None)")
    parser.add_argument('-b', '--bed_screen', action='store', 
        dest='bed_screen', default=False, 
        help="bed formated file of TE positions, file will be screened for TE type (default: Off)")
    parser.add_argument('-t', '--te_type', action='store', 
        dest='te_type', default='LINE', 
        help="TE type used to screen bed file (default: LINE)")
    parser.add_argument('--filter_alt_chrom', action='store', 
        dest='filter_alt_chrom', default=True, 
        help="filter alternative chromosomoes such as random, unknown, fix, hap (default: True)")
    parser.add_argument('-o', '--output_file', action='store', 
        dest='output_file', default="filter_output.txt",
        help="ouput file (default: filter_output.txt)")
    parser.add_argument('--stats_file', action='store', 
        dest='stats_file', default="filter_stats.txt",
        help="filter statistics output file (default: filter_stats.txt)")
    parser.add_argument('-f', '--filter', action='store', 
        dest='filter', default="ceu",
        choices = ['ceu','custom','normal'],
        help="polymorphic subtraction filtering criteria, see README.md for more info (default: ceu)")
    parser.add_argument('--ini_filter', action='store', 
        dest='ini_filter', default="ceu",
        choices = ['ceu','stringent','custom','normal','ceu_p'],
        help="initial filtering criteria, see README.md for more info (default: ceu)")
    parser.add_argument('--qual_thresh', action='store', 
        dest='qual_thresh', default=0.85, type=float,
        help="initial filter quality threshold for clipped and discordant read alignments. (default: 0.85)")
    parser.add_argument('--pm_qual_thresh', action='store', 
        dest='pm_qual_thresh', default=0.75, type=float,
        help="initial filter quality threshold for clipped and discordant read alignments. (default: 0.75)")
    parser.add_argument('--insert_size_est', action='store', dest='insert_size', type=int, default=340, 
        help='insert size estimate (default: 340)')
    parser.add_argument('--read_length', action='store', dest='read_length', type=int, default=100, 
        help='read length (default: 100)')
    parser.add_argument('--no_polyA_info', action='store', dest='no_polyA_info', default=False, 
        help='for backward compatibility with some old results files. (default: False)')
    parser.add_argument('-v', '--verbose', action='store', 
        dest='verbose', default=1, type=int,
        help="verbose level, set to 0 for quiet. (default: 1)")
    parser._action_groups.reverse()
    return parser


if __name__ == '__main__':
    import argparse 
    from io_functions import eprint
    import filter_functions as flt_fun
    main() 