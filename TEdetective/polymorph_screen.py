
import argparse 
from TEdetective.io_functions import eprint
import TEdetective.polymorph_screen_functions as pm_fun
    
def main():
    parser = argparse.ArgumentParser()
    arg_parser = polymorph_setup_arg_parser(parser)
    args = arg_parser.parse_args()
    exec_polymorph(args)
    
def exec_polymorph(args):
    if args.verbose > 0:
        eprint("reading input file: " + args.input_file) 
    if args.verbose > 0:
        eprint("filtering input file.") 
    filter_results = pm_fun.filter_input_file(args.input_file, args.filter, args.qual_thresh)
    filter_cnt = 0   
    filter_file_names = dict()
    if args.screen_file_list != 'None':
        for file_name in args.screen_file_list.split(','):
            filter_cnt += 1
            filter_file_names[filter_cnt] = file_name
            if args.verbose > 0:
                eprint("reading filter file: " + file_name) 
            filter_results = pm_fun.add_filter_data(filter_results, file_name, filter_cnt, args.pm_qual_thresh)
    if args.verbose > 0:
        eprint("writing results")
    total_initial_pass,total_pass,total_not_found,total_initial_predictions,total_filtered,results = pm_fun.calc_filter_results(args.input_file, filter_cnt, filter_results)
    pm_fun.write_stats(total_initial_pass,total_pass,total_not_found,total_initial_predictions,total_filtered,filter_file_names,args.stats_file)
    pm_fun.write_results(results,filter_cnt,filter_file_names,args.output_file)
    return


def polymorph_setup_arg_parser(parser):
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument('-i', '--input_file', action='store', 
        dest='input_file', required=True, 
        help="prediction list to be filtered (output from analyze module)")
    parser.add_argument('-s', '--screen_file(s)', action='store', 
        dest='screen_file_list', default='None', 
        help="comma separated list of prediction file(s) for filtering, set to None to skip filter (default: None)")
    parser.add_argument('-o', '--output_file', action='store', 
        dest='output_file', default="filter_output.txt",
        help="ouput file (default: filter_output.txt)")
    parser.add_argument('--stats_file', action='store', 
        dest='stats_file', default="filter_stats.txt",
        help="filter statistics output file (default: filter_stats.txt)")
    parser.add_argument('-f', '--filter', action='store', 
        dest='filter', default="ceu",
        help="initial filtering criteria, see README.md for more info (default: ceu)")
    parser.add_argument('--qual_thresh', action='store', 
        dest='qual_thresh', default=0.85, type=float,
        help="initial filter quality threshold for clipped and discordant read alignments. (default: 0.85)")
    parser.add_argument('--pm_qual_thresh', action='store', 
        dest='pm_qual_thresh', default=0.75, type=float,
        help="initial filter quality threshold for clipped and discordant read alignments. (default: 0.75)")
    parser.add_argument('-v', '--verbose', action='store', 
        dest='verbose', default=1, type=int,
        help="verbose level, set to 0 for quiet. (default: 1)")
    parser._action_groups.reverse()
    return parser


if __name__ == '__main__':
    import argparse 
    from io_functions import eprint
    import polymorph_screen_functions as pm_fun
    main() 